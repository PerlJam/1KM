#! env python
"""Update small molecule csv 
     python small_molecule_update.py  ../tests/small_molecule_initial.csv out.csv

   Author: Samuel A. Inverso, Wyss Institute, 2014-06-11
"""

# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch
import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from Bio import SeqIO 
from Bio import Entrez
import pandas as pd
import numpy as np
import csv

from collections import deque

from copy import copy

from pprint import pprint

from progressbar import ProgressBar

import logging

import pubchem_fetch
import chebi_fetch
import hmdb_fetch
import kegg_fetch
import drugbank_fetch
import convert_id_1km
import fix_value

logger = logging.getLogger('1km')
logger.setLevel(logging.DEBUG)

class SmallMolecule:
    """Class represents a small molecule table row 
        updates missing elements from public databases
    """

    SEPARATOR = ";" # The separator for multiple values in a cell
    # The ID Columns, numbers are priority
    ID_COLUMNS = {'PubChem_CID':1,  'ChEBI_ID':2, 'KEGG_ID':3, 'BiGG_ID':4, 
        'HMDB_ID':5, 'DrugBank_ID':6}

    # THe columns, true indicates it's an ID column, false otherwise
    COLUMN_NAMES = {'Name':False, 'IUPAC_Name':False,  'AlternativeNames':False, 
        'Metabolic_Network_ID':False, 
        'PubChem_CID':True, 'ChEBI_ID':True, 'KEGG_ID':True, 'BiGG_ID':True, 
        'HMDB_ID':True,
        'DrugBank_ID':True, 'InChi_Parent':False, 'InChi_Key_Parent':False,
        'SMILES_Parent':False,
        'Protein_Associations':False,  'Molecular_Mass_Canonical':False,    
        'Molecular_Formula_Canonical':False, 'Formal_Charge':False, 
        'Physiological_Charge':False }


    def __init__(self,row):
        """Takes a dataframe row representing a SmallMolecule"""
        #self.row_to_class(row)
        self._row = row
        self._cache = {} 


    @staticmethod
    def compare_db_ids(id1,id2):
        """comparator for sorted"""
        return SmallMolecule.ID_COLUMNS[id1] - SmallMolecule.ID_COLUMNS[id2] 

    
    def resolve_missing_ids(self,attempt=0):
        """attempts to find missing ids

            returns True if we added any missing ids, false otherwise
        """
        result = False
        MAX_ATTEMPT = 1
        ids_seen = [] # id's we've already gone to get data from

        # get the missing ids and ids_to_check for the missing ones
        #   note ids_to_check is a dequeu and we popleft 
        (missing_ids, ids_to_check) = self.get_ids_missing_and_to_check(ids_seen)
        num_missing = len(missing_ids)
        # while are still missing ids and there are ids we can check
        while missing_ids and ids_to_check:

            # Get the id to check 
            id_to_check = ids_to_check.popleft()

            ids_seen.append(id_to_check)

            # Update the ids, different attempts try from in different ways
            logger.debug('Resolving Missing ID: Attempt %d, %s %s', attempt, id_to_check, "; ".join(missing_ids))

            if attempt == 0:
                # use the external database converion tool, though should
                # be fast ish
                self.update_missing_ids_from_1km(id_to_check,missing_ids)
            elif attempt == 1:
                # Actually pull the full record from the database and
                # try to resolve the ids - this is slower
                self.update_from_db(id_to_check,missing_ids)
            else:
                break

            # update the missing_ids and to check
            (missing_ids, ids_to_check) = self.get_ids_missing_and_to_check(ids_seen)
            if num_missing > len(missing_ids):
                result = True # we updated some ids
                num_missing = len(missing_ids)


        # if there are still missing ids try another attempt
        if missing_ids and attempt < MAX_ATTEMPT:
            # we may have already update missing ids, so we want
            # to remain true if we are true
            result = self.resolve_missing_ids(attempt+1) or result

        return result

    
    def update_missing_ids_from_1km(self,id_to_check,missing_ids):
        """Update missing ids from the convert_id_1km which is a pass through to fiehnlab database
            id_to_check is the database column us to resolve the missing_ids database columns
            e.g. PubChem_CID  [ChEBI_ID]
        """
        ## loop through the missing id columns and try to convert them
        for to_db in missing_ids:
            id_value = str(self._row[id_to_check].iloc[0])
            converted = convert_id_1km.convert_id(id_to_check,to_db,id_value)
            # if we actually got something back then set it
            if converted:
                self._row[to_db] = converted[0]


    def get_ids_missing_and_to_check(self,ids_seen):
        """Helper function for resolve_missing_ids
            gets the missing ids from the row, and returns the missing and
            the ones we still need to check, i.e. those present in the row
            that are not in ids_seen

            Returns:
                missing_ids: list the missing ids
                ids_to_check: deque the ids that still need to be checked
        """
        missing_cols = self.missing_columns()

        missing_ids = []
        present_ids = []
        for col in missing_cols:
            # test if it's an id
            if SmallMolecule.COLUMN_NAMES[col]:
                if missing_cols[col].bool():
                    missing_ids.append(col)
                else:
                    present_ids.append(col)

        # only check the ones we have not seen yet
        ids_to_check = list(set(present_ids) - set(ids_seen))
    
        # Make sure the ids are in priority order 
        ids_to_check = sorted(ids_to_check,cmp=SmallMolecule.compare_db_ids)

        return (missing_ids, deque(ids_to_check))


    def fix_values(self):
        """Fixes the values so they correspond to our format
             e.g. InChI= prepended to inchi
        """
        fix_value.fix_all(self._row)


    def update(self):
        """Updates missing pieces"""
        # Test if there are missing values
        #   if there are missing values continue
        #   else we are done

        # first resolve any missing ids
        updated_any =  self.resolve_missing_ids()
        # determine if anything is still missing,
        missing_cols = self.missing_columns()
        # we do not want to recheck ids so remove them
        for id_col in SmallMolecule.ID_COLUMNS:
            if missing_cols[id_col].values.any():
                missing_cols[id_col] = False

        # get the missing column names
        missing_col_names =  missing_cols.columns[missing_cols.values[0]]

        # we have at least one missing
        if missing_cols.values.any():

            # now in database priority order, resolve the 
            # the missing columns

            id_columns = copy(SmallMolecule.ID_COLUMNS)
            id_columns = sorted(id_columns,cmp=SmallMolecule.compare_db_ids)

            while id_columns:
                db_id = id_columns.pop()
                self.update_from_db(db_id,missing_col_names)
        
        # check if we updated any columns
        if len(missing_cols) > len(self.missing_columns()):
            updated_any = True

        # if we've updated any then we want to make sure they are correct values
        # so fix them
        self.fix_values()

        return updated_any


    def update_from_db(self,db_id,cols):
        """Calls the correct function to update from the given database
        """
        if db_id == 'PubChem_CID':
            self.update_from_pubchem(cols)
        elif db_id == 'ChEBI_ID':
            self.update_from_chebi(cols)
        elif db_id == 'KEGG_ID':
            self.update_from_kegg(cols)
        elif db_id == 'HMDB_ID':
            self.update_from_hmdb(cols)
        elif db_id == 'DrugBank_ID':
            self.update_from_drugbank(cols)


    def get_record(self,db_name,db_id):
        """Checks the cache for the record otherwise goes out to the 
           databse to get it
           # note we only cache one record per database.
           """
        result = None
        if db_name in self._cache:
            result = self._cache[db_name]
            # TODO Test t osee if the cached one has the same id
        else:
            record_id = self._row.iloc[0].get(db_id,None)
            # Test that the record_id exists, that it's not None (implicit),
            # and that if it's a numpy number that it is not an nan
            if record_id and not (isinstance(record_id,float) and np.isnan(record_id)):
                # if it is a float we need to make it an int
                if 'PUBCHEM' == db_name:
                    result = pubchem_fetch.from_id(record_id)
                elif 'CHEBI' == db_name:
                    result = chebi_fetch.from_id(record_id)
                elif 'KEGG' == db_name:
                    result = kegg_fetch.from_id(record_id)
                elif 'HMDB' == db_name:
                    result = hmdb_fetch.from_id(record_id)
                elif 'DRUGBANK' == db_name:
                    result = drugbank_fetch.from_id(record_id)
                else:
                    raise Exception("Unknown database name: " + db_name)
                self._cache[db_name] = result

        return result


    def update_from_pubchem(self,cols):
        """Update row from pubchem"""
        db = self.get_record('PUBCHEM','PubChem_CID')
        if db is None:
            # we were unable to get the pubchem database
            # we're going to return
            logger.warning('Unable to get PubChem database.')
            return

        for col in cols:
            if 'Name' == col:
                # use the first alternte name
                self._row[col] = db.synonyms[0]
            elif 'IUPAC_Name' == col:
                self._row[col] = db.iupac_name
            elif 'AlternativeNames' == col:
                self._row[col] = SmallMolecule.SEPARATOR.join(db.synonyms)
            elif 'Metabolic_Network_ID' == col:
                # unknown
                pass
            elif  'PubChem_CID' == col: 
                # we already have it if we got this far
                pass
            elif 'ChEBI_ID' == col:
                # we get this from the synonyms
                for syn in db.synonyms:
                    if syn.startswith('CHEBI:'):
                        self._row[col]  = syn
                        break
            elif 'KEGG_ID' == col:
                # unknown
                pass
            elif 'BiGG_ID' == col:
                # unknown
                pass
            elif 'HMDB_ID' == col:
                # unknown
                pass
            elif 'DrugBank_ID' == col:
                # unknown
                pass
            elif 'InChi_Parent' == col:
                self._row[col] = db.inchi
            elif 'InChi_Key_Parent' == col:
                self._row[col] = db.inchikey
            elif 'SMILES_Parent' == col:
                self._row[col] = db.canonical_smiles
            elif 'Protein_Associations' == col:
                # unknown
                pass
            elif 'Molecular_Mass_Canonical' == col:
                self._row[col] = db.molecular_weight
            elif 'Molecular_Formula_Canonical' == col:
                self._row[col] = db.molecular_formula
            elif 'Formal_Charge' == col:
                self._row[col] = db.charge
            elif 'Physiological_Charge' == col:
                pass
            else:
                # we ignore oneswe can't deal with
                pass


    def update_from_chebi(self,cols):
        """Update row from chebi"""
        db = self.get_record('CHEBI', 'ChEBI_ID')

        if db is None:
            # we were unable to get the  database
            # we're going to return
            logger.warning('Unable to get ChEBI database.')
            return

        # it's faster if we pre pull the database link, otherwise 
        #   we have to loop through it multiple times
        kegg = None
        drugbank = None
        hmdb = None
        if 'DatabaseLinks' in db:
            for dbl in db.DatabaseLinks:
                if "KEGG COMPOUND accession"  == dbl['type']:
                    kegg = dbl['data']
                elif "DrugBank accession" == dbl['type']:
                    drugbank = dbl['data']
                elif "HMDB accession" == dbl['type']:
                    hmdb = dbl['data']

        for col in cols:
            if 'Name' == col and 'chebiAsciiName' in db:
                self._row[col] = db.chebiAsciiName
                pass
            elif 'IUPAC_Name' == col and 'IupacNames' in db:
                self._row[col] = SmallMolecule.SEPARATOR.join(
                    [x['data'] for x in db.IupacNames] )
            elif 'AlternativeNames' == col and 'Synonyms' in db:
                # don't tknow this mapping
                self._row[col] =  SmallMolecule.SEPARATOR.join(
                        [x['data'] for x in db.Synonyms]  )
            elif 'Metabolic_Network_ID' == col:
                pass
            elif  'PubChem_CID' == col: 
                pass
            elif  'ChEBI_ID' == col:
                # this is chebi
                pass
            elif 'KEGG_ID' == col and kegg is not None:
                self._row[col] = kegg
            elif 'BiGG_ID' == col:
                pass
            elif 'HMDB_ID' == col and hmdb is not None:
                self._row[col] = hmdb
            elif 'DrugBank_ID' == col and drugbank is not None:
                self._row[col] = drugbank
            elif 'InChi_Parent' == col and 'inchi' in db:
                self._row[col] = db.inchi
            elif 'InChi_Key_Parent' == col and 'inchiKey' in db:
                self._row[col] = db.inchiKey
            elif 'SMILES_Parent' == col and 'smiles' in db:
                self._row[col] = db.smiles
            elif 'Protein_Associations' == col:
                pass
            elif 'Molecular_Mass_Canonical' == col and 'mass' in db:
                self._row[col] = db.mass
            elif 'Molecular_Formula_Canonical' == col and 'Formulae' in db:
                self._row[col] = SmallMolecule.SEPARATOR.join(
                    [x['data'] for x in db.Formulae])
            elif 'Formal_Charge' == col and 'charge' in db:
                self._row[col] = db.charge
            elif 'Physiological_Charge' == col:
                pass
            else:
                # we ignore oneswe can't deal with
                pass


    def update_from_kegg(self,cols):
        """update row from kegg """
        db = self.get_record('KEGG', 'KEGG_ID')

        db_d = kegg_fetch.get_drug_record(db) # the drug record that corresponds to this cid 

        if db is None:
            # we were unable to get the  database
            # we're going to return
            logger.warning('Unable to get KEGG database.')
            return

        # names comes as a list
        names = db['name']
        for col in cols:
            if 'Name' == col:
                if names:
                    self._row[col] = names[0] 
                pass
            elif 'IUPAC_Name' == col:
                pass
            elif 'AlternativeNames' == col:
                if len(names) > 1:
                    self._row[col] = SmallMolecule.SEPARATOR.join(names[1:]) 
            elif 'Metabolic_Network_ID' == col:
                pass
            elif  'PubChem_CID' == col: 
                if 'PubChem' in db['dblinks']:
                    sid = db['dblinks']['PubChem']
                    # kegg returns an sid, conver to cid
                    cid = pubchem_fetch.sid_to_cid(sid)
                    if cid:
                        self._row[col] = cid
            elif  'ChEBI_ID' == col:
                if 'ChEBI' in db['dblinks']:
                    self._row[col] =  db['dblinks']['ChEBI']
                pass
            elif 'KEGG_ID' == col:
                pass
            elif 'BiGG_ID' == col:
                pass
            elif 'HMDB_ID' == col:
                pass
            elif 'DrugBank_ID' == col:
                if db_d and 'DrugBank' in db_d['dblinks']:
                    self._row[col] = db_d['dblinks']['DrugBank']
            elif 'InChi_Parent' == col:
                pass
            elif 'InChi_Key_Parent' == col:
                pass
            elif 'SMILES_Parent' == col:
                pass
            elif 'Protein_Associations' == col:
                pass
            elif 'Molecular_Mass_Canonical' == col and 'exact_mass' in db:
                self._row[col] = db['exact_mass'] 
            elif 'Molecular_Formula_Canonical' == col:
                self._row[col] = db['formula']
            elif 'Formal_Charge' == col:
                pass
            elif 'Physiological_Charge' == col:
                pass
            else:
                # we ignore oneswe can't deal with
                pass


    def _assign_if_not_empty(self,col,val):
        """ This is mostly just for update_from_hmdb and update_from_drugbank
            the value is always a list 
            if it's one element we just unwrap it. if it's 
            multiple elements then we combine them with the SEPARATOR 
        """
        # check if it's not empty or none
        if val:
            if len(val) == 1:
                self._row[col] = val[0] # just one element getit 
            else:
                self._row[col] = SmallMolecule.SEPARATOR.join( val )


    def update_from_hmdb(self,cols):
        """ Update row from hmdb
        """
        db = self.get_record('HMDB', 'HMDB_ID')

        if db is None:
            # we were unable to get the  database
            # we're going to return
            logger.warning('Unable to get HMDB database.')
            return

        for col in cols:
            if 'Name' == col:
                self._assign_if_not_empty(col,db['name'])
            elif 'IUPAC_Name' == col:
                self._assign_if_not_empty(col,db['iupac_name'])
            elif 'AlternativeNames' == col:
                # don't tknow this mapping
                pass
            elif 'Metabolic_Network_ID' == col:
                pass
            elif  'PubChem_CID' == col: 
                self._assign_if_not_empty(col,db['pubchem_compound_id'])
            elif  'ChEBI_ID' == col:
                self._assign_if_not_empty(col,db['chebi_id'])
            elif 'KEGG_ID' == col:
                self._assign_if_not_empty(col,db['kegg_id'])
            elif 'BiGG_ID' == col:
                self._assign_if_not_empty(col,db['bigg_id'])
            elif 'HMDB_ID' == col:
                pass
            elif 'DrugBank_ID' == col:
                self._assign_if_not_empty(col,db['drugbank'])
            elif 'InChi_Parent' == col:
                self._assign_if_not_empty(col,db['inchi'])
            elif 'InChi_Key_Parent' == col:
                self._assign_if_not_empty(col,db['inchikey'])
            elif 'SMILES_Parent' == col:
                self._assign_if_not_empty(col,db['smiles'])
            elif 'Protein_Associations' == col:
                self._assign_if_not_empty(col,db['protein_associations//protein//uniprot_id'])
            elif 'Molecular_Mass_Canonical' == col:
                self._assign_if_not_empty(col,db['average_molecular_weight'])
            elif 'Molecular_Formula_Canonical' == col:
                self._assign_if_not_empty(col,db['chemical_formula'])
            elif 'Formal_Charge' == col:
                self._assign_if_not_empty(col,db[
                    "predicted_properties//property/.[kind='formal_charge']//value"])
            elif 'Physiological_Charge' == col:
                self._assign_if_not_empty(col,db[
                    "predicted_properties//property/.[kind='physiological_charge']//value"])      
            else:
                # we ignore oneswe can't deal with
                pass



    def update_from_drugbank(self,cols):
        """ Update row from hmdb
        """
        db = self.get_record('DRUGBANK', 'DrugBank_ID')

        if db is None:
            # we were unable to get the  database
            # we're going to return
            logger.warning('Unable to get DrugBank database.')
            return

        for col in cols:
            if 'Name' == col:
                self._assign_if_not_empty(col,db.name())
            elif 'IUPAC_Name' == col:
                self._assign_if_not_empty(col,db.iupac_name())
            elif 'AlternativeNames' == col:
                # don't tknow this mapping
                pass
            elif 'Metabolic_Network_ID' == col:
                pass
            elif  'PubChem_CID' == col: 
                self._assign_if_not_empty(col,db.pubchem_compound_id())
            elif  'ChEBI_ID' == col:
                self._assign_if_not_empty(col,db.chebi())
            elif 'KEGG_ID' == col:
                self._assign_if_not_empty(col,db.kegg_compound())
            elif 'BiGG_ID' == col:
                pass
            elif 'HMDB_ID' == col:
                pass
            elif 'DrugBank_ID' == col:
                pass
            elif 'InChi_Parent' == col:
                self._assign_if_not_empty(col,db.inchi())
            elif 'InChi_Key_Parent' == col:
                self._assign_if_not_empty(col,db.inchi_key())
            elif 'SMILES_Parent' == col:
                self._assign_if_not_empty(col,db.smiles())
            elif 'Protein_Associations' == col:
                pass
            elif 'Molecular_Mass_Canonical' == col:
                self._assign_if_not_empty(col,db.molecular_weight())
            elif 'Molecular_Formula_Canonical' == col:
                self._assign_if_not_empty(col,db.formula())
            elif 'Formal_Charge' == col:
                pass
            elif 'Physiological_Charge' == col:
                self._assign_if_not_empty(col,db.physiological_charge())      
            else:
                # we ignore ones we can't deal with
                pass
            pass

    def missing_columns(self):
        """returns a pandas dataframe with columns marked True if missing
        """
        return self._row.isnull()


    def get_row(self):
        """returns a dataframe row"""
        return self._row


    def is_valid(self,row):
        """ check if this is a valid row"""
        # Do we know how to deal with all the columns?
        coldiff = set(row.columns) - set(self.COLUMN_NAMES)
        if coldiff:
            raise Exception("Row contains unknown columns: " + str(coldiff) )


    def flush_cache(self,db_id=None):
        """ removes all cached data, or just the one passed """
        if db_id:
            self._cache.pop(db_id)
        else:
            self._cache = {}


def update_row(df_row):
    """updates inplace, therefore if you passed a row from a dataframe
        it will be updated 
    """ 

    sm = SmallMolecule(df_row)
    sm.update()
    return sm.get_row()


def update(df,disp_progressbar=False):
    """Iterate over each row and update it, updates in place
    """
    rows = len(df)
    if disp_progressbar:
        pbar = ProgressBar(maxval=rows).start()

    for irow in xrange(rows):
        df_row = df.iloc[irow:(irow+1)]
        update_row(df_row).T
        if disp_progressbar:
            pbar.update(irow +1)
        df.iloc[irow] = df_row.iloc[0]

    return df


def read_csv(input_file):
    """Reads a csv file and returns a data frame
    """
    return pd.read_csv(input_file)


def write_csv(df, output_file,index=False):
    """writes the dataframe as a csv file"""
    df.to_csv(output_file,index=False)


def main(argv):
    """Main entry point of program.
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Converts a genebank file into a series of csv files for loading'
        ' into the database.\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('input', metavar='input.csv', action='store',
        help='Genbank file to process')
    parser.add_argument('output', metavar='output.csv', action='store',
        help='Output file name.')

    args = parser.parse_args()
    input_file = args.input
    output_file = args.output

    df = read_csv(input_file)

    df = update(df,True)
    #print df

    # write out the updated csv
    write_csv(df,output_file)

    
# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
