#! env python
"""Fetchs record from Kegg  
    Uses bioservices
    http://pythonhosted.org/bioservices/references.html?highlight=chebi#module-bioservices.chebi
   Author: Samuel A. Inverso, Wyss Institute, 2014-06-10

   Notes:
   KEGG's PubChem id can actually be the SID not the CID
"""


import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from pprint import pprint

import pandas as pd

from bioservices import *

def from_id(cid):
    """Pass through to KEGG get and parse
    """
    # make the parsers
    s = KeggParser()

    # get the data
    res = s.get(cid)

    # parse the data to dictionary
    db = s.parse(res)

    return db


def get_drug_id(db):
    """Finds the drug id from the db (db is result of from_id(cid))
        returns the kegg drug id
        e.g. for id=C00031  it returns D00009
    """
    drug_id  = None

    if db:
        # find the drug id
        if 'remark' in db:
            remarks  = db['remark']
            # it can be a string or a list
            # if it's a string, wrap it in a list so we can 
            # duplicate less code
            if isinstance(remarks,basestring):
                remarks = [remarks]
     
            for remark in remarks:
                if remark.startswith('Same as: '):
                    drug_id = remark[9:]
                    break

    return drug_id


def get_drug_record(db):
    """Finds the drug record from the db (db is result from_id(cid))
        conveniece function that wraps get_drug_id and from_id
    """
    result = None
    drug_id = get_drug_id(db)
    if drug_id:
        result = from_id(drug_id)

    return result


def main(argv):
    """Main entry point of program.
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Retrieves compounds from Kegg using a list of Kegg\n'
        'For example, \n'
        'python kegg_fetch.py "C00234',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('keggId', metavar='keggId', nargs='+',action='store',
        help='Kegg id to fetch.')
   # parser.add_argument('output_dir', metavar='output_dir', action='store',
   #     help='Output directory to place csvs into.')

    # parser.add_argument('-s', nargs=2, metavar=('identifier', 'type'), action='store', dest='search',
    #     help='Identifier and type to search for: e.g. -s Glucose name\n'
    #         '  Types supported: cid, name, smiles, sdf, inchi, inchikey, formula.\n')

    args = parser.parse_args()
    keggId = args.keggId
   
    s = KeggParser()
    for kId in keggId:
        res = s.get(kId)

        db = s.parse(res)
        pprint(db)

        print db['name'] # name and alternate names
        if 'ChEBI' in db['dblinks']:
            print db['dblinks']['ChEBI']
        print db['exact_mass']
        print db['formula']

        print 'KeggDrug id:', get_drug_id(db)
        db_d = get_drug_record(db)
        if 'DrugBank' in db_d['dblinks']:
            print "DrugBank", db_d['dblinks']['DrugBank']
        

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
