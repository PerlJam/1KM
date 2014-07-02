#! env python
"""Fetchs record from hmdb  
   Author: Samuel A. Inverso, Wyss Institute, 2014-06-18
"""

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from pprint import pprint

import pandas as pd

from bioservices import *

from hmdb_index import HMDB_INDEX, HMDB_DIRECTORY

import xml.etree.cElementTree as ET


class HMDB:
    def __init__(self,tree):
        self._tree = tree
        self._cache = {}


    def __getitem__(self,key):
        '''Allows access through brackets, e.g.
         foo['synonyms//synonym']
        '''
        result = None
        # check if it's tin the cache first
        if key in self._cache:
            result = self._cache[key]
        else:
            # it's not in the cache so retrieve it
            result = self._get_from_tree(key)
            # remove None values
            result = [x for x in result if x is not None]
            self._cache[key] = result

        return result


    def _get_from_tree(self,key):
        '''trys to get the key from the tree, cache's the result
        e.g. to get all the synonyms 
         key = 'synonyms//synonym' 
        '''
        result = []

        # There are some special cases 
        for node in self._tree.findall(key):
            result.append(node.text)

        return result


    def clear_cache(self):
        self._cache = {}


def from_id(cid):
    """Get's the HMDB class given an id
    """

    filename = HMDB_INDEX[cid]
    file_path = os.path.join(HMDB_DIRECTORY,filename)
    # load the xml
    tree = ET.ElementTree(file=file_path)


    return HMDB(tree)


def main(argv):
    """Main entry point of program.
    """
    parser = argparse.ArgumentParser(description=
        'Retrieves compounds from KeGG using either a list of ID e.g. HMDB00122\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('id', metavar='id', nargs='+',action='store',
        help='id(s) to fetch.')


    args = parser.parse_args()
    cids = args.id
  
    for cid in cids:

        db = from_id(cid)


        print db['name']
        print ";".join(  db['synonyms//synonym']   )

        print 'Iupac', db['iupac_name'] #
        # alternate name
        # metabolic_network_id
        # pubchem
        print db['pubchem_compound_id']
        # chebi
        print db['chebi_id']
        # kegg
        print db['kegg_id']
        # bigg
        print db['bigg_id']
        # HMDB
        # DrugBank
        print db['drugbank']
        print db['inchi'] #  
        print db['inchikey']
        print db['smiles']
        # protein associations
        print ";".join(db['protein_associations//protein//uniprot_id'])
        print db['average_molecular_weight']
        print db['chemical_formula']
        print "FormalCharge:", db["predicted_properties//property/.[kind='formal_charge']//value"]
        print "PhysiologicalCharge:", db["predicted_properties//property/.[kind='physiological_charge']//value"]



# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
