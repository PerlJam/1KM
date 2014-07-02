#! env python
"""Fetchs record from drugbank  
   Author: Samuel A. Inverso, Wyss Institute, 2014-06-30
"""

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from pprint import pprint

import pandas as pd

from bioservices import *

import xml.etree.cElementTree as ET

from drugbank_config import DRUGBANK_PATH

import logging
logger = logging.getLogger('1km')

class DrugBank:
    NAMESPACES = {'drugbank': 'http://www.drugbank.ca'}
    _top_root = None # Top root of the parsesd xml
    _cache_record = {} # entries by cid 

    def __init__(self,cid):
        """Initialze the class, load the drugbank xml file if it's not yet"""
        self._cache = {}

        # Check if we've parased the drugbank xml file yet
        if not DrugBank._top_root:
            logger.debug("DrugBank: loading: " + DRUGBANK_PATH)
            DrugBank._top_root = ET.parse(DRUGBANK_PATH).getroot()
        
        # check if we've looked for this cid before
        if cid in DrugBank._cache_record:
            self._root = DrugBank._cache_record[cid]
        else:
            self._root = self._get_drug_entry(cid)
            DrugBank._cache_record[cid] = self._root

    def formula(self):
        """The molecular formula
            Can exist in two places calculated-properties and experimental-properties
        """
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='Molecular Formula']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='Molecular Formula']/drugbank:value" ]

        return self._get_property(xpaths)
    

    def name(self):
        xpaths = ['drugbank:name']
        return self._get_property(xpaths)


    def molecular_weight(self):
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='Molecular Weight']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='Molecular Weight']/drugbank:value" ]
        return self._get_property(xpaths)


    def smiles(self):
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='SMILES']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='SMILES']/drugbank:value" ]
        return self._get_property(xpaths)


    def inchi(self):
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='InChI']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='InChI']/drugbank:value" ]
        return self._get_property(xpaths)


    def inchi_key(self):
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='InChIKey']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='InChIKey']/drugbank:value" ]
        return self._get_property(xpaths)


    def physiological_charge(self):
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='Physiological Charge']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='Physiological Charge']/drugbank:value" ]
        return self._get_property(xpaths)


    def iupac_name(self):
        xpaths = ["drugbank:calculated-properties/drugbank:property/.[drugbank:kind='IUPAC Name']/drugbank:value",
                    "drugbank:experimental-properties/drugbank:property/.[drugbank:kind='IUPAC Name']/drugbank:value" ]
        return self._get_property(xpaths)


    def pubchem_compound_id(self):
        xpaths = ["drugbank:external-identifiers/drugbank:external-identifier/.[drugbank:resource='PubChem Compound']/drugbank:identifier" ]
        return self._get_property(xpaths)


    def pubchem_substance_id(self):
        xpaths = ["drugbank:external-identifiers/drugbank:external-identifier/.[drugbank:resource='PubChem Substance']/drugbank:identifier" ]
        return self._get_property(xpaths)


    def kegg_compound(self):
        xpaths = ["drugbank:external-identifiers/drugbank:external-identifier/.[drugbank:resource='KEGG Compound']/drugbank:identifier" ]
        return self._get_property(xpaths)


    def kegg_drug(self):
        xpaths = ["drugbank:external-identifiers/drugbank:external-identifier/.[drugbank:resource='KEGG Drug']/drugbank:identifier" ]
        return self._get_property(xpaths)


    def chebi(self):
        xpaths = ["drugbank:external-identifiers/drugbank:external-identifier/.[drugbank:resource='ChEBI']/drugbank:identifier" ]
        return self._get_property(xpaths)


    def _get_property(self, xpaths):
        """Helper function for the properties, e.g. formula"""
        result = None
        for xpath in xpaths:
            result = self[xpath]
            if len(result) > 0:
                break

        return result


    def __getitem__(self,key):
        """Allows access through brackets, e.g.
         foo['drugbank:name']
         Make sure to use the namespace drugbank.
        """
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


    def _get_drug_entry(self,cid):
        """Finds the entry related this id"""
        result = None
        search_str = "drugbank:drug/drugbank:drugbank-id/.[@primary='true']/..[drugbank:drugbank-id='%s']"  % (cid) 
        result = DrugBank._top_root.find(search_str,DrugBank.NAMESPACES)
        return result


    def _get_from_tree(self,key):
        '''trys to get the key from the tree, cache's the result
        e.g. to get all the synonyms 
         key = 'synonyms//synonym' 
        '''
        result = []

        # There are some special cases 
        for node in self._root.findall(key,DrugBank.NAMESPACES):
            result.append(node.text)

        return result


    def clear_cache(self):
        """Clears all caching
            includes entryies and the parsed xml 
        """
        self._cache = {}
        DrugBank._cache_record = {}
        DrugBank._top_root = None


def from_id(cid):
    """Get's the DrugBank class given an id
    """
    return DrugBank(cid)


def main(argv):
    """Main entry point of program.
    """
    parser = argparse.ArgumentParser(description=
        'Retrieves compounds from DrugBank using either a list of ID e.g. DB02379\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('id', metavar='id', nargs='+',action='store',
        help='id(s) to fetch.')


    args = parser.parse_args()
    cids = args.id
  
    for cid in cids:
        db = from_id(cid)

        print db.name()
        # print ";".join(  db['synonyms//synonym']   )

        print db.iupac_name()
        # # alternate name
        # # metabolic_network_id
        # # pubchem
        print db.pubchem_compound_id()
        # # chebi
        print db.chebi()
        # # kegg
        print db.kegg_compound()
        # # bigg
        # print db['bigg_id']
        # # HMDB
        # # DrugBank
        # print db['drugbank']
        print db.inchi()
        print db.inchi_key() 
        print db.smiles()
        # # protein associations
        # print ";".join(db['protein_associations//protein//uniprot_id'])
        print db.molecular_weight()
        print db.formula()
        # print "FormalCharge:", db["predicted_properties//property/.[kind='formal_charge']//value"]

        print db.physiological_charge()
        print "chemical formula:",  db.formula()

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
