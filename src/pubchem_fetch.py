#! env python
"""Fetchs record from pubchem  
    basically a pass through to PubChemPy
    http://pubchempy.readthedocs.org/en/latest/guide/gettingstarted.html

   Author: Samuel A. Inverso, Wyss Institute, 2014-04-16
"""

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from pprint import pprint

import pandas as pd

from pubchempy import Compound, Substance, get_compounds

import urllib2

import logging

from pprint import pprint

# PubChemPy by default logs debug, we don't want that now.
logger = logging.getLogger('pubchempy')
logger.setLevel(level=logging.INFO)

def from_id(cid):
    """This is a pass through to pubchempy.Compound.from_cid(cid)
    """
    # sometimes cid comes in as a float, it should be an int
    if isinstance(cid,float):
        cid = int(cid)

    return Compound.from_cid(cid)


def from_sid(sid):
    """This is a pass trhough to pubchempy.Substance.from_sid(sid)"""
    # sometimes cid comes in as a float, it should be an int

    if isinstance(sid,float):
        sid = int(sid)

    return Substance.from_sid(sid)


def sid_to_cid(sid):
    """Attempts to find sid's corresponding cid, usesful the pubchem rest interface,
        sid 

        returns an int
    """
    if isinstance(sid,float):
        sid = int(sid)
    response = None
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/" + \
         str(sid) +  "/cids/TXT?cids_type=standardized"
    try:
        response = urllib2.urlopen(url).read()
    except urllib2.HTTPError as e:
        print e.reason

    return int(response)


def main(argv):
    """Main entry point of program.
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Retrieves compounds from PubChem using either a list of CID''s or'
        ' identifier and type, e.g. Glucose name .\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('cid', metavar='cid', nargs='+',action='store',
        help='CIDs to fetch.')
   # parser.add_argument('output_dir', metavar='output_dir', action='store',
   #     help='Output directory to place csvs into.')

    parser.add_argument('-s', nargs=2, metavar=('identifier', 'type'), action='store', dest='search',
        help='Identifier and type to search for: e.g. -s Glucose name\n'
            '  Types supported: cid, name, smiles, sdf, inchi, inchikey, formula.\n')

    args = parser.parse_args()
    cid = args.cid
    #output_dir = args.output_dir
    search = args.search

    c = Compound.from_cid(cid, as_dataframe=True)
    
    pprint(c.record)
    print "Synonyms", ";".join(c.synonyms)

    print c.iupac_name
    # alternate name
    # metabolic_network_id
    print c.cid
    # chebi
    # kegg
    # bigg
    # HMDB
    # DrubBank
    print c.inchi
    print c.inchikey
    print c.canonical_smiles
    # protein associations
    print c.exact_mass
    print c.molecular_weight
    print c.molecular_formula
    print c.charge

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
