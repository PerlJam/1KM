#! env python
"""Fetchs record from chebi  
    Uses bioservices
    http://pythonhosted.org/bioservices/references.html?highlight=chebi#module-bioservices.chebi
   Author: Samuel A. Inverso, Wyss Institute, 2014-06-10
"""

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from pprint import pprint

import pandas as pd

from bioservices import *

def from_id(chebiId):
    """This is a pass through to ch = ChEBI()
    res = ch.getCompleteEntity(chebiId)
    """
    ch = ChEBI()
    # Check that the ID has "CHEBI:" in front of it, if not prepend it
   # print chebiId, isinstance(chebiId,basestring)

    if not isinstance(chebiId,basestring):
        chebiId = repr(int(chebiId)) # case where it's a numpy.float64

    # make sure it begins with CHEBI:
    if not chebiId.startswith('CHEBI:'):
        chebiId = 'CHEBI:' + chebiId

    return ch.getCompleteEntity(chebiId)


def main(argv):
    """Main entry point of program.
    """
    
    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Retrieves compounds from Chebi using a list of ChebiIds\n'
        'For example, \n'
        'python chebi_fetch "CHEBI:27732"',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('chebiId', metavar='chebiId', nargs='+',action='store',
        help='CheEBI id to fetch.')
   # parser.add_argument('output_dir', metavar='output_dir', action='store',
   #     help='Output directory to place csvs into.')

    # parser.add_argument('-s', nargs=2, metavar=('identifier', 'type'), action='store', dest='search',
    #     help='Identifier and type to search for: e.g. -s Glucose name\n'
    #         '  Types supported: cid, name, smiles, sdf, inchi, inchikey, formula.\n')

    args = parser.parse_args()
    chebiId = args.chebiId
   
    ch = ChEBI()
    c = ch.getCompleteEntity(chebiId)
    print c

    print c.chebiId #
    print c.chebiAsciiName
    print ";".join([x['data'] for x in c.Synonyms]  )

    print 'Iupac', ";".join([x['data'] for x in c.IupacNames] ) #
    # alternate name
    # metabolic_network_id
    # chebi
    # kegg
    for dbl in c.DatabaseLinks:
        if "KEGG COMPOUND accession"  == dbl['type']:
            kegg = dbl['data']
            print kegg
        elif "DrugBank accession" == dbl['type']:
            drugbank = dbl['data']
            print drugbank
        elif "HMDB accession" == dbl['type']:
            hmdb = dbl['data']
            print hmdb


    # bigg
    # HMDB
    # DrubBank
    print c.inchi #  
    print c.inchiKey #
    print c.smiles #
    # protein associations
    print c.mass #
    print ";".join([x['data'] for x in c.Formulae]) #
    print c.charge #


#     print "Done!"

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
