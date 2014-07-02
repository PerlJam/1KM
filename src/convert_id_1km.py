#! env python
"""Converts ids from one public database to another 

Uses: 1km database id names
    
    Author: Samuel A. Inverso, Wyss Institute, 2014-06-25
"""

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter

from pprint import pprint

import convert_id_fiehnlab
import fix_value

DB_TO_FIEHNLAB_DB = {'PubChem_CID':'PubChem CID',
    'ChEBI_ID':'ChEBI', 'KEGG_ID':'KEGG', 'BiGG_ID': None, 
    'HMDB_ID':'Human Metabolome Database', 'DrugBank_ID':'DrugBank', 
    'InChi_Parent': 'InChI Code', 'InChi_Key_Parent':'InChIKey'}


def validate_convert_args(args):
    """Make sure all convert args are there"""
    if not args.from_db or not args.to_db or not args.id:
        sys.stderr.write('ERROR: missing from_db, to_db, and/or id\n')
        exit(-1)


def print_valid_database_ids():
    """Prints the valid database ids"""
    # get the dbs we support, one's with None have no mapping
    dbs = [ k for k,v in DB_TO_FIEHNLAB_DB.iteritems() if v is not None]
    dbs = sorted(dbs)
    print "\n".join(dbs)
   

def convert_ids(from_db,to_db,ids):
    """Takes a list of ids to convert"""
    result = {}
    for id_value in ids:
        result[id_value] = convert_id(from_db,to_db, id_value)

    return result


def convert_id(from_db,to_db,id_value):
    """Takes an id to convert"""
    if is_valid_db(from_db) and is_valid_db(to_db):
        result = convert_id_fiehnlab.convert_id(
            DB_TO_FIEHNLAB_DB[from_db],DB_TO_FIEHNLAB_DB[to_db], id_value)

        if 'ChEBI_ID' == to_db and result:
            result = map( fix_value.fix_chebi_id, result )
        elif 'InChi_Key_Parent' == to_db and result:
            result = map( fix_value.fix_inchi_parent_key, result)

    else:
        result = None

    return result


def is_valid_db(db):
    return db in DB_TO_FIEHNLAB_DB and DB_TO_FIEHNLAB_DB[db] 


def validate_db(db):
    if not is_valid_db(db):
        raise argparse.ArgumentTypeError("%s is an unsupported database" % (db))
    return db


def main(argv):
    """Main entry point of program.
    """
    # Parse arguments
    # setup two parse groups, one that is for convert, and the other to
    # list valid to and from databases
    parser = argparse.ArgumentParser( usage='%(prog)s [-h] [-l] from to id [id ...]',
        description=
        'Converts from one database id to another\n',
        formatter_class=RawTextHelpFormatter)

    group1 = parser.add_argument_group("Convert")
    group1.add_argument('from_db', metavar='from', action='store', nargs='?',
            type=validate_db,
        help='Database to convert from.')
    group1.add_argument('to_db', metavar='to', action='store', nargs='?',
        type=validate_db,
        help='Databse to convert to.')
    group1.add_argument('id', metavar='id', nargs="*", action='store',
        help='id(s) to convert')

    group2 = parser.add_argument_group("List valid database")
    group2.add_argument('-l', action='store_true',
        help='List valid from and to databases. ')

    # make sure we have from to and id if -l is not specified
    args = parser.parse_args()
    is_print_list = args.l
    if is_print_list:
        print_valid_database_ids()
    else:
        validate_convert_args(args)
        from_db = args.from_db
        to_db = args.to_db
        ids = args.id
        converted_id = convert_ids(from_db,to_db,ids)

        # display the result to the user
        pprint( converted_id )

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
