#! env python
"""Converts ids from one public database to another

Uses:
 http://cts.fiehnlab.ucdavis.edu/service/convert/pubchem%20cid/human%20metabolome%20database/122347

   From values:
   http://cts.fiehnlab.ucdavis.edu/service/convert/fromValues

   To values:
    http://cts.fiehnlab.ucdavis.edu/service/convert/toValues

    Author: Samuel A. Inverso, Wyss Institute, 2014-06-25
"""

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter

import time

import json
import urllib2

from pprint import pprint

import logging

logger = logging.getLogger('1km')

VALID_DATABASE_TO_URL = 'http://cts.fiehnlab.ucdavis.edu/service/convert/toValues'
VALID_DATABASE_FROM_URL = 'http://cts.fiehnlab.ucdavis.edu/service/convert/fromValues'
CONVERT_ID_URL =  'http://cts.fiehnlab.ucdavis.edu/service/convert/%s/%s/%s'

WAIT_TIME = 0.1

def validate_convert_args(args):
    """Make sure all convert args are there"""
    if not args.from_db or not args.to_db or not args.id:
        sys.stderr.write('ERROR: missing from_db, to_db, and/or id\n')
        exit(-1)


def print_valid_database_ids():
    """Prints the valid database ids"""
    print  "*"*10, "D A T A B A S E   F R O M:", "*"*10, "\n"
    valid_from = json.load(urllib2.urlopen(VALID_DATABASE_FROM_URL))
    print "\n".join(valid_from)

    print "\n", "*"*10, "D A T A B A S E   T O", "*"*10,"\n"
    valid_to = json.load(urllib2.urlopen(VALID_DATABASE_TO_URL))
    print "\n".join(valid_to)
   

def convert_ids(from_db,to_db,ids):
    """Takes a list of ids to convert"""
    num_ids =len(ids)
    result = {}
    for ind in xrange(num_ids):
        id_orig = ids[ind]
        result[id_orig] = (convert_id(from_db,to_db,id_orig))
        if ind < num_ids:
            time.sleep(WAIT_TIME)
    return result


def convert_id(from_db,to_db,id_orig):
    quoted_params = map(lambda x: urllib2.quote(x), (from_db,to_db,id_orig))
    url_str = CONVERT_ID_URL % tuple(quoted_params)
    logger.debug("convert_id url: " + url_str)
    new_ids = json.load(urllib2.urlopen(url_str))
    return new_ids[0]['result']


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
        help='Database to convert from.')
    group1.add_argument('to_db', metavar='to', action='store', nargs='?',
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
