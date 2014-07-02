#! env python
"""Pointed to an HMDB directory with xml files in it, it creates and index
    file of ids to xmls files. This is because one file can have more
    than on id inside of it, it takes a long time to search 
    through the files if the id happens to not be the filename
    
    The index file is a dictionary imported from hmdb_index.py

   For usage type:
   python hmdb_create_index.py -h 
   
   # from src directory
   python hmdb_create_index.py  ../tests/hmdb_metabolites 

   Author: Samuel A. Inverso, Wyss Institute, 2014-06-18
"""


# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch
import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from pprint import pprint

import xml.etree.cElementTree as ET




def main(argv):
    """Main entry point of program.
    """

    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Pointed to an HMDB directory with xml files in it, it creates and index'
        ' file of ids to xmls files. This is because one file can have more'
        ' than on id inside of it, it takes a long time to search' 
        ' through the files if the id happens to not be the filename.\n'    
        ' The index file is a dictionary imported from hmdb_index.py.\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('input_dir', action='store',
        help='Directory with xml files to index.')
    parser.add_argument('-o', metavar='output_dir', action='store',
        default = '.',
        help='Output directory to place hmdb_index.py., default=.')

    output_filename = "hmdb_index.py"

    # Parse the args
    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.o

    # loop through the directory and find the accession to build the index
    index = {} # id : file
    for fn in os.listdir(input_dir):
        file_path = os.path.join(input_dir,fn)
        if os.path.isfile(file_path) and file_path.endswith('.xml'):
            print file_path
            tree = ET.ElementTree(file=file_path)

            for elem in tree.iter():
                if elem.tag == 'accession' and elem.text.startswith('HMDB'):
                    index[elem.text] = fn

    # write out the index file
    with open(os.path.join(output_dir,output_filename),"w") as fout:
        fout.write('# HMDB index, id : filename \n') 
        fout.write('HMDB_DIRECTORY = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../dat/hmdb_metabolites")\n')
        fout.write('HMDB_INDEX = ')
        pprint(index, fout)

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
