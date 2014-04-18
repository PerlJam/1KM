#! python
"""Converts a genbank file into a table loadable to the Database

   Author: Samuel A. Inverso, Wyss Institute, 2014-04-16
"""


# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch
import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from Bio import SeqIO 
from Bio import Entrez
import pandas as pd

def getFromKeyList(lst,key,sep=":",firstOnly=False):
    """ Parases a list of strings that have a deliminator
    """
    result = []
    values = []
    for itm in lst:
        if itm.startswith(key):
            values.append(itm.partition(sep)[2])

    if firstOnly:
        if len(values) > 0:
            result = values[0]
        else:
            result = None
    else:
        result = values

    return result

def genbank_dataframe_to_database_csv(df, output_dir, sep=','):
    """Takes a genbank datafram and converts it into csv's
        That can be loaded into the database
    """
    gene_columns = ['gene_ID', 'Alternate_Name', 'Sequence' ]
    protein_columns =  ['pp_ID', 'Amino_Acid_Sequence', 'gene_ID']


    gene_file_name = "gene.csv"
    protein_file_name = "protein.csv"

    # TODO test if directory exists if not create it
    
    # Gene
    output_path = os.path.join(output_dir, gene_file_name) 
    df[gene_columns].to_csv(output_path,sep)

    # protein
    output_path = os.path.join(output_dir, protein_file_name)
    df[protein_columns].to_csv(output_path,sep)

def genbank_feature_to_dict(feature, record):
    """Take a genbank feature and put it into a dictionary with
        that can be further processed into a database table
        
        Ony processes CDS entries ignores gene 
    """
    result = {}
    
    if feature.type == 'CDS':
        # Gene
        result['gene_ID'] = feature.qualifiers['locus_tag'][0]
        result['Alternate_Name'] = feature.qualifiers['locus_tag'][0]
        result['Sequence'] = str(feature.extract(record.seq))
        result['Coding_Indicies'] = str(feature.location)

        # Protein
        result['UniProt_ID'] = getFromKeyList( 
                          feature.qualifiers['db_xref'],"UniProtKB/Swiss-Prot",
                          firstOnly=True )

        result['pp_ID'] = getFromKeyList( 
                          feature.qualifiers['db_xref'],"UniProtKB/Swiss-Prot",
                          firstOnly=True )

        if 'translation' in feature.qualifiers: 
            result['Amino_Acid_Sequence'] = feature.qualifiers['translation'][0]
        else:
            with warnings.catch_warnings(record=True) as ws:
                result['Amino_Acid_Sequence'] = \
                    feature.extract(record.seq).translate()
                for w in ws:
                    if 'warning' not in result:
                        result['warning'] = []
                    
                    result['warning'].append(w.message)


    return result 

def genbank_record_to_dataframe(record):
    """Converts a genbank Bio:Record to a dictionary that
        can be then converted into a CSV

        Uses the database mapping.
    """

    rows = []
    for feature in record.features:
        row = genbank_feature_to_dict(feature,record)
        if row:
            rows.append(row)

    return pd.DataFrame(rows)

def genbank_file_to_record(filename):
    record = SeqIO.read(filename,"genbank")
    return genbank_record_to_dataframe(record)

def main(argv):
    """Main entry point of program.
    """

    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Converts a genebank file into a series of csv files for loading'
        ' into the database.\n',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('gb', metavar='genbank', action='store',
        help='Genbank file to process')
    parser.add_argument('output_dir', metavar='output_dir', action='store',
        help='Output directory to place csvs into.')

    args = parser.parse_args()
    genbank_file = args.gb
    output_dir = args.output_dir


    # Convert the records to a dataframe
    df = genbank_file_to_record(genbank_file)

    # turn the dataframe into a series of tables
    genbank_dataframe_to_database_csv(df,output_dir)


# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
