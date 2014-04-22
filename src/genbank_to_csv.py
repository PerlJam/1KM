#! env python
"""Converts a genbank file into a table loadable to the Database
    
   For usage type:
   python genbank_file_to_csv.py -h 
   
   # from src directory
   python genbank_to_csv.py ../tests/gi_U00096.gbk ../tests -f ../tests/blattner_filter.csv

   Author: Samuel A. Inverso, Wyss Institute, 2014-04-16
"""


# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch
import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter
from Bio import SeqIO 
from Bio import Entrez
import pandas as pd
import csv

def load_filter_file(filter_file):
    """Reads in a csv filter file

        filter_file: file name, or None
        returns a DataFrame with the csv data or an empty DataFrame if 
            filter_file is none
    """
    df = pd.DataFrame()
    # if
    if filter_file:
        with open(filter_file,'rb') as csvfile:
            df = pd.read_csv(csvfile)
    return df


def get_from_key_list(lst,key,sep=":",firstOnly=False):
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


def genbank_dataframe_to_database_csv(df, output_dir, sep=',',index=False):
    """Takes a genbank datafram and converts it into csv's
        That can be loaded into the database

        index : boolean, wrete row names (index), default False
    """
    gene_columns = ['gene_ID', 'Alternate_Name', 'Sequence' ]
    protein_columns =  ['pp_ID', 'Amino_Acid_Sequence', 'gene_ID']


    gene_file_name = "gene.csv"
    protein_file_name = "protein.csv"

    # TODO test if directory exists if not create it
    
    # Gene
    output_path = os.path.join(output_dir, gene_file_name) 
    df[gene_columns].to_csv(output_path,sep,index=index)

    # protein
    output_path = os.path.join(output_dir, protein_file_name)
    df[protein_columns].to_csv(output_path,sep,index=index)


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
        result['UniProt_ID'] = get_from_key_list( 
                          feature.qualifiers['db_xref'],"UniProtKB/Swiss-Prot",
                          firstOnly=True )

        result['pp_ID'] = get_from_key_list( 
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


def genbank_record_to_dataframe(record,dfFilter=None):
    """Converts a genbank Bio:Record to a dictionary that
        can be then converted into a CSV

        Uses the database mapping.
    """
    rows = []
    for feature in record.features:
        row = genbank_feature_to_dict(feature,record)
        if row:
            rows.append(row)

    df = pd.DataFrame(rows)
    
    df = filter_dataframe(df,dfFilter)

    return df


def filter_dataframe(df,dfFilter):
    """ Match columns between df and dfFilter. If the values in 
        df appear in dfFilter then we keep that row.

        Only columns that appear in both are used for filtering
        i.e. extra columns in dfFilter that are not in df are ignored.
    """
    result = df
    if dfFilter is not None and not dfFilter.empty:
        # for some reason using the dataframe as a filter is not working
        # so we are going to loop trough the columns
        col_names = list(dfFilter.columns.values)
        df_col_names = list(df.columns.values)
        # only loop through the columns that are in df
        mask = [True] * len(df.index)
        for col in  col_names:
            if col in df:
                mask &= df[col].isin(dfFilter[col])
    
        result = df[mask]

    return  result


def genbank_file_to_record(filename,dfFilter=None):
    """ Reads a genbank file and returns a BioPython record
    """
    record = SeqIO.read(filename,"genbank")
    return genbank_record_to_dataframe(record,dfFilter)


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

    parser.add_argument('-f', metavar='filter', action='store', dest='filter',
        help='A CSV file with column header the name of a tag to filter'
            ' on and rows with values to keep.'
            ' The row header should be the name that will appear in the final'
            ' database table. For example, to filter on blattner number or an'
            ' ecoli genbank file, make a CSV file with one column'
            '     Alternate_Name\n'
            '     b0002\n'
            '     b0003\n'
            '\n'
            'python genbank_file_to_csv gi_U00096.gbk . -f blattner_filter.csv\n'
            'Limits to entries with blattner numbes b0002 and b0003.' )

    args = parser.parse_args()
    genbank_file = args.gb
    output_dir = args.output_dir
    filter_file = args.filter


    # get the filter file if it exists
    dfFilter = load_filter_file(filter_file)

    # Convert the records to a dataframe
    df = genbank_file_to_record(genbank_file, dfFilter)

    # turn the dataframe into a series of tables
    genbank_dataframe_to_database_csv(df,output_dir)

# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
