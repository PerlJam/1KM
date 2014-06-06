#! env python
"""Fetches parts and sequences from the synberc registry.
    registry.synberc.org

    entry ids are in ../test/synberc_entry_ids.txt
    
   Author: Samuel A. Inverso, Wyss Institute, 2014-06-02
"""
# http://stackoverflow.com/questions/21280978/saxparseexception-on-file-header-with-suds-consuming-soap-service-how-to-consu
from suds.client import Client
from suds.wsse import *
import logging
from suds.plugin import MessagePlugin

import argparse, os, sys, warnings
from argparse import RawTextHelpFormatter

import pandas as pd
import csv

import time

from progressbar import ProgressBar

#logging.basicConfig(level=logging.INFO)
#logging.getLogger('suds.client').setLevel(logging.DEBUG)

# The entries in the part record to keep and put into the csv
ENTRIES_TO_KEEP = ('bioSafetyLevel',
   'creator',  
   'keywords',
   'name',
   'partId', 
   'principalInvestigator',
   'selectionMarkers',
   'shortDescription',
   'type',
   'backbone', 
   'originOfReplication',
   'promoters',
   'replicatesIn',
   'recordId',
   'sequence' )


def part_to_dict(part,keep=None):
    """Converts a synberc record and sequence into a dictionary of key values
        that can be later turned into a pandas dataframe.

        Keep: tuple or list of entries to keep, e.g. ENTRIES_TO_KEEP. if None,
            then all entries are kept 
    """ 

    result = {}

    if keep is None:
      keep = ENTRIES_TO_KEEP

    for key in keep:
        if key in part:
            result[key] = part[key]
        else:
            result[key] = None

    return result


def dicts_to_dataframe(d):
    """Converts a list of dictionaries to a pandas dataframe
    """
    return pd.DataFrame(d)


class Filter(MessagePlugin):
    ''' Synberc returns an HTTP header instead of just the xml. We need to 
        parse it.
    --uuid:a775142b-b4e0-4ebf-8d94-f258caf2af70
    Content-Type: application/xop+xml; charset=UTF-8; type="text/xml";
    Content-Transfer-Encoding: binary
    Content-ID: <root.message@cxf.apache.org>

    <soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/"><soap:Body><ns2:loginResponse xmlns:ns2="https://api.registry.jbei.org/"><return>f2e4d1b8bcc0e7a298cdee272dfb8bf3d37f94b1</return></ns2:loginResponse></soap:Body></soap:Envelope>
    --uuid:a775142b-b4e0-4ebf-8d94-f258caf2af70--
    '''

    def getsoap(self,text):
        return text[text.find('<soap:Envelope'):text.rfind('>')+1]


    def received(self, context):
        context.reply = self.getsoap(context.reply)


class Synberc():
    '''Handles logging in and staying logged in while retrieving data from 
        SynBerc
    '''
    login = None
    username = None
    password = None
    wsdl_file = None
    session_id = None
    client = None
    request_wait = None # how long to wait between requests in seconds

    def __init__(self,username,password,wsdl_file,request_wait=0.5):
        '''Contructor, logins automatically
        '''

        self.username = username
        self.password = password
        self.wsdl_file = wsdl_file
        self.request_wait= request_wait

        # connect to the server
        self.client = Client(wsdl_file,plugins=[Filter()])

        # login
        self.login()
     

    def login(self):
        '''Login with username and pass, throws an exception if we can't login
        '''
        # check if we're authenticated 
        self.session_id = self.client.service.login(login=self.username,
            password=self.password)

        # check that we're logged in
        if not self.is_authenticated():
            raise Exception('Could not authenticate to Synberc')


    def is_authenticated(self):
        '''Check if we're authenticated
        '''
        return self.client.service.isAuthenticated(self.session_id)


    def get_part_and_sequence(self, entry_ids,keep): 
        ''' list of entry_ids and returns 
            a list of records with sequence imbeded
        '''
        result = []
        # check if we're still authenticaed 
        if not self.is_authenticated():
            # try to reauthenticate
            self.login()

        # get the part and sequence for each eid
        pbar = ProgressBar()
        first = True # first time through loop
        for eid in pbar(entry_ids):
            # Wait if this isn't the first time through the loop
            if not first:
                time.sleep(self.request_wait)
                first = False

            # get the part
            part = self.client.service.getPartByRecordId(self.session_id,eid)

            # wait 
            time.sleep(self.request_wait)

            # get the sequence
            sequence = self.client.service.getSequence(self.session_id,eid)

            # join the part and sequence, convert the whole thing to a dict
            # and append it to the results
            if sequence is not None:
                part['sequence'] = sequence.sequence
            else:
                part['sequence'] = ''

            result.append(part_to_dict(part,keep))
    
        return result


def test_synberc():
    ''' Quick tester to check that wsdl methods work
        Need to change password to the actual password
    ''' 
    wsdlFile = 'https://registry.synberc.org/api/RegistryAPI?wsdl'
    namespace = 'https://api.registry.jbei.org/'
    username = 'Pamela Silver'
    password =  None # need password
    entryId = '1b5e7355-4fa6-4cd5-9ae5-06ef38f841e2'
    # connect to server
    client = Client(wsdlFile,plugins=[Filter()])
    print client

    # Call login
    print "Logging in:"
    sessionId = client.service.login(login=username,password=password)
    print "Sesion ID", sessionId


    # Get part by ID
    print "Getting part by id"
    part = client.service.getPartByRecordId(sessionId,entryId)
    print part
    print part['creator']
    # get seequence

    print "Getting sequence"
    sequence = client.service.getSequence(sessionId,entryId)
    print sequence

    print "The actual sequence"
    print sequence.sequence
    print "Done!"


def load_entry_ids(filename):
    """Loads record ids from a file with one record per line
    """
    result = []

    if filename is not None:
        with open(filename) as f:
            result = f.read().splitlines()

    return result


def main(argv):
    """Main entry point of program.
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description=
        'Retrieves parts and sequences from the SynBerc registry requires the record id\n'
        'Example:\n'
        'python synberc.py -u username -p password -f ../tests/synberc_entry_ids.txt -o synberc_parts.csv',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-r', metavar='entry_id', nargs='+',action='store',
        dest='entry_id', default=[],
        help='Record ids to fetch.')

    parser.add_argument('-f', metavar=file, action='store', dest='file',
        help='File with one record id per line\n')

    parser.add_argument('-w', metavar="wsdl_file", action='store', 
        dest='wsdl_file', default='https://registry.synberc.org/api/RegistryAPI?wsdl',
        help='url or file address for wsdl')

    parser.add_argument('-u', metavar="username", action='store', 
        required=True,
        dest='username', help='Username to login to synberc.')

    parser.add_argument('-p', metavar="password", action='store',
        required=True,
        dest='password', help='password to login to synberc.')

    parser.add_argument('-o', metavar='output_file', action='store',
        dest="output_file", help='Output file for csv.')

    parser.add_argument('-s', metavar='seconds', action='store',
        default=0.1,
        dest="request_wait", help='Number of seconds to wait between requests')

    args = parser.parse_args()
   
    # print help if there were no options given
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # get args
    wsdl_file = args.wsdl_file #'https://registry.synberc.org/api/RegistryAPI?wsdl'
    username = args.username 
    password = args.password
    output_file = args.output_file
    entry_ids = args.entry_id
    entry_id_file = args.file
    request_wait = args.request_wait

    # 
    # Concate records listed on command line with any listed in the file
    #
    entry_ids += load_entry_ids(entry_id_file)


    # Login to the database
    db = Synberc(username,password,wsdl_file,request_wait)

    # 
    # get the entries
    #
    parts = db.get_part_and_sequence(entry_ids, ENTRIES_TO_KEEP)

    # make it a panda data frame for output 
    parts_df = dicts_to_dataframe(parts)
 
    # if we're supposed to write this as a csv convert it and write
    # otherwise display the result to the user as a dataframe
    if output_file:
        parts_df.to_csv(output_file,index=False, encoding='utf-8')
    else:
        print parts_df


# Call the main method if we're run as script
if __name__ == "__main__":
    main(sys.argv)
