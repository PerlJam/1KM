"""
1KM unit tests
"""
import json
import factory
import os
import logging

from django.test import TestCase
from tastypie.test import ResourceTestCase, TestApiClient

import reports.tests
from reports.tests import MetaHashResourceBootstrap, assert_obj1_to_obj2, \
        find_all_obj_in_list, find_obj_in_list
import db.models
from reports.serializers import CSVSerializer
from test.factories import *


logger = logging.getLogger(__name__)


BASE_URI = '/db/api/v1'
BASE_REPORTS_URI = '/reports/api/v1'
import db; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__))
BASE_URI_DB = '/db/api/v1'


class TestApiInit(reports.tests.TestApiInit):
    # FIXME: factor out test_2api_init
    
    def setUp(self):
        super(TestApiInit, self).setUp()
#         super(TestApiInit, self)._setUp()
#         super(TestApiInit, self).test_0bootstrap_metahash()
        # NOTE: run reports tests to set up the metahash resource
        
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')

        
    def test_2api_init(self):
        
        print '***================ super: db test_2api_init =============== '

        super(TestApiInit, self).test_2api_init()
        
        print '***================ local: db test_2api_init =============== '
        
        
        serializer=CSVSerializer() 
        # todo: doesn't work for post, see TestApiClient.post() method, it is 
        # incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=serializer) 
        
        filename = os.path.join(self.db_directory,'api_init_actions.csv')
        with open(filename) as input_file:
            api_init_actions = serializer.from_csv(input_file.read(), root=None)
            
            bootstrap_files = [
                'metahash_resource_data.csv',
                'vocabularies_data.csv']
            for action in api_init_actions:
                
                print '\n++++=========== processing action', json.dumps(action)
                command = action['command'].lower() 
                resource = action['resource'].lower()
                resource_uri = BASE_REPORTS_URI + '/' + resource
                
                if command == 'delete':
                    resp = testApiClient.delete(
                        resource_uri, authentication=self.get_credentials())
                    self.assertHttpAccepted(resp)
                
                else:
                    filename = os.path.join(self.db_directory,action['file'])
                    search_excludes = []
                    # exclude 'resource_uri' from equivalency check during 
                    # bootstrap, because resource table needs to be loaded for
                    # the uri generation
                    if action['file'] in bootstrap_files:
                        search_excludes = ['resource_uri'] 
                    logger.info(str(('+++++++++++processing file', filename)))
                    with open(filename) as data_file:
                        input_data = serializer.from_csv(data_file.read())
                        
                        if command == 'put':
                            resp = testApiClient.put(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
                            logger.debug(str(('action: ', json.dumps(action), 
                                              'response: ' , resp.status_code)))
                            self.assertTrue(
                                resp.status_code in [200], str((resp.status_code, resp)))
                            
                            # now see if we can get these objects back
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), data={ 'limit': 999 })
                            self.assertTrue(
                                resp.status_code in [200], str((resp.status_code, resp)))
                            #   self.assertValidJSONResponse(resp)
                            new_obj = self.deserialize(resp)
                            result, msgs = find_all_obj_in_list(
                                input_data['objects'], new_obj['objects'], 
                                excludes=search_excludes)
                            self.assertTrue(
                                result, str((command, 'input file', filename, 
                                             msgs, new_obj['objects'])) )
                        
                        elif command == 'patch':
                            resp = testApiClient.patch(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
#                             self.assertHttpAccepted(resp)
                            self.assertTrue(resp.status_code in [202, 204], str((
                                'response not accepted, resource_uri:', resource_uri, 
                                'response', resp)))
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), data={ 'limit': 999 } )
                            self.assertTrue(
                                resp.status_code in [200], str((resp.status_code, resp)))
                            #                             self.assertValidJSONResponse(resp)
                            new_obj = self.deserialize(resp)
                            with open(filename) as f2:
                                input_data2 = serializer.from_csv(f2.read())
                                result, msgs = find_all_obj_in_list(
                                    input_data2['objects'], new_obj['objects'], 
                                    excludes=search_excludes)
                                self.assertTrue(
                                    result, str(( command, 'input file', filename, msgs )) )
                        
                        elif command == 'post':
                            self.fail((
                                'resource entry: ' + json.dumps(action) + '; '
                                'cannot POST multiple objects to tastypie; '
                                'therefore the "post" command is invalid with '
                                'the initialization scripts'))
                        else:
                            self.fail('unknown command: ' + command + ', ' + json.dumps(action))


class LibraryContentLoadTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        print '============== LibraryContentLoadTest setup ============'
        super(LibraryContentLoadTest, self).setUp()
        super(LibraryContentLoadTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(LibraryContentLoadTest, self)._bootstrap_init_files()
        print '============== LibraryContentLoadTest setup: begin ============'
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
    def test_load_sdf(self):
        pass;
        

class SmallMoleculeTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        print '============== SmallMoleculeTest setup ============'
        super(SmallMoleculeTest, self).setUp()
        super(SmallMoleculeTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(SmallMoleculeTest, self)._bootstrap_init_files()
        print '============== SmallMoleculeTest setup: begin ============'
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
        testApiClient = TestApiClient(serializer=self.csv_serializer) 

        filename = os.path.join(self.db_directory,'metahash_fields_smallmolecule.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.smallmolecule'})

        print '============== SmallMoleculeTest setup: done ============'
        

    def test1_create_smallmolecule(self):
        logger.info(str(('==== test1_create_smallmolecule =====')))
        
        resource_uri = BASE_URI_DB + '/smallmolecule'
        
        item = SmallMoleculeFactory.attributes()
        
        logger.info(str((item)))
        resp = self.api_client.post(
            resource_uri, format='json', 
            data=item,authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        # create a second smallmolecule
        item = SmallMoleculeFactory.attributes()
        
        logger.info(str(('item', item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp)))
        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(item, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         item, new_obj['objects'])))
        logger.info(str(('item found', obj)))

    def test2_create_register(self):

        logger.info(str(('==== test2_create_register =====')))
        
        resource_uri = BASE_URI_DB + '/smallmolecule'
        
        item = SmallMoleculeFactory.attributes()
        
        item['sm_id'] = None
        
        logger.info(str(('item to post', item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        # create a second smallmolecule
        item = SmallMoleculeFactory.attributes()
        item['sm_id'] = None
        
        logger.info(str(('item to post', item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))

        result, obj = find_obj_in_list(item, new_obj['objects'],
                                       excludes=['sm_id'] )
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         item, new_obj['objects'])))
        logger.info(str(('item found', obj)))
        self.assertTrue(not obj['sm_id']==None, 
                        str(('obj returned has no sm_id', obj)))
