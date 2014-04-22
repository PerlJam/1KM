from __future__ import unicode_literals

import logging
import json

from tastypie.utils import timezone
from tastypie.exceptions import BadRequest, Unauthorized
from tastypie.utils.urls import trailing_slash
from tastypie.authorization import Authorization, ReadOnlyAuthorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication,\
    MultiAuthentication
from tastypie.constants import ALL_WITH_RELATIONS
from tastypie import fields

from django.db import models
from django.conf.urls import url

from reports.serializers import CursorSerializer, LimsSerializer, SmallMoleculeSerializer
from reports.models import MetaHash, Vocabularies, ApiLog
from reports.api import ManagedModelResource, ManagedResource, ApiLogResource, \
        SuperUserAuthorization
from db.models import SmallMolecule, Reaction, Gene, Protein


import logging

logger = logging.getLogger(__name__)

    
class SmallMoleculeResource(ManagedModelResource):

    class Meta:
        queryset = SmallMolecule.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= SuperUserAuthorization()        
        resource_name = 'smallmolecule'
        
        # NOTE: in order to patch_list, wherein 'resource_uri' is not set, 
        # the method 'put' is required.  TODO: figure out better allowed methods
#         allowed_methods = ['get', 'patch', 'delete', 'put', 'post']
#         list_allowed_methods = ['get','patch','put','delete']
        
        always_return_data = True
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        
    def __init__(self, **kwargs):
        super(SmallMoleculeResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # TODO: create a custom mapper for schema before others to avoid 
        # complicated regex
        return [
            url( (r"^(?P<resource_name>%s)/log/(?P<apilog>\d+)%s$" ) 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_smallmolecule_apilog_view'), 
                name="api_dispatch_smallmolecule_apilog_view"),
            url((r"^(?P<resource_name>%s)/(?P<sm_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
                ]

    def dispatch_smallmolecule_apilog_view(self, request, **kwargs):
        return self.dispatch('list', request, **kwargs)     

    def obj_get_list(self, bundle, **kwargs):
        
        # TODO: this override is implemented as a convenience to access the 
        # "batch" updated at one time.  Refactor so that the ApiLogResource
        # can redirect to -any- resource
        apilog_key = None
        if 'apilog' in kwargs:
            apilog_key = kwargs.pop('apilog')
        
        obj_query = super(SmallMoleculeResource, self).obj_get_list(bundle, **kwargs)
        
        # This is a hook to get items for the /log/log# uri
        if apilog_key:
            apilog = ApiLog.objects.get(pk=apilog_key)
            logger.info(str(('apilog found', apilog)))
            
            if apilog.api_action == 'PATCH_LIST' and apilog.added_keys:
                keys = json.loads(apilog.added_keys)
                logger.info(str(('query for keys', keys)))
                # FIXME: on refactor; grab the key field from the resource definition
                # FIXME: on refactor, will have to deal with composite keys
                obj_query = obj_query.filter(sm_id__in=keys)
        
        return obj_query
    
    def put_list(self, request, **kwargs):
        return super(SmallMoleculeResource, self).put_list(request, **kwargs) 
    
    def obj_create(self, bundle, **kwargs):
        logger.info(str(('===creating smallmolecule'))) #, bundle.data)))
        bundle = super(SmallMoleculeResource, self).obj_create(bundle, **kwargs)
        return bundle
    
class GeneResource(ManagedModelResource):

    class Meta:
        queryset = Gene.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= SuperUserAuthorization()        
        resource_name = 'gene'
        
        always_return_data = True
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        
    def __init__(self, **kwargs):
        super(GeneResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # TODO: create a custom mapper for schema before others to avoid 
        # complicated regex
        return [
            url( (r"^(?P<resource_name>%s)/log/(?P<apilog>\d+)%s$" ) 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_gene_apilog_view'), 
                name="api_dispatch_gene_apilog_view"),
            url((r"^(?P<resource_name>%s)/(?P<gene_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
                ]

    def dispatch_gene_apilog_view(self, request, **kwargs):
        return self.dispatch('list', request, **kwargs)     

    def obj_get_list(self, bundle, **kwargs):
        
        # TODO: this override is implemented as a convenience to access the 
        # "batch" updated at one time.  Refactor so that the ApiLogResource
        # can redirect to -any- resource
        apilog_key = None
        if 'apilog' in kwargs:
            apilog_key = kwargs.pop('apilog')
        
        obj_query = super(GeneResource, self).obj_get_list(bundle, **kwargs)
        
        # This is a hook to get items for the /log/log# uri
        if apilog_key:
            apilog = ApiLog.objects.get(pk=apilog_key)
            logger.info(str(('apilog found', apilog)))
            
            if apilog.api_action == 'PATCH_LIST' and apilog.added_keys:
                keys = json.loads(apilog.added_keys)
                logger.info(str(('query for keys', keys)))
                # FIXME: on refactor; grab the key field from the resource definition
                # FIXME: on refactor, will have to deal with composite keys
                obj_query = obj_query.filter(gene_id__in=keys)
        
        return obj_query
    
    def put_list(self, request, **kwargs):
        return super(GeneResource, self).put_list(request, **kwargs) 
    
    def obj_create(self, bundle, **kwargs):
        bundle = super(GeneResource, self).obj_create(bundle, **kwargs)
        return bundle
    
class ProteinResource(ManagedModelResource):
    
    hms_gene_id = fields.CharField('gene__gene_id')
    
    class Meta:
        queryset = Protein.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= SuperUserAuthorization()        
        resource_name = 'protein'
        
        always_return_data = True
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        
    def __init__(self, **kwargs):
        super(ProteinResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # TODO: create a custom mapper for schema before others to avoid 
        # complicated regex
        return [
            url( (r"^(?P<resource_name>%s)/log/(?P<apilog>\d+)%s$" ) 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_protein_apilog_view'), 
                name="api_dispatch_protein_apilog_view"),
            url((r"^(?P<resource_name>%s)/(?P<protein_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
                ]

    def dispatch_protein_apilog_view(self, request, **kwargs):
        return self.dispatch('list', request, **kwargs)     

    def obj_get_list(self, bundle, **kwargs):
        
        # TODO: this override is implemented as a convenience to access the 
        # "batch" updated at one time.  Refactor so that the ApiLogResource
        # can redirect to -any- resource
        apilog_key = None
        if 'apilog' in kwargs:
            apilog_key = kwargs.pop('apilog')
        
        obj_query = super(ProteinResource, self).obj_get_list(bundle, **kwargs)
        
        # This is a hook to get items for the /log/log# uri
        if apilog_key:
            apilog = ApiLog.objects.get(pk=apilog_key)
            logger.info(str(('apilog found', apilog)))
            
            if apilog.api_action == 'PATCH_LIST' and apilog.added_keys:
                keys = json.loads(apilog.added_keys)
                logger.info(str(('query for keys', keys)))
                # FIXME: on refactor; grab the key field from the resource definition
                # FIXME: on refactor, will have to deal with composite keys
                obj_query = obj_query.filter(protein_id__in=keys)
        
        return obj_query
    
    def put_list(self, request, **kwargs):
        return super(ProteinResource, self).put_list(request, **kwargs) 
    
    def obj_create(self, bundle, **kwargs):
        bundle = super(ProteinResource, self).obj_create(bundle, **kwargs)
        return bundle

    def hydrate(self, bundle):
        ''' 
        Called by full_hydrate 
        sequence is obj_create->full_hydrate(hydrate, then full)->save
        
        Our custom implementation will create an auth_user for the input; so 
        there will be a reports_userprofile.user -> auth_user.
        '''
        bundle = super(ProteinResource, self).hydrate(bundle);
        
        try:
            gene = Gene.objects.get(gene_id=bundle.data.get('hms_gene_id'))
        except Exception, e:
           logger.warn(str(('could not find gene: ', bundle.data)))
           raise e
        bundle.obj.gene = gene;
        
        return bundle

class ReactionResource(ManagedModelResource):

    class Meta:
        queryset = Reaction.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        always_return_data = True
        authorization= SuperUserAuthorization()        
        resource_name = 'reaction'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        
    def __init__(self, **kwargs):
        super(ReactionResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)/(?P<r_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]
    
    def obj_create(self, bundle, **kwargs):
        logger.info(str(('===creating reaction', bundle.data)))

        return super(ReactionResource, self).obj_create(bundle, **kwargs)
