from __future__ import unicode_literals

import logging
from tastypie.utils import timezone
from tastypie.exceptions import BadRequest
from tastypie.utils.urls import trailing_slash
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication,\
    MultiAuthentication
from tastypie.constants import ALL_WITH_RELATIONS
from tastypie import fields

from django.db import models
from django.conf.urls import url

from lims.api import CursorSerializer, LimsSerializer, SmallMoleculeSerializer
from reports.models import MetaHash, Vocabularies, ApiLog
from reports.api import ManagedModelResource, ManagedResource, ApiLogResource
from db.models import SmallMolecule, Reaction


import logging

logger = logging.getLogger(__name__)

class SmallMoleculeResource(ManagedModelResource):

    class Meta:
        queryset = SmallMolecule.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= Authorization()        
        resource_name = 'smallmolecule'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        
    def __init__(self, **kwargs):
        super(SmallMoleculeResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)/(?P<sm_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]
    
    def obj_create(self, bundle, **kwargs):
        logger.info(str(('===creating smallmolecule', bundle.data)))

        return super(SmallMoleculeResource, self).obj_create(bundle, **kwargs)


class ReactionResource(ManagedModelResource):

    class Meta:
        queryset = Reaction.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= Authorization()        
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