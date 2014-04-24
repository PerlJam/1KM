define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid'
], function($, _, Backbone, Iccbl ){
  
  var API_VERSION = 'api/v1';
  var REPORTS_API_URI = '/reports/' + API_VERSION;
  var DB_API_URI = '/db/' + API_VERSION;
  var DEBUG = true;
  

  var SchemaClass = function() {};
  SchemaClass.prototype.detailKeys = function()
  {
    var self = this;
    
    var keys = Iccbl.sortOnOrdinal(
      _.keys(self.fields), self.fields)
    var detailKeys = _(keys).filter(function(key){
      return _.has(self.fields, key) && 
          _.has(self.fields[key], 'visibility') && 
          _.contains(self.fields[key]['visibility'], 'detail');
    });
    return detailKeys;

  };  

  var AppState = Backbone.Model.extend({
    
    defaults: {

      // TODO: deprecate these variables
      // use the REPORTS_API_URI, and DB_API_URI defined below
      root_url: '/reports',  // used for the backbone history
      api_root_url: '/reports/api/v1',

      path: '',
      actualStack: [],
      
      // NOTE: in backbone, change notifications are only detected for the 
      // object identity.  For this reason, each of the following state items is
      // kept at the top level.
      current_view: 'home',
      current_resource_id: 'home',
      current_options: {},
      routing_options: {},
      current_scratch: {},
      current_details: {},
      login_information: {},

      // Define the menu data structure:
      // menu's can only be two levels deep (todo: recursive cool menu thing)
      // - the menu tree here is composed of ui_resources, defined below.
      // - the topographical structure here represents the visual structure in the app.
      menu: {
        view: 'home',
        submenus:{
          'smallmolecule': {
            view: 'list'
          },
          'reaction': {
            view: 'list'
          },

          'gene': {
            view: 'list',
            submenus: {},
          },
          'protein': {
            view: 'list',
            submenus: {},
          },
          'expressionhost': {
            view: 'list',
            submenus: {},
          },
//          'about': {
//            view: 'list',
//            submenus: {},
//          },
          'admin': {
            expanded: false,
            // TODO: menu view not implemented; would show children menu items as links in the page
            view: 'menu',  
            submenus: {
              'metahash':{
                  view: 'list',
              },
              'resource':{
                  view: 'list',
              },
              'vocabularies':{
                  view: 'list',
              },
              'apilog':{
                  view: 'list',
              },
              'user':{
                  view: 'list',
              },
              'groups':{
                  view: 'list',
              },
              'permissions':{
                  view: 'list',
              },
              'jobs':{
                view: 'list'
              }
            }
          }
        }
      },
      
      // will search ui_resources for each of the items identified as resources on the stack
      // TODO: this should also come from the API?
      ui_resources: {   
          admin: {
              title: 'Admin',
              route: '',
              view: 'HomeView',
              content_header: 'Welcome',
              description: 'Perform admin activities'
          },
          home: {
            title: 'Home',
            route: '/',
            view: 'HomeView',
            content_header: 'Welcome',
            description: 'Menu starting point'
          },
//          about: {
//            title: 'About',
//            route: 'about',
//            view: 'HomeView',
//            content_header: 'About',
//            description: 'About the 1000 Molecules project'
//          },
          jobs: {
            header_message: 'User Jobs',
            title: 'Jobs',
            route: 'jobs',
            listView: 'JobsView',
            detailView: 'JobView',
            api_resource: 'jobs',
            url_root: '/reports/api/v1',
            options: { },
            description: 'Jobs'
            
          },
          
          metahash: {
              header_message: 'Define fields for display on detail and list views',
              title: 'Field Information',
              route: 'list/metahash',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'metahash',
              url_root: '/reports/api/v1',
              options: { search: {'scope__exact': 'fields.metahash'}},
              description: 'Control field settings'
          },

          resource: {
              header_message: 'Define resources for display in the reporting interface',
              title: 'Resource Information',
              route: 'list/resource',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'resource',
              url_root: '/reports/api/v1',
              description: 'Control resource information'
          },

          vocabularies: {
              header_message: 'Define controlled vocabularies',
              title: 'Application Vocabularies',
              route: 'list/vocabularies',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'vocabularies',
              url_root: '/reports/api/v1',
              description: 'Enter controlled vocabularies'
          },

          apilog: {
              header_message: 'View change logs',
              title: 'Change logs',
              route: 'list/apilog',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'apilog',
              url_root: '/reports/api/v1',
              description: 'Change logs'
          },
          smallmolecule: {
            header_message: 'Small Molecules',
            title: 'Small Molecules',
            route: 'smallmolecule',
            list_view: 'ListView',
            detailView: 'DetailView',
            api_resource: 'smallmolecule',
            url_root: '/db/api/v1',
            description: 'Small Molecules'
          },
          reaction: {
            header_message: 'Reaction',
            title: 'Reaction',
            route: 'reaction',
            list_view: 'ListView',
            detailView: 'DetailView',
            api_resource: 'reaction',
            url_root: '/db/api/v1',
            description: 'Reaction'
          },
 
          gene: {
            header_message: 'Genes',
            title: 'Genes',
            route: 'gene',
            list_view: 'ListView',
            detailView: 'DetailView',
            api_resource: 'gene',
            url_root: '/db/api/v1',
            description: 'Genes'
          },
          protein: {
            header_message: 'Proteins',
            title: 'Proteins',
            route: 'protein',
            list_view: 'ListView',
            detailView: 'DetailView',
            api_resource: 'protein',
            url_root: '/db/api/v1',
            description: 'Proteins'
          },
          
          expressionhost: {
            header_message: 'Expression Hosts',
            title: 'Expression Hosts',
            route: 'expressionhost',
            list_view: 'ListView',
            detailView: 'DetailView',
            api_resource: 'expressionhost',
            url_root: '/db/api/v1',
            description: 'Expression Hosts'
          },
          
          user: {
            header_message: 'Users',
            title: 'Users',
            route: 'user',
            list_view: 'ListView',
            detailView: 'DetailView',
            api_resource: 'user',
            url_root: '/reports/api/v1',
            description: 'Users'
            
          }
      },
      list_defaults: {
          page: 1,
          rpp: 25,
          order: null,
          search: null,
      },
      detail_defaults: {
      }

    },

    initialize : function() {
      _.bindAll(this,'jqXHRerror', 'error');
    },
        
    start: function(callBack){
      console.log('start app_state.js');
      this.getResources(callBack);
    },

    getResources: function(callBack){
      console.log('- getResources from the server...')
      var self = this;
      // Retrieve the resource definitions from the server
      var resourceUrl = self.reportsApiUri + '/resource'
      Iccbl.getCollectionOnClient(resourceUrl, function(collection){
        
        // Store the URI for each resource.
        // TODO: store the apiUri on the resource in the server
        var schemaClass = new SchemaClass();
        var resourcesCollection = collection.toJSON();
        _.each(resourcesCollection,function(resource){
          resource.apiUri = '/' + resource.api_name + '/' + 
            self.apiVersion + '/' + resource.key;
          resource.schema = _.extend(resource.schema, schemaClass);
        });

        // 1. Create a resource hash, keyed by the resource id key
        // 2. Augment uiResources with the (api)resources
        var ui_resources = self.get('ui_resources');
        var resources = {};
        _.each(resourcesCollection, function(resource){
          resources[resource.key] = resource;
          if (_.has(ui_resources, resource.key)) {
            ui_resources[resource.key] = _.extend(
                {}, ui_resources[resource.key], resource);
          } else {
            ui_resources[resource.key] = _.extend({},resource);
          }
        });

        // 3. for "virtual" uiResources, find the underlying apiResource, 
        // and augment the uiResource with it
        _.each(_.keys(ui_resources), function(key){
          var ui_resource = ui_resources[key];
          // then augment the ui resources with their api_resource, if different
          if ( key !== ui_resource.api_resource ) {
            ui_resources[key] = _.extend(
                {}, resources[ui_resource.api_resource], ui_resource);
          }
        });
        callBack();                
      });
      
      console.log('finished getResources')
    },

    
    /**
     * Note that schema now comes from resource.schema
     * 
     * A schema: (for ResourceResource)
     * 
     * schema: {
        extraSelectorOptions: {},
        fields: {
          api_name,
          comment,
          description,
          id,
          id_attribute,
          is_restricted,
          key,
          ordinal,
          resource_uri,
          scope,
          title,
          title_attribute,
          visibility
        },
        resource_definition: {}
        },
     * 
     */
    getSchema: function(resourceId) {
      return this.getResource(resourceId).schema;
    },
    

    
    /**
     * Get a model from the server
     */
    getModel: function(resourceId, key, callBack) {
      var self = this;
      var resource = this.getResource(resourceId);
      var url = resource.apiUri + '/' + key;

      var ModelClass = Backbone.Model.extend({
        url : url,
        defaults : {}
      });
      var instance = new ModelClass();
      instance.fetch({
          success : function(model) {
            model.resource = resource;
            model.key = key;
            callBack(model);
          },
          error: this.jqXHRerror

      });
    },
    
    /**
     * A resource (ResourceResource):
     * {
          alias: "",
          api_name: "reports",
          comment: null,
          description: "Resource describes tables, queries that are available through the API",
          id: 680,
          id_attribute: [
            "key"
          ],
          is_restricted: null,
          key: "resource",
          ordinal: 0,
          resource_uri: "/reports/api/v1/resource/resource",
          schema: {},
          scope: "resource",
          title: "Resources",
          title_attribute: [
            "title"
          ],
          visibility: [
            "list",
            "detail"
          ]
        }
     *    
     */
    getResource: function(resourceId){
      var uiResources = this.get('ui_resources');
      var resources = this.get('resources');
      if(!_.has(uiResources, resourceId)) {
        if(!_.has(resources, resourceId)) {
          this.error('Unknown resource: ' + resourceId);
          throw "Unknown resource: " + resourceId;
        }
        return resources[resourceId];
      }
      return uiResources[resourceId];
    },
    
    jqXHRerror: function(model, response, options) {
      //console.log('error fetching the model: '+ model + ', response:
      // ' + JSON.stringify(response));
      var url = options.url || model.url || '';
      var msg = 'Error locating resource: ' + url;
      this.error(msg);
      
      var sep = '\n';
      if (!_.isUndefined(response.status))
          msg += sep + response.status;
      if (!_.isUndefined(response.statusText))
          msg += sep + response.statusText;
      if (!_.isEmpty(response.responseText))
          msg += sep + response.responseText;
      
      if(DEBUG) window.alert(msg);
      else console.log(msg);
      // TODO: use Bootstrap inscreen alert classed message div
    },
    
    error: function(msg){
      var msgs = this.get('messages');
      if (msgs && msgs.length > 0) {
        msgs.push(msg);
        this.set('messages', msgs);
      }else{
        this.set('messages', [msg]);
      }
    }
  });

  var appState = new AppState();
  
  appState.resources = {};   // TO be retrieved from the server  
  
  appState.apiVersion = API_VERSION;
  appState.reportsApiUri = REPORTS_API_URI;
  appState.dbApiUri = DB_API_URI;
  appState.LIST_ARGS = ['page','rpp','order','search','log'];      
  
  
  return appState;
});
