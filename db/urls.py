from tastypie.api import Api
from django.conf.urls import patterns, include, url
from django.contrib import admin

from db import views
from db.api import SmallMoleculeResource, ReactionResource

admin.autodiscover()

v1_api = Api(api_name='v1')
v1_api.register(SmallMoleculeResource())
v1_api.register(ReactionResource())

urlpatterns = patterns('',
    url(r'^$', views.main, name="home"),
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', include(admin.site.urls)),

    # Login / logout.
    # Note: login.html is actually served by the reports project:
    # that is, reports/templates/login.html version; 
    # this is done because at this time only the reports project has
    # all of the necessary css and javascript installed
    (r'^accounts/login/$', 'django.contrib.auth.views.login', {'template_name': 'login.html'}),
    url(r'^accounts/logout/$', views.logout_page, name='logout'),

    (r'^api/', include(v1_api.urls)),

    url(r'^reports/', include('reports.urls')),
)

