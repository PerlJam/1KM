# Django settings for 1km project
# Copy this file to settings.py and change local values for your installation
# NOTE: this settings file is not suitable for production: (turn off DEBUG, 
# logging, change SECRET_KEY).
try:
    from base_settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Base Settings not defined.  Please configure a version of
    base_settings.py for this site.'''
    del sys

# NOTE: the parent settings file defines the PROJECT_ROOT
print 'PROJECT_ROOT: ', PROJECT_ROOT, ', ' , os.path.join(PROJECT_ROOT, '..')

DEBUG = True
TEMPLATE_DEBUG = DEBUG
# TASTYPIE_FULL_DEBUG is less useful than it seems.  
# When it is used the error is easily discernable in the server logs (and not visible there otherwise), 
# but the client message is spammed with an html response with the error code buried somewhere inside.
#TASTYPIE_FULL_DEBUG = DEBUG

ADMINS = (
    ('Admin', 'admin@email.com'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', 
        'NAME': '1km',                      
        'USER': '1km',
        'PASSWORD': '',
        'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
        'PORT': '',                      # Set to empty string for default.
    }
}

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = []

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'US/Eastern'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1


# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = ''

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'x-=2gqw2@c)#_z_n0fpxj+y)v%n7&e#xpm+coo0(5^*@l___%8i*0'

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s:%(lineno)d %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(asctime)s: %(name)s:%(lineno)d: %(message)s'
        },
        'simple1': {
            'format': '%(levelname)s %(msecs)s: %(name)s:%(lineno)d: %(message)s'
        },
    },
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        },
        'logfile': {
            'level':'DEBUG',
            'class':'logging.handlers.RotatingFileHandler',
            'filename': os.path.join(PROJECT_ROOT, '..') +  "/logs/1km.log",
            'maxBytes': 5000000,
            'backupCount': 5,
            'formatter': 'simple',
        },
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        }, 
    },
    'loggers': {
        'django.request': {
            'handlers': ['console'],
            'level': 'ERROR',
            'propagate': False,
        },
        'db': {  #
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },        
        'lims': {  #
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },               
        'reports': {  #
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },
        'db.tests': {  
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },        
        'lims.tests': {  #
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports.tests': {  #
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'django': {  #
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },        
        'utils': {  # for SQL
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },        
        'tastypie': {  #
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },
    }    
}
