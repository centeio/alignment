import os
from .base import *  # noqa

DEBUG = False


INTERNAL_IPS = ["127.0.0.1"]

SECRET_KEY = "secret"

# DATABASE SETTINGS
# https://docs.djangoproject.com/en/1.10/ref/settings/#databases
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'development.sqlite3',
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    },
}

CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.dummy.DummyCache"
    }
}

EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"



# DATABASE SETTINGS
# https://docs.djangoproject.com/en/1.10/ref/settings/#databases
# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.',
#         'NAME': '',
#         'USER': '',
#         'PASSWORD': '',
#         'HOST': '',
#         'PORT': '',
#     },
# }

# IMPORTANT!:
# You must keep this secret, you can store it in an
# environment variable and set it with:
# export SECRET_KEY="phil-dunphy98!-bananas12"
# https://docs.djangoproject.com/en/1.10/howto/deployment/checklist/#secret-key

# WSGI SETTINGS
# https://docs.djangoproject.com/en/1.10/ref/settings/#wsgi-application
WSGI_APPLICATION = 'alignment_awesome_tool.wsgi.application'

# NOTIFICATIONS
# A tuple that lists people who get code error notifications.
# https://docs.djangoproject.com/en/1.10/ref/settings/#admins
ADMINS = (
         ('Higor Anjos', 'higor.araujo.anjos@gmail.com'),
)
MANAGERS = ADMINS

COMPRESS_ENABLED = os.environ.get('COMPRESS_ENABLED', False) 

# DJANGO-COMPRESSOR SETTINGS
# STATICFILES_FINDERS = STATICFILES_FINDERS + (
#     'compressor.finders.CompressorFinder',
# )
ALLOWED_HOSTS = ['localhost', 'www.saat.herokuapp.com']

try:
    from local_settings import * # noqa
except ImportError:
    pass
