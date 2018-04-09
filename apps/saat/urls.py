"""urlconf for the base application"""

from django.conf.urls import url

from . import views


urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^global', views.home, name='global'),
    url(r'^local', views.local, name='local'),
    
]
