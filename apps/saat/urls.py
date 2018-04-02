"""urlconf for the base application"""

from django.conf.urls import url

from . import views


urlpatterns = [
    url(r'^$', views.home, name='home'),
    url(r'^globallinear', views.globallinear, name='globallinear'),
    url(r'^globalmultiples', views.globalmultiples, name='globalmultiples'),
    url(r'^local', views.local, name='local'),
    url(r'^result', views.result, name='result'),
    
]
