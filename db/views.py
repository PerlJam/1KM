import logging
from django.contrib.auth import logout
from django.shortcuts import render
from django.http import HttpResponse
from django.http.response import HttpResponseRedirect

logger = logging.getLogger(__name__)

def main(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index_require.html', {'search': search})


def logout_page(request):
    """
    Log users out and re-direct them to the main page.
    """
    logout(request)
    return HttpResponseRedirect('/db') # todo: use system property