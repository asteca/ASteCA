# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 12:00:00 2015

@author: gabriel
"""

from functions import __version__
import urllib2


def updater():
    '''
    Checks if a new version of the code is available for download.
    '''

    try:
        # Get latest version number. Wait 3 seconds and break out if
        # there's no response from the server.
        f = urllib2.urlopen("https://raw.githubusercontent.com/asteca/"
        "asteca/master/functions/__init__.py", timeout=3)
        s = f.read().split('"')

        if s[1] != __version__:
            print "*******************************************"
            print "              IMPORTANT\n"
            print "There is a new version of ASteCA available!"
            print "   Get the latest version '{}' from:\n".format(s[1][1:])
            print "       http://asteca.github.io/"
            print "*******************************************\n"
    except:
        #import traceback
        #print traceback.format_exc()
        pass
