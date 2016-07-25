
import traceback
import urllib2
from packages._version import __version__


def check():
    '''
    Checks if a new version of the code is available for download.
    '''
    try:
        # Get latest version number. Wait 3 seconds and break out if
        # there's no response from the server.
        f = urllib2.urlopen("https://raw.githubusercontent.com/asteca/"
                            "AsteCA/master/packages/_version.py", timeout=3)
        s = f.read().split('"')

        if s[1] != __version__:
            print "*******************************************"
            print "              IMPORTANT\n"
            print "There is a new version of ASteCA available!"
            print "   Get the latest version '{}' from:\n".format(s[1][1:])
            print "       http://asteca.github.io/"
            print "*******************************************\n"
    except:
        print("  WARNING: could not check for code updates.\n")
        print traceback.format_exc()
