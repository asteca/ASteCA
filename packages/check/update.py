
import urllib2
from packages._version import __version__


def check(up_flag, **kwargs):
    '''
    Checks if a new version of the code is available for download.
    '''
    if up_flag:
        print("Checking for updates...")
        t_out = 3.
        try:
            # Get latest version number. Wait 3 seconds and break out if
            # there's no response from the server.
            f = urllib2.urlopen("https://raw.githubusercontent.com/asteca/"
                                "AsteCA/master/packages/_version.py",
                                timeout=t_out)
            s = f.read().split('"')

            if s[1] != __version__:
                print "*******************************************"
                print "              IMPORTANT\n"
                print "There is a new version of ASteCA available!"
                print "   Get the latest version '{}' from:\n".format(s[1][1:])
                print "       http://asteca.github.io/"
                print "*******************************************\n"
            else:
                print("You are running the latest version of ASteCA.\n")
        except urllib2.URLError:
            print("  WARNING: could not check for code updates.\n"
                  "  Connection timed out after {} sec.\n".format(t_out))
