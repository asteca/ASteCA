
from packages._version import __version__
# In place for #243
import sys
if sys.version_info[0] == 2:
    from urllib2 import urlopen as urlopen
    from urllib2 import URLError as URLError
else:
    from urllib.request import urlopen as urlopen
    from urllib.error import URLError as URLError


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
            f = urlopen("https://raw.githubusercontent.com/asteca/"
                        "AsteCA/master/packages/_version.py",
                        timeout=t_out)
            if sys.version_info[0] == 2:
                s = f.read().split('"')
            else:
                s = f.read().decode('utf-8').split('"')

            if s[1] != __version__:
                print("*******************************************")
                print("              IMPORTANT\n")
                print("There is a new version of ASteCA available!")
                print("  Get the latest version '{}' from:\n".format(s[1][1:]))
                print("       http://asteca.github.io/")
                print("*******************************************\n")
            else:
                print("You are running the latest version of ASteCA.\n")
        except URLError:
            print("  WARNING: could not check for code updates.\n"
                  "  Connection timed out after {} sec.\n".format(t_out))
