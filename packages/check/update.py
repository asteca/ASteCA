
import random
from packages._version import __version__
from urllib.request import urlopen as urlopen
from urllib.error import URLError as URLError


def check():
    """
    Checks if a new version of the code is available for download.
    """
    # Check every fourth run approximately.
    if random.choice([True] + [False] * 3):
        t_out = 3.
        try:
            # Get latest version number. Wait 3 seconds and break out if
            # there's no response from the server.
            f = urlopen("https://raw.githubusercontent.com/asteca/"
                        "AsteCA/master/packages/_version.py",
                        timeout=t_out)
            s = f.read().decode('utf-8').split('"')

            if s[1][1:] != __version__.strip('-dev'):
                print("*******************************************")
                print("              IMPORTANT\n")
                print("There is a new version of ASteCA available")
                print("    Upgrading is STRONGLY recommended!\n")
                print("  Get the latest version '{}' from:".format(s[1][1:]))
                print("       http://asteca.github.io/")
                print("*******************************************\n")

        except URLError:
            print("  WARNING: could not check for code updates.\n"
                  "  Connection timed out after {} sec.\n".format(t_out))
