

import sys


def check():
    """
    Check if all the necessary packages are installed.
    """
    # Catch in case 'pip' is not installed.
    try:
        import pip
        # check for globally installed packages.
        inst_packgs = pip.get_installed_distributions(local_only=False)
        inst_packgs_lst = ["%s" % (i.key) for i in inst_packgs]
        missing_pckg = []
        for pckg in ['numpy', 'matplotlib', 'scipy', 'astroml',
                     'scikit-learn']:
            if pckg not in inst_packgs_lst:
                missing_pckg.append(pckg)

        if missing_pckg:
            print "ERROR: the following packages are missing:\n"
            for p in missing_pckg:
                print " - {}".format(p)
            sys.exit("\nInstall with: pip install <package>\n")
    except ImportError:
        # Python versions 2.7.7 onward apparently have 'pip' included by
        # default, so this check should become obsolete.
        print("  WARNING: 'pip' is not present. Can't check for installed"
              " packages.\n")
        # Return empty list.
        inst_packgs_lst = []

    return inst_packgs_lst


if __name__ == "__main__":
    check()
