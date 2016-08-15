
import sys
import traceback
from packages.inp import isoch_params


def check_get(pd):
    """
    Check that all metallicity files needed are in place. To save time, we
    store the data and pass it.
    """
    # Read metallicity files.
    try:
        # Store all isochrones in all the metallicity files in isoch_list.
        ip_list = isoch_params.main(**pd)
    except:
        print traceback.format_exc()
        sys.exit("ERROR: unknown error reading metallicity files.")

    return ip_list
