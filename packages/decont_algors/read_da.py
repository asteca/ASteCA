
import numpy as np
from ..inp import data_IO


def main(cl_region, clust_name, memb_file, readda_idcol=0, readda_mpcol=-2):
    """
    Read MP values from file. Any star in the defined cluster whose ID is not
    found in the membership file will be assigned MP=0.5.

    The indexes for the ID and MPs columns are hardcoded.
    """
    print("Reading membership probabilities from file")
    # Read IDs and MPs from file.
    data = data_IO.dataRead(clust_name, memb_file, 'r')
    # Read IDs as strings since that is how they are stored in 'cl_region'
    id_list = [str(_) for _ in data.columns[readda_idcol]]
    try:
        memb_probs = data.columns[readda_mpcol]
    except IndexError:
        print("  WARNING: MPs column not found. Assigned MP=1. to all stars")
        memb_probs = np.ones(len(data))

    N_not, memb_probs_cl_region = 0, []
    # Assign probabilities read from file according to the star's IDs.
    for star in cl_region:
        if star[0] in id_list:
            # Index of star in file.
            i = id_list.index(star[0])
            # Assign the probability stored in file for this star.
            memb_probs_cl_region.append(memb_probs[i])
        else:
            # Stars not present in the list are assigned a fixed value.
            memb_probs_cl_region.append(0.5)
            N_not += 1

    if N_not > 0:
        print(("  WARNING: {} stars where not present in the membership\n" +
               "  file and were assigned MP=0.5").format(N_not))

    return memb_probs_cl_region
