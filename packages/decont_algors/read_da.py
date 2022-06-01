
import numpy as np
from ..inp import data_IO


def main(cl_region, data_file, id_ids, da_algor):
    """
    Read MP values from file. Any star in the defined cluster whose ID is not
    found in the membership file will be assigned MP=0.5.
    """
    print("Reading membership probabilities from file")
    # Read input file.
    data = data_IO.dataRead(None, data_file, 'r')

    # Read IDs as strings since that is how they are stored in 'cl_region'
    id_list = list(np.array(data[id_ids], dtype=str))
    memb_probs = data[da_algor]

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
        print(("  WARNING: {} stars where not present in the membership\n"
               + "  file and were assigned MP=0.5").format(N_not))

    return memb_probs_cl_region
