
from astropy.io import ascii


def main(cl_region, memb_file, readda_idcol=0, readda_mpcol=-2):
    """
    Read MP values from file. Any star whose ID is not in the defined cluster
    region (cl_region) will be assigned an MP of 0.01.

    The indexes for the ID and MPs columns are hardcoded.
    """
    print('Reading membership probabilities from file.')
    # Read IDs and MPs from file.
    data = ascii.read(memb_file)
    # Read IDs as strings since that is how they are stored in 'cl_region'
    id_list, memb_probs = [str(_) for _ in data.columns[readda_idcol]],\
        data.columns[readda_mpcol]

    memb_probs_cl_region = []
    # Assign probabilities read from file according to the star's IDs.
    # Those stars not present in the list are assigned a very low value.
    for indx, star in enumerate(cl_region):
        if star[0] in id_list:
            # Index of star in file.
            i = id_list.index(star[0])
            # Assign the probability stored in file for this star.
            memb_probs_cl_region.append(memb_probs[i])
        else:
            memb_probs_cl_region.append(0.01)

    return memb_probs_cl_region
