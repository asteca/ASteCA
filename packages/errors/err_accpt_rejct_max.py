

def main(cld, err_max):
    """
    Accept stars with photometric errors < e_max both in mag and in color.
    """
    # USE FIRST MAGNITUDE AND COLOR READ.
    em, ec = cld['em'][0], cld['ec'][0]

    # Initialize empty list to hold accepted stars' indexes.
    acpt_indx = []
    # Iterate through all stars
    for st_ind in range(len(em)):

        # Reject stars with at least one error >= err_max.
        if em[st_ind] >= err_max or ec[st_ind] >= err_max:
            pass
        else:
            # Accept star.
            acpt_indx.append(st_ind)

    return acpt_indx
