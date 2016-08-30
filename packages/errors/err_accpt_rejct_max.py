

def main(cld, er_params, **kwargs):
    """
    Accept stars with photometric errors < e_max both in mag and in color.
    """

    # Unpack params.
    e_max = er_params[1]
    # USE FIRST MAGNITUDE AND COLOR READ.
    em, ec = cld['em'][0], cld['ec'][0]

    # Initialize empty list to hold accepted stars' indexes.
    acpt_indx = []
    # Iterate through all stars
    for st_ind in range(len(em)):

        # Reject stars with at least one error >= e_max.
        if em[st_ind] >= e_max or ec[st_ind] >= e_max:
            pass
        else:
            # Accept star.
            acpt_indx.append(st_ind)

    # Pass empty list.
    err_plot = []
    return acpt_indx, err_plot
