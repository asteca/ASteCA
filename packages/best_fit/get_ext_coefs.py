

def main(all_syst_filters, filters, colors):
    """
    """

    # Extract names of all read filters in the order in which they are stored
    # in 'isoch_list'.
    print all_syst_filters
    print filters
    print colors
    import pdb; pdb.set_trace()  # breakpoint c8739183 //
    
    all_filts = []
    for ps in all_syst_filters:
        all_filts = all_filts + list(ps[1:])
