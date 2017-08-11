

def main(cl_region, fixedda_port):
    """
    Assign MP values to stars within the cluster region, divided into
    4 sections from brightest to dimmest.
    """

    # Cluster's sequence using the main magnitude.
    main_mag = list(zip(*zip(*cl_region)[3])[0])
    mi, mf = min(main_mag), max(main_mag)
    # Magnitude range
    deltam = mf - mi
    # Magnitude range of the bottom part of the sequence, divided by three.
    deltam2 = (1. - fixedda_port) * deltam / 3.
    # First break point (end of top sequence part)
    m1 = mi + fixedda_port * deltam
    # Second and third sequence points.
    m2 = m1 + deltam2
    m3 = m2 + deltam2

    memb_probs_cl_region = []
    for star in cl_region:
        # Main magnitude
        mmag = star[3][0]
        if mmag < m1:
            memb_probs_cl_region.append(1.)
        elif m1 <= mmag < m2:
            memb_probs_cl_region.append(.75)
        elif m2 <= mmag < m3:
            memb_probs_cl_region.append(.5)
        elif m3 <= mmag:
            memb_probs_cl_region.append(.25)

    return memb_probs_cl_region
