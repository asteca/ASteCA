

def main(isochrone, N_fc, dm):
    """
    Receives an isochrone of a given age and metallicity and modifies
    its magnitude values according to a given distance modulus.
    N_fc is the number of filters (N_fc[0]), and colors defined (N_fc[1]).

    The maximum observed magnitude cut is applied here to avoid carrying a
    lot of stars in the following functions.
    """
    Nf, Nc = N_fc

    # Move filters.
    for fi in range(Nf):
        isochrone[fi] += dm
        # Move filters with binary data.
        isochrone[Nf + Nc + 1 + fi] += dm

    return isochrone
