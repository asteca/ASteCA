def main(isochrone, m_ini_idx, dm):
    """
    Receives an isochrone of a given age and metallicity and modifies
    its magnitude values according to a given distance modulus.
    """
    # Move magnitude
    isochrone[0] += dm
    # Move binary magnitude if it exists
    if isochrone.shape[0] > m_ini_idx + 1:
        isochrone[m_ini_idx + 1] += dm

    return isochrone
