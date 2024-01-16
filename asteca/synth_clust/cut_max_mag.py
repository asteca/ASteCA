def main(isoch_moved, max_mag_syn):
    """
    Remove stars from isochrone with magnitude values larger that the maximum
    observed value.
    """
    # Discard stars in isochrone beyond max_mag_syn limit.
    msk = isoch_moved[0] < max_mag_syn
    return isoch_moved[:, msk], (~msk).sum()
