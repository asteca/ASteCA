def main(
    ext_coefs,
    rand_norm,
    rand_unif,
    DR_dist,
    DR_percentage,
    m_ini_idx,
    Av,
    dr,
    Rv,
    isochrone,
):
    """
    Modifies magnitude and color(s) according to given values for the
    total absorption Av. Using this parameter instead of the E(B-V) extinction
    reduces the correlation with Rv.

    The distance modulus was applied before this function.

    isochrone = [mag, c1, (c2), .., Mini, mag_b, c1b, .., Mini_b]
    ext_coefs = [mag_ec, c1_ec, ...]

    where:
    mag_ec = [a, b]  ; cX_ec = [[a1, b1], [a2, b2]]

    and:
    ccm_coef = a + b / Rv = ext_coefs[0] + ext_coefs[1] / Rv

    Ax = ef * Av
    m_obs = M_int + Ax + dist_mod
          = M_int + ef * R_V * E(B-V) + dist_mod
          = M_int + (a + b / Rv) * R_V * E(B-V) + dist_mod

    E(m1 - m2) = A_m1 - A_m2
               = (ef_m1 - ef_m2) * Av
               = [(a1 + b1/Rv) - (a2 + b2/Rv)] * Av
               = [(a1 - a2) + (b1 - b2)/Rv] * Av
               = (a12 + b12/Rv) * Av
               = (a12 + b12/Rv) * R_V * E(B-V)
    (m1 - m2)_obs = (m1 - m2)_int + E(m1 - m2)
    (m1 - m2)_obs = (m1 - m2)_int + (a12 + b12/Rv) * R_V * E(B-V)
    """

    if dr > 0.0:
        Ns = isochrone.shape[-1]

        if DR_dist == "uniform":
            Av_dr = (2 * rand_unif[: isochrone.shape[-1]] - 1) * dr
            # Av_dr = rand_unif[:Ns] * dr
        elif DR_dist == "normal":
            # Av_dr = abs(rand_norm[:Ns]) * dr
            Av_dr = rand_norm[:Ns] * dr

        Av_dr[rand_unif[:Ns] > DR_percentage] = 0.0
        Av = Av + Av_dr

    def colmove(ci):
        Ex = (
            (ext_coefs[ci][0][0] + ext_coefs[ci][0][1] / Rv)
            - (ext_coefs[ci][1][0] + ext_coefs[ci][1][1] / Rv)
        ) * Av
        return Ex

    # Move magnitude
    Ax_d = (ext_coefs[0][0] + ext_coefs[0][1] / Rv) * Av
    isochrone[0] += Ax_d
    # Move filters with binary data.
    if isochrone.shape[0] > m_ini_idx + 1:
        isochrone[m_ini_idx + 1] += Ax_d

    # Move colors.
    for ci in range(1, m_ini_idx):
        Ex = colmove(ci)
        isochrone[ci] += Ex
        # Move colors with binary data.
        if isochrone.shape[0] > m_ini_idx + 1:
            isochrone[m_ini_idx + 1 + ci] += Ex

    return isochrone
