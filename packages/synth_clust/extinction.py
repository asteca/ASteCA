

def main(
    isochrone, E_bv, E_BV_dr, Rv, ext_coefs, N_fc, DR_dist, DR_percentage,
        rand_norm, rand_unif, m_ini_idx):
    """
    Modifies color and magnitude values according to given values for the
    total absorption Av. Using this parameter instead of the E(B-V) extinction
    reduces the correlation with Rv.

    The distance modulus was applied before this function.
    N_fc is the number of filters (N_fc[0]), and colors defined (N_fc[1]).

                 |------Nf-----|         |------Nc-----|
    isochrone = [f1, f2, .., fNf, M_ini, c1, c2, .., cNc, M_b]

                 |------Nf-----|  |------Nc-----|
    ext_coefs = [cf1, ...., cfNf, cc1, ...., ccNc]

    where:
    cfX = [a, b]  ; ccX = [[a1, b1], [a2, b2]]
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

    # Av = E_bv * Rv
    # dr = E_BV_dr * Rv
    Av = E_bv
    dr = E_BV_dr

    if dr > 0.:
        Ns = isochrone.shape[-1]

        if DR_dist == 'uniform':
            Av_dr = (2 * rand_unif[:isochrone.shape[-1]] - 1) * dr
            # Av_dr = rand_unif[:Ns] * dr
        elif DR_dist == 'normal':
            # Av_dr = abs(rand_norm[:Ns]) * dr
            Av_dr = rand_norm[:Ns] * dr

        Av_dr[rand_unif[:Ns] > DR_percentage] = 0.

        # WHY WAS I USING THIS LINE??
        # Av = np.clip(Av + Av_dr, a_min=0., a_max=np.inf)
        Av = Av + Av_dr

    Nf, Nc = N_fc

    def magmove(fi):
        Ax = (ext_coefs[fi][0] + ext_coefs[fi][1] / Rv) * Av
        return Ax

    def colmove(ci):
        Ex = ((ext_coefs[ci][0][0] + ext_coefs[ci][0][1] / Rv)
              - (ext_coefs[ci][1][0] + ext_coefs[ci][1][1] / Rv)) * Av
        return Ex

    # Move filters.
    for fi in range(Nf):
        Ax_d = magmove(fi)
        isochrone[fi] += Ax_d
        # Move filters with binary data.
        isochrone[Nf + Nc + 1 + fi] += Ax_d

    # Move colors.
    for ci in range(Nf, Nf + Nc):
        Ex = colmove(ci)
        isochrone[ci] += Ex
        # Move colors with binary data.
        isochrone[m_ini_idx + Nf + ci] += Ex

    return isochrone
