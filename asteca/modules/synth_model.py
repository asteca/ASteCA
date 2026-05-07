import numpy as np


def generate(
    params: dict,
    met_age_dict: dict,
    def_params: dict,
    m_ini_idx: int,
    theor_tracks: np.ndarray,
    ext_law: str,
    ext_coefs: list | np.ndarray,
    rand_floats: dict[str, np.ndarray],
    DR_distribution: str,
    st_dist_mass: list,
    max_mag_syn: float,
    err_dist_synth: list[np.ndarray],
    N_synth_stars: int,
    return_flag: str = "array",
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Generate a synthetic cluster.

    The synthetic cluster is generated according to the parameters given in
    the ``params`` dictionary and the already calibrated
    :py:class:`Synthetic` object.

    :param params: Dictionary containing the values for the fundamental parameters.
        The dictionary must include values for all the parameters, e.g.:
        ``params = {met: 0.0152, loga: 8.1, alpha: 0.1, beta: 1, Av: 0.2, DR: 0., Rv: 3.1, dm: 9.7}``
    :type params: dict
    :param met_age_dict: Metallicity and age grid metadata.
    :type met_age_dict: dict
    :param def_params: Dictionary of default model parameters.
    :type def_params: dict
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param theor_tracks: Processed theoretical tracks.
    :type theor_tracks: np.ndarray
    :param ext_law: Extinction law.
    :type ext_law: str
    :param ext_coefs: Extinction coefficients.
    :type ext_coefs: list | np.ndarray
    :param rand_floats: Pre-generated random values used during synthesis.
    :type rand_floats: dict[str, np.ndarray]
    :param DR_distribution: Differential reddening distribution.
    :type DR_distribution: str
    :param st_dist_mass: Sampled IMF mass distributions.
    :type st_dist_mass: list
    :param max_mag_syn: Maximum synthetic magnitude allowed.
    :type max_mag_syn: float
    :param err_dist_synth: Error distributions used to perturb synthetic photometry.
    :type err_dist_synth: list[np.ndarray]
    :param N_synth_stars: Number of synthetic stars to generate.
    :type N_synth_stars: int
    :param return_flag: Return mode: ``array``, ``isoch``, or ``isoch+array``.
    :type return_flag: str

    :return: Synthetic array for ``return_flag='array'`` or ``'isoch'``; tuple
        ``(isochrone, synthetic_array)`` for ``return_flag='isoch+array'``.
        The synthetic cluster contains the data ``[mag, c1, (c2), mass, mass_b]``,
        where ``mag`` is the magnitude, ``c1`` is the color, ``c2`` is the
        optional second color, and ``mass, mass_b`` are the masses of the single
        and secondary components of the binary systems, respectively
        (if generated). If the system is a single star, then ``mass_b==np.nan``.
    :rtype: np.ndarray | tuple[np.ndarray, np.ndarray]
    """

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    met, loga, alpha, beta, av, dr, rv, dm, ml, mh, al, ah = properModel(
        met_age_dict, def_params, params
    )

    # If (z, a) are both fixed, use the single processed isochrone
    if ml == al == mh == ah == 0:
        # The np.array() is important to avoid overwriting 'theor_tracks'
        isochrone = np.array(theor_tracks[0][0])
    else:
        # Generate a weighted average isochrone from the (z, log(age)) values in
        # the 'model'.
        isochrone = zaWAverage(
            theor_tracks,
            met_age_dict,
            m_ini_idx,
            met,
            loga,
            ml,
            mh,
            al,
            ah,
        )

    binar_flag = True
    if alpha == 0.0 and beta == 0.0:
        binar_flag = False

        # # TODO: this was not tested thoroughly (April 2025)
        # # Remove binary photometry
        # isochrone = isochrone[: self.m_ini_idx + 2]
        # # TODO: this was not tested thoroughly

    # Move theoretical isochrone using the distance modulus
    isoch_moved = move_isochrone(isochrone, binar_flag, m_ini_idx, dm)

    # Apply extinction correction
    isoch_extin = extinction(
        ext_law,
        ext_coefs,
        rand_floats["norm"][0],
        rand_floats["unif"][0],
        DR_distribution,
        m_ini_idx,
        binar_flag,
        av,
        dr,
        rv,
        isoch_moved,
    )

    # Remove isochrone stars beyond the maximum magnitude
    isoch_cut = cut_max_mag(isoch_extin, max_mag_syn)
    if isoch_cut.size == 0:
        if return_flag == "isoch+array":
            return np.array([]), np.array([])
        return np.array([])

    # Return the isochrone only
    if return_flag == "isoch":
        return isoch_cut

    # Interpolate IMF's sampled masses into the isochrone.
    isoch_mass = mass_interp(
        isoch_cut,
        m_ini_idx,
        st_dist_mass[ml][al],
        N_synth_stars,
        # binar_flag,
    )
    if isoch_mass.size == 0:
        if return_flag == "isoch+array":
            return np.array([]), np.array([])
        return np.array([])

    # import matplotlib.pyplot as plt
    # # plt.title('2000, steps 1')
    # plt.title('2000, steps 2')
    # plt.scatter(isoch_mass[1], isoch_mass[0], alpha=.25)
    # plt.scatter(isoch_mass[5], isoch_mass[4], alpha=.25)
    # plt.gca().invert_yaxis()
    # plt.show()

    # Assignment of binarity.
    isoch_binar = binarity(
        alpha,
        beta,
        binar_flag,
        m_ini_idx,
        rand_floats["unif"][1],
        isoch_mass,
    )

    # Assign errors according to errors distribution.
    synth_clust = add_errors(isoch_binar, err_dist_synth)

    # Return both the isochrone and the full synthetic cluster
    if return_flag == "isoch+array":
        return isoch_cut, synth_clust

    return synth_clust


def properModel(
    met_age_dict: dict, def_params: dict, fit_params: dict
) -> tuple[float, float, float, float, float, float, float, float, int, int, int, int]:
    """Define the 'proper' model with values for (z, a) taken from its grid,
    and filled values for those parameters that are fixed.

    :param met_age_dict: Dictionary of metallicity and age values.
    :type met_age_dict: dict
    :param def_params: Dictionary of default parameters.
    :type def_params: dict
    :param fit_params: Dictionary of fitted parameters.
    :type fit_params: dict

    :returns: A tuple containing all fundamental parameters, and the indexes of the
     (z, a) values in the grid that define the box that enclose the proper (z, a)
     values.
    :rtype: tuple[float, float, float, float, float, float, float, float, int, int, int, int]

    """
    # Combine `fit_params` with default parameter values in `def_params`, updating
    # duplicated values with those in `fit_params`
    fit_params = def_params | fit_params

    ml = mh = 0
    if len(met_age_dict["met"]) > 1:
        mh = min(
            len(met_age_dict["met"]) - 1,
            np.searchsorted(met_age_dict["met"], fit_params["met"]),
        )
        ml = mh - 1

    al = ah = 0
    if len(met_age_dict["loga"]) > 1:
        ah = min(
            len(met_age_dict["loga"]) - 1,
            np.searchsorted(met_age_dict["loga"], fit_params["loga"]),
        )
        al = ah - 1

    return (
        float(fit_params["met"]),
        float(fit_params["loga"]),
        float(fit_params["alpha"]),
        float(fit_params["beta"]),
        float(fit_params["Av"]),
        float(fit_params["DR"]),
        float(fit_params["Rv"]),
        float(fit_params["dm"]),
        int(ml),
        int(mh),
        int(al),
        int(ah),
    )


def zaWAverage(
    theor_tracks: np.ndarray,
    met_age_dict: dict,
    m_ini_idx: int,
    z_model: float,
    a_model: float,
    ml: int,
    mh: int,
    al: int,
    ah: int,
) -> np.ndarray:
    """Generate a new "weighted" isochrone from the four closest points in the
    (z, a) grid.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, (c2), M_ini, M_b, f1b,.., c1b, (c2b)]
    where:
    fX:  individual filters (mags)
    cX:  colors
    M_ini: initial mass
    M_b:  binary masses
    fXb: filters with binary data
    cXb: colors with the binary data

    It is important that the returned isochrone is a *copy* of the
    theor_tracks[mx][ax] array if no averaging is done. Otherwise the
    'theor_tracks' are modified.

    :param theor_tracks: Array of theoretical isochrones.
    :type theor_tracks: np.ndarray
    :param met_age_dict: Dictionary of metallicity and age values.
    :type met_age_dict: dict
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param z_model: Metallicity value.
    :type z_model: float
    :param a_model: Age value.
    :type a_model: float
    :param ml: Index of the lower metallicity.
    :type ml: int
    :param mh: Index of the higher metallicity.
    :type mh: int
    :param al: Index of the lower age.
    :type al: int
    :param ah: Index of the higher age.
    :type ah: int

    :returns: Weighted isochrone.
    :rtype: np.ndarray
    """
    # Distances between the (z, a) points in the 'model', and the four
    # points in the (z, a) grid that contain the model point.
    # Fast euclidean distance: https://stackoverflow.com/a/47775357/1391441

    # The four points in the (z, age) grid that define the box that contains
    # the model value (z_model, a_model)
    z1 = met_age_dict["met"][ml]
    z2 = met_age_dict["met"][mh]
    a1 = met_age_dict["loga"][al]
    a2 = met_age_dict["loga"][ah]

    # OLD Feb 2026
    # pts = np.array([(z1, a1), (z1, a2), (z2, a1), (z2, a2)])
    # a_min_b = np.array([(z_model, a_model)]) - pts
    # # Don't take the square root, it's not necessary
    # dist = np.einsum("ij,ij->i", a_min_b, a_min_b)

    # Use direct arithmetic instead of array operations
    # Calculate squared distances (skip sqrt since we only need relative values)
    dz1 = z_model - z1
    dz2 = z_model - z2
    da1 = a_model - a1
    da2 = a_model - a2
    dist = np.array(
        [
            dz1 * dz1 + da1 * da1,  # (z1, a1)
            dz1 * dz1 + da2 * da2,  # (z1, a2)
            dz2 * dz2 + da1 * da1,  # (z2, a1)
            dz2 * dz2 + da2 * da2,  # (z2, a2)
        ]
    )

    # Index to the minimum dist value
    idx_dist_min = np.argmin(dist)
    ma = ((ml, al), (ml, ah), (mh, al), (mh, ah))
    # Index to the corresponding (z, a) indexes
    m0, a0 = ma[idx_dist_min]

    # If the model has a 0. distance to the closest isochrone, return that isochrone
    if dist[idx_dist_min] == 0.0:
        return np.array(theor_tracks[m0][a0])

    # Weighted average by the (inverse) distance to the four (z, a) grid points
    inv_d = 1.0 / dist
    weights = inv_d / inv_d.sum()
    # Stack the four isochrones
    isoch_stack = np.array(
        [
            theor_tracks[ml][al],
            theor_tracks[ml][ah],
            theor_tracks[mh][al],
            theor_tracks[mh][ah],
        ]
    )
    # Weighted average: einsum is faster than manual weight * array operations
    # isochrone = (
    #     theor_tracks[ml][al] * weights[0]
    #     + theor_tracks[ml][ah] * weights[1]
    #     + theor_tracks[mh][al] * weights[2]
    #     + theor_tracks[mh][ah] * weights[3]
    # )
    isochrone = np.einsum("i,i...->...", weights, isoch_stack)

    # DO NOT average the masses or their distribution will be lost. We use the
    # mass values from the closest isochrone.
    isochrone[m_ini_idx] = theor_tracks[m0][a0][m_ini_idx]
    # Now for the secondary masses
    isochrone[m_ini_idx + 1] = theor_tracks[m0][a0][m_ini_idx + 1]

    # temp_plot(
    #     "OLD",
    #     theor_tracks,
    #     isochrone,
    #     m_ini_idx,
    #     ml,
    #     mh,
    #     al,
    #     ah,
    #     z_model,
    #     a_model,
    #     pts,
    #     weights,
    # )

    return np.ascontiguousarray(isochrone)


def temp_plot(
    title,
    theor_tracks,
    isochrone,
    m_ini_idx,
    ml,
    mh,
    al,
    ah,
    z_model,
    a_model,
    pts,
    weights,
):
    """ """
    import matplotlib.pyplot as plt

    print("ml, mh, al, ah:", ml, mh, al, ah)
    print("z_model, a_model:", z_model, a_model)
    print("(z, loga) grid:", *pts)
    print("weights:", weights)

    plt.suptitle(title)

    plt.subplot(221)
    N_pts = len(theor_tracks[ml][al][m_ini_idx])
    plt.scatter(
        # theor_tracks[ml][al][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[ml][al][0],
        c="b",
        alpha=0.5,
        label=f"mlal (w={weights[0]:.3f})",
    )
    plt.scatter(
        # theor_tracks[ml][ah][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[ml][ah][0],
        c="r",
        alpha=0.5,
        label=f"mlah (w={weights[1]:.3f})",
    )
    plt.scatter(
        # theor_tracks[mh][al][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[mh][al][0],
        c="cyan",
        alpha=0.5,
        label=f"mhal (w={weights[2]:.3f})",
    )
    plt.scatter(
        # theor_tracks[mh][ah][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[mh][ah][0],
        c="orange",
        alpha=0.5,
        label=f"mhah (w={weights[3]:.3f})",
    )
    plt.xlabel("mass")
    plt.ylabel("mag")
    # plt.scatter(isochrone[m_ini_idx], isochrone[0], marker="x", c="k")
    plt.scatter(np.arange(0, len(isochrone[0])), isochrone[0], marker="x", c="k")
    plt.legend()

    #
    plt.subplot(222)
    plt.scatter(
        # theor_tracks[ml][al][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[ml][al][1],
        c="b",
        alpha=0.5,
        label="mlal",
    )
    plt.scatter(
        # theor_tracks[ml][ah][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[ml][ah][1],
        c="r",
        alpha=0.5,
        label="mlah",
    )
    plt.scatter(
        # theor_tracks[mh][al][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[mh][al][1],
        c="cyan",
        alpha=0.5,
        label="mhal",
    )
    plt.scatter(
        # theor_tracks[mh][ah][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[mh][ah][1],
        c="orange",
        alpha=0.5,
        label="mhah",
    )
    # plt.scatter(isochrone[m_ini_idx], isochrone[1], marker="x", c="k")
    plt.scatter(np.arange(0, len(isochrone[0])), isochrone[1], marker="x", c="k")
    plt.xlabel("mass")
    plt.ylabel("color")
    plt.legend()

    # First color
    plt.subplot(223)
    plt.scatter(
        theor_tracks[ml][al][1], theor_tracks[ml][al][0], c="b", alpha=0.2, label="mlal"
    )
    plt.scatter(
        theor_tracks[ml][ah][1], theor_tracks[ml][ah][0], c="r", alpha=0.2, label="mlah"
    )
    plt.scatter(
        theor_tracks[mh][al][1],
        theor_tracks[mh][al][0],
        c="cyan",
        alpha=0.2,
        label="mhal",
    )
    plt.scatter(
        theor_tracks[mh][ah][1],
        theor_tracks[mh][ah][0],
        c="orange",
        alpha=0.2,
        label="mhah",
    )
    plt.scatter(isochrone[1], isochrone[0], c="k", marker="x", alpha=0.5, label="avrg")
    plt.gca().invert_yaxis()
    plt.legend()
    plt.xlabel("color")
    plt.ylabel("mag")

    # Plot aligned masses
    plt.subplot(224)
    mmin = theor_tracks[ml][al][m_ini_idx].min()
    mmax = theor_tracks[ml][al][m_ini_idx].max()
    plt.plot(theor_tracks[ml][al][m_ini_idx], label=f"mlal, [{mmin:.2f}, {mmax:.2f}]")
    mmin = theor_tracks[ml][ah][m_ini_idx].min()
    mmax = theor_tracks[ml][ah][m_ini_idx].max()
    plt.plot(theor_tracks[ml][ah][m_ini_idx], label=f"mlah, [{mmin:.2f}, {mmax:.2f}]")
    mmin = theor_tracks[mh][al][m_ini_idx].min()
    mmax = theor_tracks[mh][al][m_ini_idx].max()
    plt.plot(theor_tracks[mh][al][m_ini_idx], label=f"mhal, [{mmin:.2f}, {mmax:.2f}]")
    mmin = theor_tracks[mh][ah][m_ini_idx].min()
    mmax = theor_tracks[mh][ah][m_ini_idx].max()
    plt.plot(theor_tracks[mh][ah][m_ini_idx], label=f"mhah, [{mmin:.2f}, {mmax:.2f}]")
    plt.legend()
    plt.xlabel("aligned points")
    plt.ylabel("mass")

    plt.show()

    # Second color
    # plt.subplot(133)
    # plt.scatter(theor_tracks[ml][al][2], theor_tracks[ml][al][0], c='b')
    # plt.scatter(theor_tracks[ml][ah][2], theor_tracks[ml][ah][0], c='r')
    # plt.scatter(theor_tracks[mh][al][2], theor_tracks[mh][al][0], c='cyan')
    # plt.scatter(theor_tracks[mh][ah][2], theor_tracks[mh][ah][0], c='orange')
    # plt.scatter(isochrone[2], isochrone[0], c='g', ls='--')
    # plt.gca().invert_yaxis()


def move_isochrone(
    isochrone: np.ndarray, binar_flag: bool, m_ini_idx: int, dm: float
) -> np.ndarray:
    """Receives an isochrone of a given age and metallicity and modifies
    its magnitude values according to a given distance modulus.

    :param isochrone: Isochrone array.
    :type isochrone: np.ndarray
    :param binar_flag: Binary flag
    :type binar_flag: bool
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param dm: Distance modulus.
    :type dm: float

    :returns: Modified isochrone array.
    :rtype: np.ndarray
    """
    isochrone = np.ascontiguousarray(isochrone)
    # Move magnitude
    isochrone[0] += dm
    # Move binary magnitude if required
    # if isochrone.shape[0] > m_ini_idx + 1:
    if binar_flag:
        # m_ini_idx + 1 --> Secondary binary mass
        # m_ini_idx + 2 --> Binary magnitude
        isochrone[m_ini_idx + 2] += dm

    return isochrone


def extinction(
    ext_law: str,
    ext_coefs: list | np.ndarray,
    rand_norm: np.ndarray,
    rand_unif: np.ndarray,
    DR_distribution: str,
    m_ini_idx: int,
    binar_flag: bool,
    Av: float,
    dr: float,
    Rv: float,
    isochrone: np.ndarray,
) -> np.ndarray:
    """Modifies magnitude and color(s) according to given values for the
    total absorption Av. Using this parameter instead of the E(B-V) extinction
    reduces the correlation with Rv.

    The distance modulus was applied before this function.

    isochrone = [mag, c1, (c2), Mini, Mini_b, mag_b, c1b, (c2b)]

    For the CCMO model:

    ext_coefs = [mag_ec, c1_ec, ...]

    where:
    mag_ec = [a, b]  ; cX_ec = [[a1, b1], [a2, b2]]

    and:
    ccm_coef = a + b / Rv = ext_coefs[i][0] + ext_coefs[i][1] / Rv

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

    :param ext_law: Extinction law to be used.
    :type ext_law: str
    :param ext_coefs: List of extinction coefficients.
    :type ext_coefs: list | np.ndarray
    :param rand_norm: Array of random normal values.
    :type rand_norm: np.ndarray
    :param rand_unif: Array of random uniform values.
    :type rand_unif: np.ndarray
    :param DR_distribution: Distribution of the differential reddening.
    :type DR_distribution: str
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param binar_flag: Flag to indicate if binarity is being used.
    :type binar_flag: bool
    :param Av: Total absorption.
    :type Av: float
    :param dr: Differential reddening.
    :type dr: float
    :param Rv: Total-to-selective extinction ratio.
    :type Rv: float
    :param isochrone: Isochrone array.
    :type isochrone: np.ndarray

    :raises ValueError: If the extinction law is not recognized.

    :returns: Modified isochrone array.
    :rtype: np.ndarray
    """
    isochrone = np.ascontiguousarray(isochrone)

    Av_dr = Av
    if dr > 0.0:
        Ns = isochrone.shape[-1]

        if DR_distribution == "uniform":
            # Transform uniform random values from [0, 1] to [-1, 1], then scale by dr
            dr_arr = (2.0 * rand_unif[:Ns] - 1.0) * dr
        elif DR_distribution == "normal":
            dr_arr = rand_norm[:Ns] * dr
        else:
            dr_arr = 0

        # In place in case I ever want to implement the percentage of stars affected.
        # Without this, all stars are affected by the DR.
        # dr_arr[rand_unif[:Ns] > DR_percentage] = 0.0

        # Clip at 0
        Av_dr = np.clip(Av + dr_arr, a_min=0, a_max=None)

    if ext_law == "CCMO":
        # Magnitude
        ec_mag = ext_coefs[0][0] + ext_coefs[0][1] / Rv
        # First color
        ec_col1 = (ext_coefs[1][0][0] + ext_coefs[1][0][1] / Rv) - (
            ext_coefs[1][1][0] + ext_coefs[1][1][1] / Rv
        )
    elif ext_law == "GAIADR3":
        # If this model is used the first color is always expected to be BP-RP
        # BP_RP = isochrone[1]
        ec_mag, ec_col1 = gaiadr3_extlaw(isochrone[1], Av_dr)
    else:
        raise ValueError(f"Unknown extinction law: {ext_law}")

    Ax = ec_mag * Av_dr
    isochrone[0] += Ax  # Magnitude
    Ex1 = ec_col1 * Av_dr
    isochrone[1] += Ex1  # First color

    # Move binary data.
    if binar_flag:
        isochrone[m_ini_idx + 2] += Ax  # Magnitude (binary system)
        isochrone[m_ini_idx + 3] += Ex1  # First color (binary system)

    # Second color
    if len(ext_coefs) > 2:
        ec_col2 = (ext_coefs[2][0][0] + ext_coefs[2][0][1] / Rv) - (
            ext_coefs[2][1][0] + ext_coefs[2][1][1] / Rv
        )
        Ex2 = ec_col2 * Av_dr
        isochrone[2] += Ex2
        # Move color with binary data.
        if binar_flag:
            isochrone[m_ini_idx + 4] += Ex2  # Second color (binary system)

    return isochrone


def gaiadr3_extlaw(
    X_: np.ndarray, Av_dr: float | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    The 'coeffs' values are the main sequence values taken from:
    https://www.cosmos.esa.int/web/gaia/edr3-extinction-law

    The order of the coefficients is:
    Intercept   X   X2  X3  A   A2  A3  XA  AX2 XA2

    :param X_: Array of BP-RP colors.
    :type X_: np.ndarray
    :param Av_dr: Total absorption (eventually containing differential reddening).
    :type Av_dr: float | np.ndarray

    :returns: Extinction coefficients for G and BP-RP.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    coeffs = {
        "G": (
            0.995969721536602,
            -0.159726460302015,
            0.0122380738156057,
            0.00090726555099859,
            -0.0377160263914123,
            0.00151347495244888,
            -2.52364537395142e-05,
            0.0114522658102451,
            -0.000936914989014318,
            -0.000260296774134201,
        ),
        "BP": (
            1.15363197483424,
            -0.0814012991657388,
            -0.036013023976704,
            0.0192143585568966,
            -0.022397548243016,
            0.000840562680547171,
            -1.31018008013549e-05,
            0.00660124080271006,
            -0.000882247501989453,
            -0.000111215755291684,
        ),
        "RP": (
            0.66320787941067,
            -0.0179847164933981,
            0.000493769449961458,
            -0.00267994405695751,
            -0.00651422146709376,
            3.30179903473159e-05,
            1.57894227641527e-06,
            -7.9800898337247e-05,
            0.000255679812110045,
            1.10476584967393e-05,
        ),
    }

    X_2 = X_ * X_
    X_3 = X_2 * X_
    Av_2 = Av_dr * Av_dr
    Av_3 = Av_2 * Av_dr

    def ext_coeff(k):
        c = coeffs[k]
        ay = (
            c[0]  # a1
            + c[1] * X_  # a2*X
            + c[2] * X_2  # a3*X^2
            + c[3] * X_3  # a4*X^3
            + c[4] * Av_dr  # a5*A0
            + c[5] * Av_2  # a6*A0^2
            + c[6] * Av_3  # a7*A0^3
            + c[7] * X_ * Av_dr  # a8*A0*X
            + c[8] * X_2 * Av_dr  # a9*A0*X^2
            + c[9] * X_ * Av_2  # a10*X*A0^2
        )
        return ay

    ec_G = ext_coeff("G")
    ec_BPRP = ext_coeff("BP") - ext_coeff("RP")

    return ec_G, ec_BPRP


def cut_max_mag(isoch_moved: np.ndarray, max_mag_syn: float) -> np.ndarray:
    """Remove stars from isochrone with magnitude values larger that the maximum
    observed value.

    :param isoch_moved: Isochrone array.
    :type isoch_moved: np.ndarray
    :param max_mag_syn: Maximum magnitude value.
    :type max_mag_syn: float

    :returns: Isochrone array with stars above the maximum magnitude removed.
    :rtype: np.ndarray
    """
    # Discard stars in isochrone beyond max_mag_syn limit.
    return isoch_moved[:, isoch_moved[0] < max_mag_syn]


def mass_interp(
    isoch_cut: np.ndarray,
    m_ini_idx: int,
    st_dist_mass: np.ndarray,
    N_synth_stars: int,
    # binar_flag: bool,
) -> np.ndarray:
    """For each mass in the sampled IMF mass distribution, interpolate its value
    (and those of all the sub-arrays in 'isoch_cut') into the isochrone.

    Masses that fall outside of the isochrone's mass range have been previously
    rejected.

    :param isoch_cut: Isochrone array.
    :type isoch_cut: np.ndarray
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param st_dist_mass: Array of sampled masses.
    :type st_dist_mass: np.ndarray
    :param N_synth_stars: Number of observed stars.
    :type N_synth_stars: int

    :returns: Interpolated isochrone array.
    :rtype: np.ndarray
    """
    isoch_cut = np.ascontiguousarray(isoch_cut)

    # ==IMPORTANT==
    # Assumes `mass_ini=isoch_cut[m_ini_idx]` is sorted min to max. This is checked
    # when loading the isochrones in isochrones_priv.interp_isochrones()
    # ==IMPORTANT==
    mass_ini = isoch_cut[m_ini_idx]

    # Filter masses in the IMF sampling that are outside of the mass
    # range given by 'isoch_cut' (st_dist_mass[0]: sampled masses from IMF)

    # # msk_min = (st_dist_mass[0] >= mass_ini.min())
    # # msk_max = (st_dist_mass[0] <= mass_ini.max())
    # # (~msk_min).sum(): stars lost below the minimum mass (photometric)
    # # (~msk_max).sum(): stars lost above the maximum mass (evolutionary)
    # msk_m = (st_dist_mass >= mass_ini.min()) & (st_dist_mass <= mass_ini.max())

    # Use direct indexing instead of min/max calls
    # Combined boolean mask in single operation
    msk_m = (st_dist_mass >= mass_ini[0]) & (st_dist_mass <= mass_ini[-1])
    # Extract up to 'N_synth_stars' masses sampled from an IMF, within the mass range
    mass_dist = st_dist_mass[msk_m][:N_synth_stars]

    if mass_dist.size == 0:
        return np.array([])

    mass_dist = np.ascontiguousarray(mass_dist, dtype=mass_ini.dtype)
    # Interpolate the sampled stars (masses) into the isochrone
    isoch_mass = interp_mass_isoch(isoch_cut, mass_ini, mass_dist)

    return isoch_mass


def interp_mass_isoch(
    isoch_cut: np.ndarray,
    mass_ini: np.ndarray,
    mass_dist: np.ndarray,
    # m_ini_idx: int,
    # binar_flag: bool,
) -> np.ndarray:
    """Find where in the original data, the values to interpolate would be inserted.

    NOTE: I already tried to speed this block up using numba (@jit(nopython=True))
    but it does not help. The code runs even slower.

    :param isoch_cut: Isochrone array.
    :type isoch_cut: np.ndarray
    :param mass_ini: Array of initial masses.
    :type mass_ini: np.ndarray
    :param mass_dist: Array of sampled masses.
    :type mass_dist: np.ndarray

    :returns: Interpolated isochrone array.
    :rtype: np.ndarray
    """

    # *****************************************************************************
    # This is a stripped down version of scipy.inter1d (kind='linear') which
    # is now legacy code (https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html)
    # Scipy now recommends using np.interp for simple linear interpolation but it
    # is quite slower.
    #
    # Note: If x_new[n] == x[m], then m is returned by searchsorted.
    x_new_indices = np.searchsorted(mass_ini, mass_dist)
    #
    # Calculate the slope of regions that each x_new value falls in.
    lo = x_new_indices - 1
    x_lo = mass_ini[lo]
    x_hi = mass_ini[x_new_indices]
    # y_lo = isoch_cut[:, lo]
    # y_hi = isoch_cut[:, x_new_indices]
    y_lo = np.take(isoch_cut, lo, axis=1)
    y_hi = np.take(isoch_cut, x_new_indices, axis=1)

    # Compute slope and result in fewer steps
    x_diff = mass_dist - x_lo
    # Calculate the actual value for each entry in x_new.
    # slope = (y_hi - y_lo) / (x_hi - x_lo)
    isoch_mass = y_lo + (y_hi - y_lo) * (x_diff / (x_hi - x_lo))
    # *****************************************************************************

    # *****************************************************************************
    # This is slower for more than 5 dimensions, and *very* slightly faster for lower
    # dimensions
    # isoch_mass = np.empty([isoch_cut.shape[0], mass_dist.size])
    # for i, arr in enumerate(isoch_cut):
    #     isoch_mass[i] = np.interp(mass_dist, mass_ini, arr)

    # # Same performance as above
    # isoch_mass = []
    # for arr in isoch_cut:
    #     isoch_mass.append(np.interp(mass_dist, mass_ini, arr))
    # isoch_mass = np.array(isoch_mass)
    # *****************************************************************************

    # *****************************************************************************
    # This is about 2x faster than interp1d() but it comes at the cost of a more
    # coarse distribution of stars throughout the isochrone. To fix this, we have to
    # apply a random noise to the magnitude(s) proportional to the percentage that
    # the masses sampled differ from those in the isochrone. This lowers the
    # performance to an increase of ~22%
    #
    # idx = np.searchsorted(mass_ini, mass_dist)
    # isoch_mass = isoch_cut[:, idx]

    # Random noise applied to magnitudes. Does not really work that good
    # m_perc = (isoch_mass[m_ini_idx] - mass_dist) / mass_dist
    # # Random noise to magnitude
    # isoch_mass[0] += isoch_mass[0] * m_perc
    # if binar_flag:
    #     # Random noise to binary magnitude
    #     isoch_mass[m_ini_idx + 2] += isoch_mass[m_ini_idx + 2] * m_perc
    # *****************************************************************************

    # import matplotlib.pyplot as plt
    # # plt.subplot(121)
    # # plt.title("interp1d")
    # plt.scatter(isoch_mass0[1], isoch_mass0[0], alpha=0.25, label="interp1d")
    # # plt.gca().invert_yaxis()
    # # plt.subplot(122)
    # # plt.title("searchsorted")
    # plt.scatter(isoch_mass[1], isoch_mass[0], alpha=0.25, label="searchsorted")
    # plt.gca().invert_yaxis()
    # plt.legend()
    # plt.show()
    # # breakpoint()

    return isoch_mass


def binarity(
    alpha: float,
    beta: float,
    binar_flag: bool,
    m_ini_idx: int,
    rand_unif_vals: np.ndarray,
    isoch_mass: np.ndarray,
) -> np.ndarray:
    """Select a fraction of stars to be binaries, given a chosen method.

    :param alpha: Alpha parameter for the binary fraction.
    :type alpha: float
    :param beta: Beta parameter for the binary fraction.
    :type beta: float
    :param binar_flag: Flag to indicate if binarity is being used.
    :type binar_flag: bool
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param rand_unif_vals: Array of random uniform values.
    :type rand_unif_vals: np.ndarray
    :param isoch_mass: Isochrone array.
    :type isoch_mass: np.ndarray

    :returns: Isochrone array with binary stars.
    :rtype: np.ndarray
    """
    # No binarity process defined
    if binar_flag is False:
        # Update the binary systems' masses so that the secondary masses for
        # SINGLE systems are identified with a 'nan' value.
        isoch_mass[m_ini_idx + 1] = np.nan
        return isoch_mass[: m_ini_idx + 2]

    Ns = isoch_mass.shape[-1]

    # Distribution of probability of binarity (multiplicity fraction) versus
    # primary masses.

    # Offner et al. (2022); Fig 1 (left), Table 1
    # b_p = np.clip(alpha + beta * np.log(mass), a_min=0, a_max=1)
    # b_p = np.clip(alpha + beta * np.arctan(mass), a_min=0, a_max=1)
    b_p = np.clip(alpha + beta / (1.0 + 1.4 / isoch_mass[m_ini_idx]), 0.0, 1.0)

    # Stars (masses) with the largest binary probabilities are selected
    # proportional to their probability
    bin_indxs = b_p > rand_unif_vals[:Ns]

    # Index of the binary magnitude: mag_binar = m_ini_idx + 1
    # Update array with new values of magnitudes and colors for the binary
    # systems. This does not change the primary mass values, just the magnitude and
    # color(s).
    if bin_indxs.any():
        for i in range(m_ini_idx):
            # isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 2 + i][bin_indxs]
            # Use np.copyto with where parameter (faster than boolean indexing)
            np.copyto(isoch_mass[i], isoch_mass[m_ini_idx + 2 + i], where=bin_indxs)

    # Update the binary systems' masses so that the secondary masses for
    # SINGLE systems are identified with a 'nan' value.
    isoch_mass[m_ini_idx + 1, ~bin_indxs] = np.nan

    # Return [mag, c1, (c2), mass, mass_b]
    return isoch_mass[: m_ini_idx + 2]


def add_errors(isoch_binar: np.ndarray, err_dist: list[np.ndarray]) -> np.ndarray:
    """Add random synthetic uncertainties to the magnitude and color(s)

    :param isoch_binar: Isochrone array.
    :type isoch_binar: np.ndarray
    :param err_dist: List of arrays of error distributions.
    :type err_dist: list[np.ndarray]

    :returns: Isochrone array with added errors.
    :rtype: np.ndarray
    """
    if not err_dist:
        return isoch_binar

    N = isoch_binar.shape[1]
    mag_sort = np.argsort(-isoch_binar[0])
    for i, sigma in enumerate(err_dist):
        isoch_binar[i, mag_sort] += sigma[:N]

    return isoch_binar
