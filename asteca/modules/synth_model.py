import numpy as np

from . import HAS_NUMBA as _HAS_NUMBA
from . import njit

# The 'coeffs' values are the main sequence values taken from:
# https://www.cosmos.esa.int/web/gaia/edr3-extinction-law
_G_COEFFS = (
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
)
# "BP" minus "RP" coefficients
_BPRP_COEFFS = (
    0.49042409542357,
    -0.0634165826723407,
    -0.036506793426665456,
    0.02189430261385411,
    -0.01588332677592224,
    0.0008075446901998552,
    -1.468074307777017e-05,
    0.006681041701047307,
    -0.001137927314099498,
    -0.0001222634137884233,
)
# Range of intrinsic BP-RP colors for which the coefficients above were fitted. The
# coefficients should not be extrapolated beyond it (the polynomials diverge), so
# the colors used to evaluate them are clamped to this range
_X_MIN, _X_MAX = -0.06, 2.5


def main(N_repeat: int = 25_000):
    """
    Main function used to test the performance of the `generate` function. It
    generates a synthetic cluster N_repeat times, and profiles the execution time
    of the function.

    :param N_repeat: Number of times to repeat the generation of the synthetic cluster.
    :type N_repeat: int
    """
    import asteca

    print(asteca.__version__)
    import pandas as pd
    from pyinstrument import Profiler

    # Isochrones parameters
    isochs = asteca.Isochrones(
        model="parsec",
        isochs_path="../../docs/_static/isochrones/parsec/gaia",
        mag="Gmag",
        color=("G_BPmag", "G_RPmag"),
        verbose=3,
    )

    # Synthetic clusters parameters
    synthcl = asteca.Synthetic(
        isochs,
        ext_law="GAIADR3",
        max_mass=40000,
        seed=42519,
        verbose=3,
    )

    obs_df = pd.read_csv("../../docs/_static/cluster.csv")
    my_cluster = asteca.Cluster(
        ra=np.array(obs_df["RA_ICRS"]),
        dec=np.array(obs_df["DE_ICRS"]),
        mag=np.array(obs_df["Gmag"]),
        e_mag=np.array(obs_df["e_Gmag"]),
        color=np.array(obs_df["BP-RP"]),
        e_color=np.array(obs_df["e_BP-RP"]),
        verbose=2,
    )
    synthcl.calibrate(my_cluster)

    rng = np.random.default_rng(12345)
    met_min, met_max = synthcl.met_age_dict["met"][0], synthcl.met_age_dict["met"][-1]
    loga_min, loga_max = (
        synthcl.met_age_dict["loga"][0],
        synthcl.met_age_dict["loga"][-1],
    )

    # Start profiler
    profiler = Profiler()
    profiler.start()
    for _ in range(N_repeat):
        # Random parameter values within reasonable ranges
        _av = rng.uniform(0.0, 3.0)
        params = {
            "met": rng.uniform(met_min, met_max),
            "loga": rng.uniform(loga_min, loga_max),
            "alpha": rng.uniform(0.0, 0.2),
            "beta": rng.uniform(0.5, 1.5),
            "Av": _av,
            "DR": rng.uniform(0.0, _av),
            "Rv": rng.uniform(2.5, 4.0),
            "dm": rng.uniform(8.0, 14.0),
        }
        _ = get_model(
            params,
            synthcl.met_age_dict,
            synthcl.met_tol,
            synthcl.age_tol,
            synthcl.def_params,
            synthcl.m_ini_idx,
            synthcl.theor_tracks,
            synthcl.ext_law,
            synthcl.ext_coefs,
            synthcl.rand_floats,
            synthcl.st_dist_mass,
            synthcl.mass_idx,
            synthcl.max_mag_syn_obs,
            synthcl.err_dist_obs,
            synthcl.N_stars_obs,
        )
    profiler.stop()
    profiler.open_in_browser()


def get_model(
    params: dict,
    met_age_dict: dict,
    met_tol: float,
    age_tol: float,
    def_params: dict,
    m_ini_idx: int,
    theor_tracks: np.ndarray,
    ext_law: str,
    ext_coefs: list | np.ndarray,
    rand_floats: dict[str, np.ndarray],
    st_dist_mass: list,
    mass_idx_lookup: dict,
    max_mag_syn: float,
    err_dist_synth: list[np.ndarray],
    N_synth_stars: int,
    show_errors: bool = False,
    isoch_flag: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a synthetic cluster.

    The synthetic cluster is generated according to the parameters given in
    the ``params`` dictionary and the already calibrated
    :py:class:`Synthetic` object.

    :param params: Dictionary containing the values for the fundamental parameters.
    :type params: dict
    :param met_age_dict: Metallicity and age grid metadata.
    :type met_age_dict: dict
    :param met_tol: Tolerance for metallicity values matching
    :type met_tol: float
    :param age_tol: Tolerance for log(age) values matching
    :type age_tol: float
    :param def_params: Dictionary of default model parameters. This dictionary must
        include values for all the parameters.
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
    :param st_dist_mass: Sampled IMF mass distributions.
    :type st_dist_mass: list
    :param mass_idx_lookup: Pre-computed search indexes of the sampled IMF masses
        into the mass grid of the theoretical isochrones, generated by
        :py:func:`synth_cluster_priv.mass_idx_lookup`.
    :type mass_idx_lookup: dict
    :param max_mag_syn: Maximum synthetic magnitude allowed.
    :type max_mag_syn: float
    :param err_dist_synth: Error distributions used to perturb synthetic photometry.
    :type err_dist_synth: list[np.ndarray]
    :param N_synth_stars: Number of synthetic stars to generate.
    :type N_synth_stars: int
    :param show_errors: If ``True``, print the reason when an empty synthetic cluster
        is returned
    :type show_errors: bool
    :param isoch_flag: ``True`` if the isochrone should be returned along with the
        synthetic cluster, ``False`` otherwise.
    :type isoch_flag: bool

    :return: Synthetic isochrone plus synthetic array; if ``isoch_flag=False``,
        an empty array is returned instead of the isochrone.
        The synthetic array contains the data ``[mag, c1, (c2), mass, mass_b]``,
        where ``mag`` is the magnitude, ``c1`` is the color, ``c2`` is the
        optional second color, and ``mass, mass_b`` are the masses of the single
        and secondary components of the binary systems, respectively
        (if generated). If the system is a single star, then ``mass_b==np.nan``.
    :rtype: tuple[np.ndarray, np.ndarray]
    """

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    met, loga, alpha, beta, av, dr, rv, dm, ml, mh, al, ah = properModel(
        met_age_dict, def_params, params, met_tol, age_tol
    )

    # If (z, a) are both fixed, use the single processed isochrone
    if ml == al == mh == ah == 0:
        # The np.array() is important to avoid overwriting 'theor_tracks'
        isochrone = np.array(theor_tracks[0][0])
        # Grid node that the masses are taken from
        m0, a0 = 0, 0
    else:
        # Generate a weighted average isochrone from the (z, log(age)) values in
        # the 'model'. (m0, a0) is the closest grid node, from which the masses
        # are taken
        isochrone, m0, a0 = zaWAverage(
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

    isoch_extin_dr_0 = np.array([])
    if isoch_flag:
        # Use DR=0 to avoid disrupting the isochrone. The 'dr' variable is not
        # modified because it is used below for the synthetic cluster
        isoch_extin_dr_0 = extinction(
            ext_law,
            ext_coefs,
            rand_floats["dr_frac"],
            m_ini_idx,
            binar_flag,
            av,
            0.0,
            rv,
            isoch_moved.copy(),  # IMPORTANT, avoids double perturbing
        )

    # Apply extinction correction.
    # NOTE: the differential reddening is drawn for each isochrone point, and the
    # synthetic stars are interpolated between two points in mass_interp(). A star's
    # absorption is thus a weighted average of two independent draws, which makes the
    # realized spread of the absorption ~sqrt(2/3) (~18% narrower) that of the
    # selected DR model, and correlates the absorption of stars with similar masses.
    # Applying the DR per star (after the interpolation and binarity) gives the
    # nominal distributions, but was found to be ~1.4-1.6X slower per model
    isoch_extin = extinction(
        ext_law,
        ext_coefs,
        rand_floats["dr_frac"],
        m_ini_idx,
        binar_flag,
        av,
        dr,
        rv,
        isoch_moved,
    )

    # Identify the isochrone points below the maximum magnitude. The isochrone
    # itself is not cut (which would copy it): the points are selected through
    # their indexes in mass_interp()
    idx_cut, msk_cut = cut_max_mag(isoch_extin, max_mag_syn)
    if idx_cut.size == 0:
        if show_errors:
            print(
                f"No stars left after max magnitude cut (met={met:.5f},"
                + f"log(age)={loga:.3f}, Av={av:.3f}, DR={dr:.3f}, dm={dm:.3f})"
            )
        return isoch_extin_dr_0, np.array([])

    # Pre-computed search indexes for this IMF sample and grid node
    mass_idx = mass_idx_lookup[(ml, al, m0, a0)]

    # Interpolate IMF's sampled masses into the isochrone.
    isoch_mass = mass_interp(
        isoch_extin,
        idx_cut,
        msk_cut,
        m_ini_idx,
        st_dist_mass[ml][al],
        N_synth_stars,
        mass_idx,
        # binar_flag,
    )
    if isoch_mass.size == 0:
        if show_errors:
            print(
                f"No stars left after IMF interpolation (met={met:.5f},"
                + f"log(age)={loga:.3f}, Av={av:.3f}, DR={dr:.3f}, dm={dm:.3f})"
            )
        return isoch_extin_dr_0, np.array([])

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

    return isoch_extin_dr_0, synth_clust


def properModel(
    met_age_dict: dict,
    def_params: dict,
    fit_params: dict,
    met_tol: float,
    age_tol: float,
) -> tuple[float, float, float, float, float, float, float, float, int, int, int, int]:
    """Define the 'proper' model with values for (z, a) taken from its grid,
    and filled values for those parameters that are fixed.

    :param met_age_dict: Dictionary of metallicity and age values.
    :type met_age_dict: dict
    :param def_params: Dictionary of default parameters.
    :type def_params: dict
    :param fit_params: Dictionary of fitted parameters.
    :type fit_params: dict
    :param met_tol: Tolerance for metallicity values matching
    :type met_tol: float
    :param age_tol: Tolerance for log(age) values matching
    :type age_tol: float

    :raises ValueError: If (z, a) are outside of the isochrones' ranges, if any of
     ``alpha, beta, Av, DR, dm`` is negative, or if ``Rv <= 0``.

    :returns: A tuple containing all fundamental parameters, and the indexes of the
     (z, a) values in the grid that define the box that enclose the proper (z, a)
     values.
    :rtype: tuple[float, float, float, float, float, float, float, float, int, int, int, int]

    """
    # Combine `fit_params` with default parameter values in `def_params`, updating
    # duplicated values with those in `fit_params`
    fit_params = def_params | fit_params

    def _check_param(name, tol):
        if len(met_age_dict[name]) == 1:
            if abs(fit_params[name] - met_age_dict[name][0]) > tol:
                raise ValueError(
                    f"Parameter {name}={fit_params[name]:.5f} != fixed value"
                    + f" {met_age_dict[name][0]}"
                )
        elif not (met_age_dict[name][0] <= fit_params[name] <= met_age_dict[name][-1]):
            raise ValueError(
                f"Parameter {name}={fit_params[name]:.5f} outside isochrone"
                + f" range [{met_age_dict[name][0]}, {met_age_dict[name][-1]}]"
            )

    _check_param("met", met_tol)
    _check_param("loga", age_tol)

    # Parameters that can not be negative (Rv can not be zero either)
    for name in ("alpha", "beta", "Av", "DR", "dm"):
        if fit_params[name] < 0.0:
            raise ValueError(f"Parameter {name}={fit_params[name]:.5f} must be >= 0")
    if fit_params["Rv"] <= 0.0:
        raise ValueError(f"Parameter Rv={fit_params['Rv']:.5f} must be > 0")

    ml = mh = 0
    if len(met_age_dict["met"]) > 1:
        mh = min(
            len(met_age_dict["met"]) - 1,
            np.searchsorted(met_age_dict["met"], fit_params["met"]),
        )
        # prevent negative index
        ml = max(0, mh - 1)

    al = ah = 0
    if len(met_age_dict["loga"]) > 1:
        ah = min(
            len(met_age_dict["loga"]) - 1,
            np.searchsorted(met_age_dict["loga"], fit_params["loga"]),
        )
        al = max(0, ah - 1)

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
) -> tuple[np.ndarray, int, int]:
    """Generate a new weighted isochrone using bilinear interpolation over the
    orthogonal (z, loga) grid space.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [mag, c1, (c2), M_ini, M_b, magb, c1b, (c2b)]
    where:
    mag:  magnitude
    cX:  colors
    M_ini: initial mass
    M_b:  binary masses
    magb: binary magnitude
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

    :returns: Weighted isochrone, and the (metallicity, age) indexes of the closest
     grid node (the masses of the weighted isochrone are taken from this node).
    :rtype: tuple[np.ndarray, int, int]
    """

    # The four points in the (z, age) grid that define the box that contains
    # the model value (z_model, a_model)
    z1 = met_age_dict["met"][ml]
    z2 = met_age_dict["met"][mh]
    a1 = met_age_dict["loga"][al]
    a2 = met_age_dict["loga"][ah]

    # Compute normalized coordinates (t, u) within the local cell bounds
    # Safe division handling for edge cases where grid width might be zero
    dz = z2 - z1
    da = a2 - a1

    t = (z_model - z1) / dz if dz > 1e-15 else 0.0
    u = (a_model - a1) / da if da > 1e-15 else 0.0

    # Compute bilinear basis weights
    # Let om_t = (1 - t) and om_u = (1 - u)
    om_t = 1.0 - t
    om_u = 1.0 - u

    # Weights
    w11 = om_t * om_u  # Weight for (ml, al)
    w12 = om_t * u  # Weight for (ml, ah)
    w21 = t * om_u  # Weight for (mh, al)
    w22 = t * u  # Weight for (mh, ah)

    # Identify the closest grid point using a geometric boundary shortcut.
    # Since t and u map strictly to [0, 1], crossing 0.5 determines the closest node.
    m0 = mh if t >= 0.5 else ml
    a0 = ah if u >= 0.5 else al

    if _HAS_NUMBA:
        # Compiled version: identical results, see _zaWAverage_nb()
        isochrone = _zaWAverage_nb(
            theor_tracks, ml, mh, al, ah, w11, w12, w21, w22, m0, a0, m_ini_idx
        )
        return isochrone, m0, a0

    # Core interpolation via in-place accumulation
    isochrone = theor_tracks[ml][al] * w11
    isochrone += theor_tracks[ml][ah] * w12
    isochrone += theor_tracks[mh][al] * w21
    isochrone += theor_tracks[mh][ah] * w22

    # DO NOT average the masses or their distribution will be lost.
    # Overwrite mutated mass values with the exact discrete distribution from the closest node
    closest_iso = theor_tracks[m0][a0]
    isochrone[m_ini_idx] = closest_iso[m_ini_idx]
    # Now for the secondary masses
    isochrone[m_ini_idx + 1] = closest_iso[m_ini_idx + 1]

    # temp_plot(
    #     theor_tracks,
    #     isochrone,
    #     m_ini_idx,
    #     ml,
    #     mh,
    #     al,
    #     ah,
    #     z_model,
    #     a_model
    # )

    return isochrone, m0, a0


if _HAS_NUMBA:

    @njit(cache=True)
    def _zaWAverage_nb(
        theor_tracks: np.ndarray,
        ml: int,
        mh: int,
        al: int,
        ah: int,
        w11: float,
        w12: float,
        w21: float,
        w22: float,
        m0: int,
        a0: int,
        m_ini_idx: int,
    ) -> np.ndarray:
        """Compiled (numba) version of the weighted average in
        :py:func:`zaWAverage`. A single pass reads the four grid isochrones and
        writes the weighted sum, with the same operations in the same order as
        the NumPy version (identical results), without temporary arrays. The
        two mass rows are copied from the closest node (m0, a0) instead of being
        averaged and then overwritten.
        """
        t11 = theor_tracks[ml, al]
        t12 = theor_tracks[ml, ah]
        t21 = theor_tracks[mh, al]
        t22 = theor_tracks[mh, ah]
        t00 = theor_tracks[m0, a0]
        N_rows, N_pts = t11.shape
        isochrone = np.empty((N_rows, N_pts), dtype=t11.dtype)
        for r in range(N_rows):
            if r == m_ini_idx or r == m_ini_idx + 1:
                # DO NOT average the masses (see zaWAverage())
                for k in range(N_pts):
                    isochrone[r, k] = t00[r, k]
            else:
                for k in range(N_pts):
                    v = t11[r, k] * w11
                    v += t12[r, k] * w12
                    v += t21[r, k] * w21
                    v += t22[r, k] * w22
                    isochrone[r, k] = v
        return isochrone


def temp_plot(
    theor_tracks: np.ndarray,
    isochrone: np.ndarray,
    m_ini_idx: int,
    ml: int,
    mh: int,
    al: int,
    ah: int,
    z_model: float,
    a_model: float,
) -> None:
    """
    Temporary plotting function to visualize the bilinear interpolation of the
    theoretical tracks and the resulting isochrone. This function is intended for
    debugging and verification purposes.

    :param theor_tracks: Theoretical isochrones
    :type theor_tracks: np.ndarray
    :param isochrone: Interpolated isochrone
    :type isochrone: np.ndarray
    :param m_ini_idx: Index of the initial mass
    :type m_ini_idx: int
    :param ml: Index of the lower metallicity grid value
    :type ml: int
    :param mh: Index of the upper metallicity grid value
    :type mh: int
    :param al: Index of the lower age grid value
    :type al: int
    :param ah: Index of the upper age grid value
    :type ah: int
    :param z_model: Metallicity of the interpolated isochrone
    :type z_model: float
    :param a_model: Age of the interpolated isochrone
    :type a_model: float
    """
    import matplotlib.pyplot as plt

    print("ml, mh, al, ah:", ml, mh, al, ah)
    print("z_model, a_model:", z_model, a_model)

    plt.subplot(221)
    N_pts = len(theor_tracks[ml][al][m_ini_idx])
    plt.scatter(
        # theor_tracks[ml][al][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[ml][al][0],
        c="b",
        alpha=0.5,
        label="mlal",
    )
    plt.scatter(
        # theor_tracks[ml][ah][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[ml][ah][0],
        c="r",
        alpha=0.5,
        label="mlah",
    )
    plt.scatter(
        # theor_tracks[mh][al][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[mh][al][0],
        c="cyan",
        alpha=0.5,
        label="mhal",
    )
    plt.scatter(
        # theor_tracks[mh][ah][m_ini_idx],
        np.arange(0, N_pts),
        theor_tracks[mh][ah][0],
        c="orange",
        alpha=0.5,
        label="mhah",
    )
    plt.xlabel("mass")
    plt.ylabel("mag")
    # plt.scatter(isochrone[m_ini_idx], isochrone[0], marker="x", c="k")
    plt.scatter(np.arange(0, len(isochrone[0])), isochrone[0], marker="x", c="k")
    plt.legend()

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
    dr_frac: np.ndarray,
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

    With differential reddening (dr>0) the total absorption of each point is:

    A_i = Av + f_i * dr

    where 'f_i' (in [0, 1]) are the fractions in 'dr_frac', drawn from the
    distribution selected with 'DR_distribution' (see
    :py:func:`synth_cluster_priv.dr_fractions`). Hence Av is the minimum and
    Av + dr the maximum absorption.

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
    :param dr_frac: Fractions of the differential reddening applied to each point.
    :type dr_frac: np.ndarray
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

    if _HAS_NUMBA and ext_law == "GAIADR3":
        # Compiled version: identical results, see _extinction_gaiadr3_nb()
        # dr_kind: 0 -> no differential reddening (scalar A), 1 -> A = Av + f*dr
        dr_kind = 1 if dr > 0.0 else 0
        _extinction_gaiadr3_nb(
            isochrone, dr_frac, dr_kind, Av, Av, dr, m_ini_idx, binar_flag
        )
        return isochrone

    if _HAS_NUMBA and ext_law == "CCMO":
        # Compiled version: identical results, see _extinction_ccmo_nb()
        ec_mag, ec_col1, ec_col2 = ccmo_coefs(ext_coefs, Rv)
        _extinction_ccmo_nb(
            isochrone,
            dr_frac,
            1 if dr > 0.0 else 0,
            Av,
            dr,
            m_ini_idx,
            binar_flag,
            ec_mag,
            ec_col1,
            ec_col2,
            len(ext_coefs) > 2,
        )
        return isochrone

    Av_dr = Av
    if dr > 0.0:
        Ns = isochrone.shape[-1]
        # A_i = Av + f_i * dr, in [Av, Av + dr]
        Av_dr = dr_frac[:Ns] * dr
        Av_dr += Av

    if ext_law == "CCMO":
        # Magnitude
        ec_mag = ext_coefs[0][0] + ext_coefs[0][1] / Rv
        # First color
        ec_col1 = (ext_coefs[1][0][0] + ext_coefs[1][0][1] / Rv) - (
            ext_coefs[1][1][0] + ext_coefs[1][1][1] / Rv
        )
        Ax = ec_mag * Av_dr
        Ex1 = ec_col1 * Av_dr
    elif ext_law == "GAIADR3":
        # If this model is used the first color is always expected to be BP-RP
        # BP_RP = isochrone[1]
        Ax, Ex1 = gaiadr3_ext(isochrone[1], Av_dr)
    else:
        raise ValueError(f"Unknown extinction law: {ext_law}")

    isochrone[0] += Ax  # Magnitude
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


if _HAS_NUMBA:

    @njit(cache=True)
    def _extinction_gaiadr3_nb(
        isochrone: np.ndarray,
        dr_frac: np.ndarray,
        dr_kind: int,
        A_scalar: float,
        Av: float,
        dr: float,
        m_ini_idx: int,
        binar_flag: bool,
    ) -> None:
        """Compiled (numba) version of the GAIADR3 branch of
        :py:func:`extinction`, modifying 'isochrone' in place.

        A single loop over the isochrone points computes the (differentially
        reddened) absorption, A_G and E(BP-RP) (same nested form as
        :py:func:`_ext_horner_XA` and :py:func:`_ext_horner_X`), and adds them to
        the magnitude and color rows. The floating point operations are the same
        and in the same order as in the NumPy version, so the results are
        identical.

        dr_kind: 0 -> scalar absorption 'A_scalar' (no differential reddening),
        1 -> differential reddening A = Av + f * dr, with f from 'dr_frac'.
        """
        g = _G_COEFFS
        b = _BPRP_COEFFS

        # Scalar absorption: coefficients of the cubic polynomial in X
        A = A_scalar
        gk0 = A * (g[0] + A * (g[4] + A * (g[5] + A * g[6])))
        gk1 = A * (g[1] + A * (g[7] + A * g[9]))
        gk2 = A * (g[2] + A * g[8])
        gk3 = A * g[3]
        bk0 = A * (b[0] + A * (b[4] + A * (b[5] + A * b[6])))
        bk1 = A * (b[1] + A * (b[7] + A * b[9]))
        bk2 = A * (b[2] + A * b[8])
        bk3 = A * b[3]

        i_mb = m_ini_idx + 2
        i_cb = m_ini_idx + 3
        for i in range(isochrone.shape[1]):
            # Clamped to the range of validity of the coefficients, as in
            # gaiadr3_ext()
            X = isochrone[1, i]
            if X < _X_MIN:
                X = _X_MIN
            elif X > _X_MAX:
                X = _X_MAX

            if dr_kind == 0:
                Ax = gk3 * X
                Ax += gk2
                Ax *= X
                Ax += gk1
                Ax *= X
                Ax += gk0
                Ex = bk3 * X
                Ex += bk2
                Ex *= X
                Ex += bk1
                Ex *= X
                Ex += bk0
            else:
                A = dr_frac[i] * dr
                A += Av

                # G: A * [P0(X) + A * (P1(X) + A * (P2(X) + c6 * A))]
                R = g[9] * X
                R += g[5]
                R += g[6] * A
                P = g[8] * X
                P += g[7]
                P *= X
                P += g[4]
                R *= A
                R += P
                P = X * g[3]
                P += g[2]
                P *= X
                P += g[1]
                P *= X
                P += g[0]
                R *= A
                R += P
                R *= A
                Ax = R
                # BP-RP
                R = b[9] * X
                R += b[5]
                R += b[6] * A
                P = b[8] * X
                P += b[7]
                P *= X
                P += b[4]
                R *= A
                R += P
                P = X * b[3]
                P += b[2]
                P *= X
                P += b[1]
                P *= X
                P += b[0]
                R *= A
                R += P
                R *= A
                Ex = R

            isochrone[0, i] += Ax  # Magnitude
            isochrone[1, i] += Ex  # First color
            if binar_flag:
                isochrone[i_mb, i] += Ax  # Magnitude (binary system)
                isochrone[i_cb, i] += Ex  # First color (binary system)

    @njit(cache=True)
    def _extinction_ccmo_nb(
        isochrone: np.ndarray,
        dr_frac: np.ndarray,
        dr_kind: int,
        Av: float,
        dr: float,
        m_ini_idx: int,
        binar_flag: bool,
        ec_mag: float,
        ec_col1: float,
        ec_col2: float,
        col2_flag: bool,
    ) -> None:
        """Compiled (numba) version of the CCMO branch of :py:func:`extinction`,
        modifying 'isochrone' in place. The coefficients ``ec_*`` (``a + b / Rv``
        values, see :py:func:`ccmo_coefs`) multiply the (differentially reddened)
        absorption of each point. The floating point operations are the same and in
        the same order as in the NumPy version, so the results are identical.

        dr_kind: 0 -> scalar absorption Av (no differential reddening),
        1 -> differential reddening A = Av + f * dr, with f from 'dr_frac'.
        """
        i_mb = m_ini_idx + 2
        i_cb = m_ini_idx + 3
        i_c2b = m_ini_idx + 4
        A = Av
        for i in range(isochrone.shape[1]):
            if dr_kind == 1:
                A = dr_frac[i] * dr
                A += Av
            Ax = ec_mag * A
            Ex1 = ec_col1 * A
            isochrone[0, i] += Ax  # Magnitude
            isochrone[1, i] += Ex1  # First color
            if binar_flag:
                isochrone[i_mb, i] += Ax  # Magnitude (binary system)
                isochrone[i_cb, i] += Ex1  # First color (binary system)
            if col2_flag:
                Ex2 = ec_col2 * A
                isochrone[2, i] += Ex2  # Second color
                if binar_flag:
                    isochrone[i_c2b, i] += Ex2  # Second color (binary system)


def ccmo_coefs(ext_coefs: list | np.ndarray, Rv: float) -> tuple[float, float, float]:
    """Coefficients of the CCMO extinction law for a given ``Rv``: the ratios
    ``A_mag / Av`` and ``E(color) / Av`` for the magnitude and the color(s), see
    :py:func:`extinction`.

    :param ext_coefs: List of extinction coefficients, ``[ec_mag, ec_col1, ...]``
    :type ext_coefs: list | np.ndarray
    :param Rv: Total-to-selective extinction ratio.
    :type Rv: float

    :returns: Coefficients for the magnitude, the first color and the second color
        (``0`` if there is no second color).
    :rtype: tuple[float, float, float]
    """
    ec_mag = ext_coefs[0][0] + ext_coefs[0][1] / Rv
    ec_col1 = (ext_coefs[1][0][0] + ext_coefs[1][0][1] / Rv) - (
        ext_coefs[1][1][0] + ext_coefs[1][1][1] / Rv
    )
    ec_col2 = 0.0
    if len(ext_coefs) > 2:
        ec_col2 = (ext_coefs[2][0][0] + ext_coefs[2][0][1] / Rv) - (
            ext_coefs[2][1][0] + ext_coefs[2][1][1] / Rv
        )
    return float(ec_mag), float(ec_col1), float(ec_col2)


def gaiadr3_ext(X: np.ndarray, A: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Gaia DR3 extinction law for main sequence stars, returning the extinction
    in G and the color excess in BP-RP directly:

    A_G = ec_G(X, A) * A ; E(BP-RP) = ec_BPRP(X, A) * A

    where ec_G, ec_BPRP are the coefficients of the law: polynomials in X and A with
    the coefficients in ``_G_COEFFS`` and ``_BPRP_COEFFS``,

    ec = c0 + c1*X + c2*X^2 + c3*X^3 + c4*A + c5*A^2 + c6*A^3 + c7*A*X + c8*A*X^2
         + c9*X*A^2

    The result is equal (to floating point rounding, ~1e-15 mag) to evaluating these
    polynomials and multiplying them by ``A``, but it is evaluated faster by
    writing the polynomial in nested (Horner) form, both in X and A, with
    in-place operations to avoid allocating temporary arrays. If ``A`` is a scalar
    (no differential reddening) the powers of ``A`` are absorbed into the
    coefficients and only a cubic polynomial in X is evaluated.

    The coefficients are evaluated with the colors clamped to the range for which
    they were fitted, ``[_X_MIN, _X_MAX]``: stars outside of it are assigned the
    coefficients of the closest valid color.

    :param X: Array of (intrinsic) BP-RP colors.
    :type X: np.ndarray
    :param A: Total absorption (eventually containing differential reddening).
    :type A: float | np.ndarray

    :returns: Extinction in G and color excess in BP-RP.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    X = np.clip(X, _X_MIN, _X_MAX)
    if isinstance(A, np.ndarray):
        return _ext_horner_XA(_G_COEFFS, X, A), _ext_horner_XA(_BPRP_COEFFS, X, A)
    return _ext_horner_X(_G_COEFFS, X, A), _ext_horner_X(_BPRP_COEFFS, X, A)


def _ext_horner_XA(c: tuple, X: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Evaluate ``A * ec(X, A)`` for array ``A``, written as:

    A * [P0(X) + A * (P1(X) + A * (P2(X) + c6 * A))]

    with P0 = c0 + c1*X + c2*X^2 + c3*X^3 ; P1 = c4 + c7*X + c8*X^2 ;
    P2 = c5 + c9*X

    :param c: The 10 coefficients of the extinction law polynomial
    :type c: tuple
    :param X: Array of (intrinsic) BP-RP colors
    :type X: np.ndarray
    :param A: Array of absorptions, one per color
    :type A: np.ndarray

    :return: Extinction (or color excess) for each color
    :rtype: np.ndarray
    """
    # P2(X) + c6 * A
    R = c[9] * X
    R += c[5]
    R += c[6] * A
    # P1(X)
    P = c[8] * X
    P += c[7]
    P *= X
    P += c[4]
    R *= A
    R += P
    # P0(X), reusing the 'P' buffer
    np.multiply(X, c[3], out=P)
    P += c[2]
    P *= X
    P += c[1]
    P *= X
    P += c[0]
    R *= A
    R += P
    R *= A
    return R


def _ext_horner_X(c: tuple, X: np.ndarray, A: float) -> np.ndarray:
    """Evaluate ``A * ec(X, A)`` for scalar ``A`` as a cubic polynomial in X:

    k0 + k1*X + k2*X^2 + k3*X^3

    where the 'k' coefficients contain all the powers of A.

    :param c: The 10 coefficients of the extinction law polynomial
    :type c: tuple
    :param X: Array of (intrinsic) BP-RP colors
    :type X: np.ndarray
    :param A: Absorption, the same for all colors
    :type A: float

    :return: Extinction (or color excess) for each color
    :rtype: np.ndarray
    """
    k0 = A * (c[0] + A * (c[4] + A * (c[5] + A * c[6])))
    k1 = A * (c[1] + A * (c[7] + A * c[9]))
    k2 = A * (c[2] + A * c[8])
    k3 = A * c[3]
    R = k3 * X
    R += k2
    R *= X
    R += k1
    R *= X
    R += k0
    return R


def cut_max_mag(
    isoch_moved: np.ndarray, max_mag_syn: float
) -> tuple[np.ndarray, np.ndarray]:
    """Identify the isochrone points with magnitude values smaller than the maximum
    observed value.

    The isochrone is not cut here: ``isoch_moved[:, idx]`` would copy the entire
    array, which is not needed since :py:func:`mass_interp` only reads the
    points it interpolates. Instead, the indexes of the kept points are returned.

    :param isoch_moved: Isochrone array.
    :type isoch_moved: np.ndarray
    :param max_mag_syn: Maximum magnitude value.
    :type max_mag_syn: float

    :returns: Indexes of the points kept (sorted), and the boolean mask of the
     points kept.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    # Discard stars in isochrone beyond max_mag_syn limit.
    msk = isoch_moved[0] < max_mag_syn
    return np.flatnonzero(msk), msk


def mass_interp(
    isochrone: np.ndarray,
    idx_cut: np.ndarray,
    msk_cut: np.ndarray,
    m_ini_idx: int,
    st_dist_mass: np.ndarray,
    N_synth_stars: int,
    mass_idx: np.ndarray,
    # binar_flag: bool,
) -> np.ndarray:
    """For each mass in the sampled IMF mass distribution, interpolate its value
    (and those of all the sub-arrays in the isochrone) into the isochrone.

    Only the isochrone points given by 'idx_cut' (those below the maximum
    magnitude, see :py:func:`cut_max_mag`) are used. The result is identical to
    interpolating into the cut isochrone ``isochrone[:, idx_cut]``, without
    having to copy it.

    Masses that fall outside of the (cut) isochrone's mass range are rejected.

    :param isochrone: Isochrone array (not cut).
    :type isochrone: np.ndarray
    :param idx_cut: Sorted indexes of the isochrone points kept by
        :py:func:`cut_max_mag`.
    :type idx_cut: np.ndarray
    :param msk_cut: Boolean mask of the isochrone points kept by
        :py:func:`cut_max_mag`.
    :type msk_cut: np.ndarray
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int
    :param st_dist_mass: Array of sampled masses.
    :type st_dist_mass: np.ndarray
    :param N_synth_stars: Number of observed stars.
    :type N_synth_stars: int
    :param mass_idx: Pre-computed search indexes of 'st_dist_mass' into the full
        (before the maximum magnitude cut) mass grid of the isochrone, see
        :py:func:`synth_cluster_priv.mass_idx_lookup`.
    :type mass_idx: np.ndarray

    :returns: Interpolated isochrone array.
    :rtype: np.ndarray
    """
    if _HAS_NUMBA:
        # Compiled version: identical results, see _mass_interp_nb()
        isoch_mass = _mass_interp_nb(
            isochrone,
            idx_cut,
            msk_cut,
            m_ini_idx,
            st_dist_mass,
            N_synth_stars,
            mass_idx,
        )
        if isoch_mass.shape[1] == 0:
            return np.array([])
        return isoch_mass

    # ==IMPORTANT==
    # Assumes `mass_ini=isochrone[m_ini_idx]` is sorted min to max. This is checked
    # when loading the isochrones in isochrones_priv.interp_isochrones()
    # ==IMPORTANT==
    mass_ini = isochrone[m_ini_idx]

    # Filter masses in the IMF sampling that are outside of the mass
    # range given by the cut isochrone (st_dist_mass[0]: sampled masses from IMF)

    # # msk_min = (st_dist_mass[0] >= mass_ini.min())
    # # msk_max = (st_dist_mass[0] <= mass_ini.max())
    # # (~msk_min).sum(): stars lost below the minimum mass (photometric)
    # # (~msk_max).sum(): stars lost above the maximum mass (evolutionary)
    # msk_m = (st_dist_mass >= mass_ini.min()) & (st_dist_mass <= mass_ini.max())

    # Extract up to 'N_synth_stars' masses sampled from an IMF, within the mass range
    # idx = np.nonzero(msk_m)[0][:N_synth_stars]
    idx = first_in_range(
        st_dist_mass, mass_ini[idx_cut[0]], mass_ini[idx_cut[-1]], N_synth_stars
    )

    mass_dist = st_dist_mass[idx]

    if mass_dist.size == 0:
        return np.array([])

    mass_dist = np.ascontiguousarray(mass_dist, dtype=mass_ini.dtype)

    # Indexes into the *cut* mass grid, i.e.:
    # np.searchsorted(mass_ini[idx_cut], mass_dist)
    # Obtained by converting the pre-computed indexes into the full mass grid, into
    # indexes into the cut grid: the number of grid points with a mass smaller than
    # the sampled mass that survived the cut
    N_kept = np.empty(msk_cut.size + 1, dtype=np.intp)
    N_kept[0] = 0
    np.cumsum(msk_cut, out=N_kept[1:])
    hi_cut = N_kept[mass_idx[idx]]

    # Map back to indexes into the full isochrone. Note that hi_cut=0 gives
    # lo=idx_cut[-1], exactly as indexing the cut isochrone with -1 would
    x_new_indices = idx_cut[hi_cut]
    lo = idx_cut[hi_cut - 1]

    # Interpolate the sampled stars (masses) into the isochrone
    isoch_mass = interp_mass_isoch(isochrone, mass_ini, mass_dist, x_new_indices, lo)

    return isoch_mass


if _HAS_NUMBA:

    @njit(cache=True)
    def _mass_interp_nb(
        isochrone: np.ndarray,
        idx_cut: np.ndarray,
        msk_cut: np.ndarray,
        m_ini_idx: int,
        st_dist_mass: np.ndarray,
        N_synth_stars: int,
        mass_idx: np.ndarray,
    ) -> np.ndarray:
        """Compiled (numba) version of :py:func:`mass_interp` (with the
        pre-computed 'mass_idx' search indexes). It performs the same floating
        point operations in the same order, so the results are identical to the
        NumPy version, but it avoids most temporary arrays and Python calls:

        1. Select the first 'N_synth_stars' sampled masses within the mass range
           of the cut isochrone, stopping as soon as they are found. The loop is
           written without branches, which is faster.
        2. Bracketing points: for each point 'k' of the full mass grid, the next
           kept point at or after 'k' and the last kept point before 'k' are
           tabulated once. With 'mass_idx' this replaces the cumulative count
           used in the NumPy version.
        3. Interpolate the isochrone linearly, one row at a time (rows are
           contiguous in memory).
        """
        mass_ini = isochrone[m_ini_idx]
        N_cut = idx_cut.size
        v_min = mass_ini[idx_cut[0]]
        v_max = mass_ini[idx_cut[N_cut - 1]]

        # 1. First 'N_synth_stars' masses within [v_min, v_max]. The index is
        # always written and the counter only advances for masses in range.
        sel = np.empty(st_dist_mass.size, dtype=np.intp)
        N = 0
        for i in range(st_dist_mass.size):
            m = st_dist_mass[i]
            sel[N] = i
            N += (m >= v_min) & (m <= v_max)
            if N >= N_synth_stars:
                break
        N = min(N, N_synth_stars)

        # 2. Next kept point at or after k (nxt), and last kept point before k
        # (prv). If there is no kept point before k, prv=idx_cut[-1], as indexing
        # the cut isochrone with -1 does in the NumPy version.
        N_grid = msk_cut.size
        nxt = np.empty(N_grid + 1, dtype=np.intp)
        prv = np.empty(N_grid + 1, dtype=np.intp)
        last = idx_cut[N_cut - 1]
        for k in range(N_grid + 1):
            prv[k] = last
            if k < N_grid and msk_cut[k]:
                last = k
        nxt_k = N_grid
        for k in range(N_grid, -1, -1):
            if k < N_grid and msk_cut[k]:
                nxt_k = k
            nxt[k] = nxt_k

        hi = np.empty(N, dtype=np.intp)
        lo = np.empty(N, dtype=np.intp)
        w = np.empty(N, dtype=isochrone.dtype)
        for j in range(N):
            i = sel[j]
            k = mass_idx[i]
            hi[j] = nxt[k]
            lo[j] = prv[k]
            x_lo = mass_ini[lo[j]]
            w_j = st_dist_mass[i] - x_lo
            w_j /= mass_ini[hi[j]] - x_lo
            w[j] = w_j

        # 3. Linear interpolation, y_lo + (y_hi - y_lo) * w
        N_rows = isochrone.shape[0]
        isoch_mass = np.empty((N_rows, N), dtype=isochrone.dtype)
        for r in range(N_rows):
            row = isochrone[r]
            out = isoch_mass[r]
            for j in range(N):
                y_lo = row[lo[j]]
                out[j] = (row[hi[j]] - y_lo) * w[j] + y_lo

        return isoch_mass


def first_in_range(
    arr: np.ndarray, v_min: float, v_max: float, N_max: int
) -> np.ndarray:
    """Return the indexes of the first (at most) 'N_max' elements of 'arr' that
    fall within the closed range [v_min, v_max].

    Equivalent to ``np.nonzero((arr >= v_min) & (arr <= v_max))[0][:N_max]`` but
    the mask is evaluated over progressively larger chunks of 'arr', stopping as
    soon as 'N_max' elements are found. Since the IMF sampled array is usually
    much larger than the number of required stars (e.g.: ~55000 vs ~1500), this
    avoids masking the entire array in most cases.

    :param arr: Array of values.
    :type arr: np.ndarray
    :param v_min: Minimum value of the range.
    :type v_min: float
    :param v_max: Maximum value of the range.
    :type v_max: float
    :param N_max: Maximum number of indexes to return.
    :type N_max: int

    :returns: Array of indexes.
    :rtype: np.ndarray
    """
    N_tot = arr.size
    # The initial chunk size of 4*N_max was found to be a good compromise for
    # typical in-range fractions; the chunk size is doubled in each iteration
    i_end = min(N_tot, 4 * N_max)
    i_start, N_found, idx_lst = 0, 0, []
    while True:
        chunk = arr[i_start:i_end]
        idx = np.flatnonzero((chunk >= v_min) & (chunk <= v_max))
        if i_start > 0:
            idx += i_start
        idx_lst.append(idx)
        N_found += idx.size
        if N_found >= N_max or i_end == N_tot:
            break
        i_start, i_end = i_end, min(N_tot, 2 * i_end)

    if len(idx_lst) == 1:
        return idx[:N_max]
    return np.concatenate(idx_lst)[:N_max]


# import _interp_mass_isoch as _c

# def interp_mass_isoch(
#     isoch_cut: np.ndarray,
#     mass_ini: np.ndarray,
#     mass_dist: np.ndarray,
# ) -> np.ndarray:
#     # C extension requires float64 C-contiguous arrays
#     return _c.interp_mass_isoch(
#         np.ascontiguousarray(isoch_cut, dtype=np.float64),
#         np.ascontiguousarray(mass_ini,  dtype=np.float64),
#         np.ascontiguousarray(mass_dist, dtype=np.float64),
#     )


def interp_mass_isoch(
    isoch_cut: np.ndarray,
    mass_ini: np.ndarray,
    mass_dist: np.ndarray,
    x_new_indices: np.ndarray | None = None,
    lo: np.ndarray | None = None,
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
    :param x_new_indices: Result of ``np.searchsorted(mass_ini, mass_dist)`` if
        already known; computed here otherwise.
    :type x_new_indices: np.ndarray | None
    :param lo: Indexes of the lower bracketing points, ``x_new_indices - 1`` if
        not given. Used when the isochrone points are selected through indexes
        (see :py:func:`mass_interp`).
    :type lo: np.ndarray | None

    :returns: Interpolated isochrone array.
    :rtype: np.ndarray
    """

    # *****************************************************************************
    # This is a stripped down version of scipy.inter1d (kind='linear') which
    # is now legacy code
    # (https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html)
    # Scipy now recommends using np.interp for simple linear interpolation but it
    # is quite slower.
    #
    # Note: If x_new[n] == x[m], then m is returned by searchsorted.
    # x_new_indices = np.searchsorted(mass_ini, mass_dist)
    # #
    # # Calculate the slope of regions that each x_new value falls in.
    # lo = x_new_indices - 1
    # x_lo = mass_ini[lo]
    # x_hi = mass_ini[x_new_indices]
    # # y_lo = isoch_cut[:, lo]
    # # y_hi = isoch_cut[:, x_new_indices]
    # y_lo = np.take(isoch_cut, lo, axis=1)
    # y_hi = np.take(isoch_cut, x_new_indices, axis=1)

    # # Compute slope and result in fewer steps
    # x_diff = mass_dist - x_lo
    # # Calculate the actual value for each entry in x_new.
    # # slope = (y_hi - y_lo) / (x_hi - x_lo)
    # isoch_mass = y_lo + (y_hi - y_lo) * (x_diff / (x_hi - x_lo))
    # # *****************************************************************************

    if x_new_indices is None:
        x_new_indices = np.searchsorted(mass_ini, mass_dist)
    #
    # Calculate the slope of regions that each x_new value falls in.
    if lo is None:
        lo = x_new_indices - 1
    x_lo = mass_ini[lo]
    # Interpolation weights: (x_new - x_lo) / (x_hi - x_lo)
    w = mass_dist - x_lo
    w /= mass_ini[x_new_indices] - x_lo

    # Calculate the actual value for each entry in x_new as:
    # y_lo + (y_hi - y_lo) * w
    # using in-place operations to avoid allocating temporary (D, N) arrays
    y_lo = np.take(isoch_cut, lo, axis=1)
    isoch_mass = np.take(isoch_cut, x_new_indices, axis=1)
    isoch_mass -= y_lo
    isoch_mass *= w
    isoch_mass += y_lo

    # *****************************************************************************
    # # This approach is 3X slower
    # isoch_T = np.ascontiguousarray(isoch_cut.T)   # shape (N_iso, D)

    # isoch_mass = np.empty((len(mass_dist), isoch_cut.shape[0]))  # (N_samp, D)
    # for i in range(isoch_cut.shape[0]):           # iterate over D rows
    #     isoch_mass[:, i] = np.interp(mass_dist, mass_ini, isoch_T[:, i])

    # isoch_mass = isoch_mass.T                     # back to (D, N_samp)
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

    if _HAS_NUMBA:
        # Compiled version: identical results, see _binarity_nb()
        _binarity_nb(alpha, beta, m_ini_idx, rand_unif_vals, isoch_mass)
        return isoch_mass[: m_ini_idx + 2]

    # Distribution of probability of binarity (multiplicity fraction) versus
    # primary masses.

    # Offner et al. (2022); Fig 1 (left), Table 1
    # b_p = np.clip(alpha + beta / (1.0 + 1.4 / isoch_mass[m_ini_idx]), 0.0, 1.0)
    # Same operations (hence identical values) as the line above, but in place.
    # Clipping to [0, 1] is not needed: the probabilities are only compared below
    # with uniform random values in [0, 1), so a value < 0 (> 1) gives the same
    # result as 0 (1)
    b_p = 1.4 / isoch_mass[m_ini_idx]
    b_p += 1.0
    np.divide(beta, b_p, out=b_p)
    b_p += alpha

    # Stars (masses) with the largest binary probabilities are selected
    # proportional to their probability
    bin_indxs = b_p > rand_unif_vals[:Ns]

    # Index of the binary magnitude: mag_binar = m_ini_idx + 1
    # Update array with new values of magnitudes and colors for the binary
    # systems. This does not change the primary mass values, just the magnitude and
    # color(s).
    for i in range(m_ini_idx):
        # isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 2 + i][bin_indxs]
        # Use np.copyto with where parameter (faster than boolean indexing)
        np.copyto(isoch_mass[i], isoch_mass[m_ini_idx + 2 + i], where=bin_indxs)

    # Update the binary systems' masses so that the secondary masses for
    # SINGLE systems are identified with a 'nan' value. Using np.copyto() on the
    # row is faster than the 2D boolean indexing isoch_mass[m_ini_idx + 1, ~bin_indxs]
    np.copyto(isoch_mass[m_ini_idx + 1], np.nan, where=~bin_indxs)

    # Return [mag, c1, (c2), mass, mass_b]
    return isoch_mass[: m_ini_idx + 2]


if _HAS_NUMBA:

    @njit(cache=True)
    def _scatter_add_nb(row: np.ndarray, order: np.ndarray, sigma: np.ndarray) -> None:
        """Compiled equivalent of ``row[order] += sigma[:order.size]``, used by
        :py:func:`add_errors`"""
        for j in range(order.size):
            row[order[j]] += sigma[j]


if _HAS_NUMBA:

    @njit(cache=True)
    def _binarity_nb(
        alpha: float,
        beta: float,
        m_ini_idx: int,
        rand_unif_vals: np.ndarray,
        isoch_mass: np.ndarray,
    ) -> None:
        """Compiled (numba) version of :py:func:`binarity` (with binaries),
        modifying 'isoch_mass' in place. A single loop over the stars computes
        the binary probability (same operations, in the same order, as the
        NumPy version: identical results), compares it with the uniform random
        value, and either replaces the magnitude and color(s) with those of the
        binary system, or flags the secondary mass of a single star with 'nan'.
        """
        mass = isoch_mass[m_ini_idx]
        mass_b = isoch_mass[m_ini_idx + 1]
        for j in range(isoch_mass.shape[1]):
            # alpha + beta / (1.0 + 1.4 / mass)
            b_p = 1.4 / mass[j]
            b_p += 1.0
            b_p = beta / b_p
            b_p += alpha
            if b_p > rand_unif_vals[j]:
                for i in range(m_ini_idx):
                    isoch_mass[i, j] = isoch_mass[m_ini_idx + 2 + i, j]
            else:
                mass_b[j] = np.nan


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
    # Indexes that sort the stars from the largest to the smallest magnitude, the
    # same order used by the error distributions (see error_distribution())
    mag_sort = np.argsort(-isoch_binar[0])

    if _HAS_NUMBA:
        # Compiled scatter-add, see _scatter_add_nb(). The sorting above is left
        # to NumPy: numba's argsort, and a compiled bucket sort, were not faster
        for i, sigma in enumerate(err_dist):
            _scatter_add_nb(isoch_binar[i], mag_sort, sigma)
        return isoch_binar

    for i, sigma in enumerate(err_dist):
        # Index the row first: fancy indexing a 1D view is noticeably faster
        # than the mixed (row, indexes) indexing of the 2D array
        row = isoch_binar[i]
        row[mag_sort] += sigma[:N]

    return isoch_binar


if __name__ == "__main__":
    main()
