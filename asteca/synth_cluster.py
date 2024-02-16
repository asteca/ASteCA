import numpy as np
from .modules.imfs import invTrnsfSmpl, sampleInv
from .modules import synth_clust_funcs as scf

#
# from .synth_clust import imfs
# from .synth_clust import zaWAverage
# from .synth_clust import move_isochrone
# from .synth_clust import cut_max_mag
# from .synth_clust import mass_interp
# from .synth_clust import extinction
# from .synth_clust import binarity
# from .synth_clust import add_errors


def add_binarity(self) -> np.ndarray:
    """
    For each theoretical isochrone defined.
        1. Draw random secondary masses for *all* stars.
        2. Interpolate magnitude values for these masses
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.

    theor_tracks.shape = (Nz, Na, N_cols, N_interp)

    If binarity is processed:
        N_cols: magnitude + color(s) + initial masses + unresolved mag +
            unresolved color(s) + secondary masses
    else:
        N_cols: magnitude + colors + initial mass

    """
    if self.alpha <= 0.:
        print("'alpha=0': binary systems are not generated")
        return self.isochs.theor_tracks
    print("Generating binary data")

    mag_idx = 0
    N_colors = len(self.isochs.color_filter_name)
    m_ini_idx = mag_idx + N_colors + 1

    # Extend to accommodate binary data
    Nz, Na, Nd, Ni = self.isochs.theor_tracks.shape
    theor_tracks = np.concatenate(
        (self.isochs.theor_tracks, np.zeros([Nz, Na, Nd, Ni])), axis=2
    )

    # For each metallicity defined.
    for mx, met in enumerate(theor_tracks):
        # For each age defined.
        for ax, isoch in enumerate(met):
            # Extract initial masses for this isochrone.
            mass_ini = isoch[m_ini_idx]

            # Mass-ratio distribution
            mass_ratios = scf.qDistribution(mass_ini, self.gamma)

            # Calculate random secondary masses of these binary stars
            # between bin_mass_ratio*m1 and m1, where m1 is the primary
            # mass.
            m2 = mass_ratios * mass_ini

            # Calculate unresolved binary magnitude
            mag = isoch[mag_idx]
            mag_m2 = np.interp(m2, mass_ini, mag)
            theor_tracks[mx][ax][m_ini_idx + 1] = scf.mag_combine(mag, mag_m2)

            # Calculate unresolved color for each color defined.
            for ic, color in enumerate(self.isochs.color_filter_name):
                mc1 = self.isochs.color_filters[mx][ax][color[0]]
                mc2 = self.isochs.color_filters[mx][ax][color[1]]
                f_m1 = np.interp(m2, mass_ini, mc1)
                f1 = scf.mag_combine(mc1, f_m1)
                f_m2 = np.interp(m2, mass_ini, mc2)
                f2 = scf.mag_combine(mc2, f_m2)
                theor_tracks[mx][ax][m_ini_idx + 1 + 1 + ic] = f1 - f2

            # Secondary masses
            theor_tracks[mx][ax][-1] = m2

    return theor_tracks


def extinction_coeffs(self) -> list:
    """
    Obtain extinction coefficients for all the observed filters and colors,
    in the order in which they are stored in theor_tracks.

    ext_coefs = [ec_mag, ec_col1, ...]
    """
    print("Obtaining extinction coefficients")

    # For the magnitude.
    # Effective wavelength in Armstrong.
    eff_wave = self.isochs.mag_color_lambdas[self.isochs.mag_filter_name]
    # CCM coefficient. Effective wavelength in inverse microns.
    ext_coefs = [scf.ccm_model(10000.0 / eff_wave)]

    # For colors.
    for color in self.isochs.color_filter_name:
        c_filt1, c_filt2 = color
        # Effective wavelength in Armstrong.
        eff_wave1 = self.isochs.mag_color_lambdas[c_filt1]
        eff_wave2 = self.isochs.mag_color_lambdas[c_filt2]
        ext_coefs.append(
            [scf.ccm_model(10000.0 / eff_wave1), scf.ccm_model(10000.0 / eff_wave2)]
        )

    return ext_coefs


def sample_imf(self, Nmets: int, Nages: int) -> list:
    """
    Returns the number of stars per interval of mass for the selected IMF.

    Parameters
    ----------
    IMF_name : str
      Name of the IMF to be used.
    max_mass: float
      Maximum mass defined.

    Returns
    -------
    st_dist_mass : list
      Tuple that contains: a given number of stars sampled from the selected
      IMF, that (approximately) sum to the associated total mass; the
      cumulative sum of those masses.
      One list per metallicity and age values defined is returned. This is to add some
      variety to the sampled IMF.
    """
    print(f"Sampling selected IMF ({self.IMF_name})")
    inv_cdf = invTrnsfSmpl(self.IMF_name)

    st_dist_mass = []
    for _ in range(Nmets):
        met_lst = []
        for _ in range(Nages):
            sampled_IMF = sampleInv(self.max_mass, inv_cdf)
            met_lst.append(sampled_IMF)
        st_dist_mass.append(met_lst)

    return st_dist_mass


def randVals(theor_tracks: np.ndarray, st_dist_mass: list) -> dict:
    """
    Generate lists of random values used by the synthetic cluster generating
    function.
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_isoch, N_mass = theor_tracks.shape[-1], 0
    for sdm in st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass, N_isoch)

    # Used by `move_isochrone()` and `add_errors`
    rand_norm_vals = np.random.normal(0.0, 1.0, (2, N_mass))

    # Used by `move_isochrone()`, `binarity()`
    rand_unif_vals = np.random.uniform(0.0, 1.0, (2, N_mass))

    rand_floats = {'norm': rand_norm_vals, 'unif': rand_unif_vals}
    return rand_floats


def get_lowest_mass(
    theor_tracks: np.ndarray, max_mag_syn: float, m_ini_idx: int, dm_min: float
) -> list:
    """ """
    min_masses = []
    for met_arr in theor_tracks:
        met_lst = []
        for age_arr in met_arr:
            mag, mass = age_arr[0], age_arr[m_ini_idx]
            i = np.argmin(abs(max_mag_syn - (mag + dm_min)))
            met_lst.append(mass[i])
        min_masses.append(met_lst)
    return min_masses


def remove_low_masses(self, my_cluster, dm_min: float) -> tuple[list, list]:
    """ """
    min_masses = get_lowest_mass(
        self.theor_tracks, my_cluster.cluster_dict['max_mag_syn'],
        my_cluster.cluster_dict['m_ini_idx'], dm_min)

    st_dist_mass_lmass = []
    for i, met_arr in enumerate(self.st_dist_mass):
        met_lst = []
        for j, mass_sample in enumerate(met_arr):
            min_mass = min_masses[i][j]
            msk = mass_sample > min_mass
            sampled_IMF = mass_sample[msk]
            met_lst.append(sampled_IMF)
        st_dist_mass_lmass.append(met_lst)

    # Make copy of original array, used for mass estimation in cluster() class
    st_dist_mass_full = self.st_dist_mass.copy()
    st_dist_mass = st_dist_mass_lmass
    return st_dist_mass_full, st_dist_mass


def synthcl_generate(
    self, model_fit: dict, model_fixed: dict, cluster_dict: dict, full_synth_arr=False
) -> np.ndarray:
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    * full_synth_arr = False

    synth_clust = [mag, c1, (c2)]

    where c1 and c2 colors defined.

    * full_synth_arr = True

    synth_clust = [mag, c1, (c2), m_ini_1, mag_b, c1_b, (c2_b), m_ini_2]

    where 'm_ini_1, m_ini_2' are the primary and secondary masses of the
    binary systems. The single systems only have a '0' stored in 'm_ini_2'.
    """
    # m_ini_idx = cluster_dict["N_colors"] + 1

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    met, loga, beta, av, dr, rv, dm, ml, mh, al, ah = scf.properModel(
        self.met_age_dict, model_fixed, model_fit
    )

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = scf.zaWAverage(
        self.theor_tracks,
        self.met_age_dict,
        cluster_dict['m_ini_idx'],
        met,
        loga,
        ml,
        mh,
        al,
        ah,
    )

    # Move theoretical isochrone using the distance modulus
    isoch_moved = scf.move_isochrone(isochrone, cluster_dict['m_ini_idx'], dm)

    # Apply extinction correction
    isoch_extin = scf.extinction(
        self.ext_coefs,
        self.rand_floats['norm'][0],
        self.rand_floats['unif'][0],
        self.DR_dist,
        self.DR_percentage,
        cluster_dict['m_ini_idx'],
        av,
        dr,
        rv,
        isoch_moved,
    )

    # Remove isochrone stars beyond the maximum magnitude
    isoch_cut, _ = scf.cut_max_mag(isoch_extin, cluster_dict["max_mag_syn"])
    if not isoch_cut.any():
        return np.array([])

    # Interpolate IMF's sampled masses into the isochrone.
    isoch_mass = scf.mass_interp(
        isoch_cut, cluster_dict['m_ini_idx'], self.st_dist_mass[ml][al],
        cluster_dict["N_obs_stars"])
    if not isoch_mass.any():
        return np.array([])

    # Assignment of binarity.
    isoch_binar = scf.binarity(
        self.alpha, beta, cluster_dict['m_ini_idx'],
        self.rand_floats['unif'][1], isoch_mass
    )

    # Assign errors according to errors distribution.
    synth_clust = scf.add_errors(
        isoch_binar, cluster_dict["err_lst"], self.rand_floats['norm'][1]
    )

    if full_synth_arr:
        return synth_clust
    return synth_clust[:cluster_dict['m_ini_idx']]
