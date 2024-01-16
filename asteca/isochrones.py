import os
import numpy as np
import pandas as pd
from . import isochs_read


def isochs_load(self):
    """

    isochrones.shape = (N_photsyst, N_z, N_a, N_interp, N_cols)
    met_age_arr.shape = (N_photsyst, N_z*N_a, 2)

    theor_tracks.shape = (N_z, N_a, N_cols, N_interp)
    """

    # Read and interpolate all metallicities an ages
    print(f"Reading {self.model} isochrones...")

    cols_keep, met_col, age_col, mass_col = get_columns(self, self.model)

    f_paths = extract_paths(self.isochs_path, self.model)

    if self.model == "PARSEC":
        isochrones, met_age_arr = isochs_read.read_PARSEC_files(
            f_paths, met_col, age_col, cols_keep, self.N_interp
        )
    elif self.model == "MIST":
        isochrones, met_age_arr = isochs_read.read_MIST_files(
            f_paths, met_col, age_col, cols_keep, self.N_interp
        )
    elif self.model == "BASTI":
        isochrones, met_age_arr = isochs_read.read_BASTI_files(
            f_paths, met_col, age_col, cols_keep, self.N_interp
        )

    N_ps, Nz, Na = len(isochrones), len(isochrones[0]), len(isochrones[0][0])
    print(f"Processing isochrones for {N_ps} photometric systems, N_z={Nz}, N_a={Na}")
    isochrones = m_ini_check(mass_col, isochrones)
    met_age_dict = met_ages_check(met_age_arr)

    theor_tracks, mags_cols_intp = order_isochrones(
        self.mag_filter_name, self.color_filter_name, mass_col, isochrones
    )

    isochs_dict = {
        "theor_tracks": theor_tracks,
        "met_age_dict": met_age_dict,
        "mags_cols_intp": mags_cols_intp,
        "filter_lambdas": self.mag_color_lambdas,
        "mag_filter_name": self.mag_filter_name,
        "color_filter_name": self.color_filter_name,
    }

    return isochs_dict


def get_columns(self, model):
    """ """
    mass_col = self.initial_mass_col
    met_col = self.metallicity_col
    age_col = self.age_col

    if model == "PARSEC":
        if mass_col is None:
            mass_col = "Mini"
        if met_col is None:
            met_col = "Zini"
        if age_col is None:
            age_col = "logAge"
    elif model == "MIST":
        if mass_col is None:
            mass_col = "initial_mass"
        if met_col is None:
            met_col = "Zinit"
        if age_col is None:
            age_col = "log10_isochrone_age_yr"
    elif model == "BASTI":
        if mass_col is None:
            mass_col = "M/Mo(ini)"
        if met_col is None:
            met_col = "Z ="
        if age_col is None:
            age_col = "Age (Myr) ="

    # Select columns to keep
    cols_keep = list(self.mag_color_lambdas.keys()) + [mass_col]

    return cols_keep, met_col, age_col, mass_col


def extract_paths(isochs_path, model):
    """
    Extract files from 'isochs_path'.

    What is in a file for each model (for a given photometric system):

    PARSEC: multiple metallicities and multiple ages
    MIST  : single metallicity and multiple ages
    BASTI : single metallicity and single age

    > : folder
    >>: file

    > PARSEC
         |---> phot_syst_1
                   |-------->> mets_ages
         |---> phot_syst_2
                   |-------->> ...

    > MIST
         |---> phot_syst_1
                   |-------->> met_1_ages
                   |-------->> met_2_ages
                   |-------->> ...
         |---> phot_syst_2
                   |-------->> ...

    > BASTI
         |---> phot_syst_1
                   |--------> met_1
                              |---->> age_1
                              |---->> age_2
                              |---->> ...
                   |--------> met_2
                              |---->> age_1
                              |---->> ...
         |---> phot_syst_2
                   |--------> ...

    """
    f_paths = []
    # For each photometric system
    for ps in os.listdir(isochs_path):
        ps_path = os.path.join(isochs_path, ps)

        if model == "PARSEC":
            # A single file per photometric system is expected here
            if len(os.listdir(ps_path)) > 1:
                raise ValueError("One file per photometric system is expected")
            for met in os.listdir(ps_path):
                f_paths.append(os.path.join(ps_path, met))

        elif model == "MIST":
            met_files = []
            for met in os.listdir(ps_path):
                met_files.append(os.path.join(ps_path, met))
            f_paths.append(met_files)

        elif model == "BASTI":
            met_files = []
            for met_folder in os.listdir(ps_path):
                met_folder_path = os.path.join(ps_path, met_folder)
                met_lst = []
                for met in os.listdir(met_folder_path):
                    met_lst.append(os.path.join(met_folder_path, met))
                met_files.append(met_lst)
            f_paths.append(met_files)

    # Check that the number of files matches across photometric systems and
    # metallicities for BASTI
    if model == "MIST":
        nmets = []
        for ps in f_paths:
            nmets.append(len(ps))
        if len(list(set(nmets))) > 1:
            raise ValueError(
                "The number of files for each photometric system must match"
            )
    elif model == "BASTI":
        nmets = []
        for ps in f_paths:
            nmets.append(len(ps))
        if len(list(set(nmets))) > 1:
            raise ValueError(
                "The number of files for each photometric system must match"
            )
        # Now check the number of age files for each metallicity
        for ps in f_paths:
            nages = []
            for met in ps:
                nages.append(len(met))
            if len(list(set(nages))) > 1:
                raise ValueError(
                    "The number of files for each metallicity system must match"
                )

        # Flatten the list once checked
        f_paths_flat = []
        for ps in f_paths:
            met_f_ps = []
            for met_f in ps:
                met_f_ps += met_f
            f_paths_flat.append(met_f_ps)
        f_paths = f_paths_flat

    return f_paths


def m_ini_check(initial_mass, isochrones):
    """
    The isochrone parameter 'initial_mass' is assumed to be equal across
    photometric systems, for a given metallicity and age. We check here that
    this is the case.

    Combine photometric systems if more than one was used.
    """
    if len(isochrones) == 1:
        isochrones_merged = isochrones[0]
        return isochrones_merged

    isochrones_merged = []
    phot_syst_0 = isochrones[0]
    for phot_syst in isochrones[1:]:
        for i, met in enumerate(phot_syst_0):
            met_vals = []
            for j, df_age in enumerate(met):
                df_x = phot_syst[i][j]

                if (
                    abs(df_age[initial_mass].values - df_x[initial_mass].values).sum()
                    > 0.001
                ):
                    raise ValueError(
                        "Initial mass values are not equal across photometric systems"
                    )

                # Merge and drop a mass_ini column
                df_x = df_x.drop(["Mini"], axis=1)
                met_vals.append(pd.concat([df_age, df_x], axis=1))
            isochrones_merged.append(met_vals)

        phot_syst_0 = phot_syst

    return isochrones_merged


def met_ages_check(met_age_dict):
    """ """
    phot_syst_0 = met_age_dict[0]
    for phot_syst in met_age_dict[1:]:
        for i, met_age in enumerate(phot_syst_0):
            met_age_x = phot_syst[i]

            if abs(np.array(met_age) - np.array(met_age_x)).sum() > 0.001:
                raise ValueError(
                    "Metallicities and ages are not equal across photometric systems"
                )

        phot_syst_0 = phot_syst

    # Select a single array
    met_age_dict = phot_syst_0

    # Discard duplicate values
    all_z, all_a = np.array(met_age_dict).T
    all_z = np.array(list(dict.fromkeys(all_z))).astype(float)
    all_a = np.array(list(dict.fromkeys(all_a))).astype(float)

    # Store in dictionary to use later
    met_age_dict = {"z": all_z, "a": all_a}

    return met_age_dict


def order_isochrones(mag_filter_name, color_filter_name, initial_mass, isochrones):
    """
    Return list structured as:

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f,.., c1, c2,.., Mini, fb,.., c1b, c2b,.., Minib]

    where:
    N     : number of metallicities in grid
    M     : number of ages in grid
    f     : magnitude
    cX    : colors
    Mini  : initial mass
    fb    : binary magnitude --
    cXb   : binary colors     |
    Minib : binary masses -----

    theor_tracks.shape = (Nz, Na, Nd, Ni)
    Nz: number of metallicities
    Na: number of log(age)s
    Nd: number of data columns
    Ni: number of interpolated values
    """

    # isochrones.shape = Nz, Na, Ni, Nd
    # Nz: number of metallicities
    # Na: number of log(age)s
    # Ni: number of interpolated values
    # Nd: number of data columns --> mini + mag + 2 * colors (one mag per color)
    Nz, Na, Ni, Nd = np.shape(isochrones)

    # m_ini_col = isoch_data['initial_mass']
    # mag_filter = cluster_data['mag_filter_name']
    # col_filter = cluster_data['color_filter_name']
    N_colors = len(color_filter_name)

    # Array that will store all the interpolated tracks
    theor_tracks = np.zeros([Nz, Na, (1 + 1 + N_colors), Ni])
    # Store the magnitudes to generate the colors separately. Used only by the
    # binarity process
    mags_cols_intp = []

    for i, met in enumerate(isochrones):
        met_lst = []
        for j, df_age in enumerate(met):
            # Store magnitude
            theor_tracks[i][j][0] = df_age[mag_filter_name]
            # Store colors
            cols_dict = {}
            for k, color_filts in enumerate(color_filter_name):
                f1 = df_age[color_filts[0]]
                f2 = df_age[color_filts[1]]
                # Store color
                theor_tracks[i][j][1 + k] = f1 - f2

                # Individual filters for colors, used for binarity
                cols_dict[color_filts[0]] = f1
                cols_dict[color_filts[1]] = f2

            met_lst.append(cols_dict)

            theor_tracks[i][j][k + 2] = df_age[initial_mass]
        mags_cols_intp.append(met_lst)

    return theor_tracks, mags_cols_intp


def isochs_minmax(isochs_dict):
    """ """
    zmin = isochs_dict['met_age_dict']['z'].min()
    zmax = isochs_dict['met_age_dict']['z'].max()
    amin = isochs_dict['met_age_dict']['a'].min()
    amax = isochs_dict['met_age_dict']['a'].max()
    return zmin, zmax, amin, amax
