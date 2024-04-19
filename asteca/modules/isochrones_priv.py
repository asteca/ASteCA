import os
import numpy as np
import pandas as pd


# These are the column names for the initial mass, metallicity, and age for the
# supported photometric system's isochrones files.
phot_syst_col_names = {
    "PARSEC": {"mass_col": "Mini", "met_col": "Zini", "age_col": "logAge"},
    "MIST": {
        "mass_col": "initial_mass",
        "met_col": "Zinit",
        "age_col": "log10_isochrone_age_yr",
    },
    "BASTI": {"mass_col": "M/Mo(ini)", "met_col": "Z =", "age_col": "Age (Myr) ="},
}


def load(self) -> tuple[np.ndarray, list, dict]:
    r"""Load the theoretical isochrones and return a dictionary with the data.

    Returns
    -------
    theor_tracks: np.ndarray
        Array of isochrones. Updated by `synth_clusters.synthcl_load()`
    color_filters : list
        Individual filters for each color defined. Used by
        `synth_clusters.synthcl_load.add_binarity()`
    met_age_dict : dict
        Used by `synth_clusters.synthcl_generate()`

    """
    cols_keep, mass_col, met_col, age_col = get_columns(self)

    f_paths = extract_paths(self)

    # isochrones.shape = (N_photsyst, N_z, N_a, N_interp, N_cols)
    # met_age_arr.shape = (N_photsyst, N_z*N_a, 2)
    isochrones, met_age_arr = read(
        self.model,
        self.N_interp,
        self.parsec_rm_stage_9,
        f_paths,
        met_col,
        age_col,
        cols_keep,
    )
    N_ps, Nz, Na = len(isochrones), len(isochrones[0]), len(isochrones[0][0])
    print(
        f"Processing {self.model} isochrones for {N_ps} photometric systems, "
        + f"N_met={Nz}, N_age={Na}..."
    )

    isochrones = m_ini_check(self, mass_col, isochrones)
    met_age_dict = met_ages_check(self, met_age_arr)

    # theor_tracks.shape = (N_z, N_a, N_cols, N_interp)
    theor_tracks, color_filters = order_isochrones(self, mass_col, isochrones)

    return theor_tracks, color_filters, met_age_dict


def get_columns(self) -> tuple[list, str, str, str]:
    """ """
    if self.column_names is None:
        mass_col = phot_syst_col_names[self.model]["mass_col"]
        met_col = phot_syst_col_names[self.model]["met_col"]
        age_col = phot_syst_col_names[self.model]["age_col"]
    else:
        mass_col = self.column_names["mass_col"]
        met_col = self.column_names["met_col"]
        age_col = self.column_names["age_col"]

    # Select columns to keep
    cols_keep = list(self.mag_color_lambdas.keys()) + [mass_col]

    return cols_keep, mass_col, met_col, age_col


def extract_paths(self) -> list:
    """Extract isochrone files from `isochs_path`.
    """
    f_paths = []
    # For each photometric system
    for ps in os.listdir(self.isochs_path):
        ps_path = os.path.join(self.isochs_path, ps)

        # Only folders beyond this point
        if os.path.isdir(ps_path) is False:
            continue

        if self.model == "PARSEC":
            # Remove possible hidden files inside this folder
            in_folder = [_ for _ in os.listdir(ps_path) if not _.startswith('.')]
            # A single file per photometric system is expected here
            if len(in_folder) > 1:
                raise ValueError(
                    "One file per photometric system is expected for the PARSEC service"
                )
            f_paths.append(os.path.join(ps_path, in_folder[0]))

        elif self.model == "MIST":
            met_files = []
            for met in os.listdir(ps_path):
                if not met.startswith('.'):
                    met_files.append(os.path.join(ps_path, met))
            f_paths.append(met_files)

        elif self.model == "BASTI":
            met_folders = [_ for _ in os.listdir(ps_path) if not _.startswith('.')]
            met_files = []
            for met_folder in met_folders:
                met_folder_path = os.path.join(ps_path, met_folder)
                met_lst = []
                for met in os.listdir(met_folder_path):
                    if not met.startswith('.'):
                        met_lst.append(os.path.join(met_folder_path, met))
                met_files.append(met_lst)
            f_paths.append(met_files)

    # Check that the number of files matches across photometric systems and
    # metallicities for BASTI
    if self.model == "MIST":
        nmets = []
        for ps in f_paths:
            nmets.append(len(ps))
        if len(list(set(nmets))) > 1:
            raise ValueError(
                "The number of files for each photometric system must match"
            )
    elif self.model == "BASTI":
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


def m_ini_check(self, initial_mass: str, isochrones: list) -> list:
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


def met_ages_check(self, met_age_arr: list) -> dict:
    """ """
    phot_syst_0 = met_age_arr[0]
    for phot_syst in met_age_arr[1:]:
        for i, met_age in enumerate(phot_syst_0):
            met_age_x = phot_syst[i]

            if abs(np.array(met_age) - np.array(met_age_x)).sum() > 0.001:
                raise ValueError(
                    "Metallicities and ages are not equal across photometric systems"
                )

        phot_syst_0 = phot_syst

    # Discard duplicate values
    all_z, all_a = np.array(phot_syst_0).T
    all_z = np.array(list(dict.fromkeys(all_z))).astype(float)
    all_a = np.array(list(dict.fromkeys(all_a))).astype(float)

    met_age_dict = {"met": all_z, "a": all_a}

    return met_age_dict


def order_isochrones(
    self, initial_mass: str, isochrones: list
) -> tuple[np.ndarray, list]:
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
    N_colors = len(self.color_filter_name)

    # Array that will store all the interpolated tracks
    theor_tracks = np.zeros([Nz, Na, (1 + 1 + N_colors), Ni])

    # Store the magnitudes to generate the colors separately. Used only by the
    # binarity process
    color_filters = []

    for i, met in enumerate(isochrones):
        met_lst = []
        for j, df_age in enumerate(met):
            # Store magnitude
            theor_tracks[i][j][0] = df_age[self.mag_filter_name]
            # Store colors
            cols_dict = {}
            for k, color_filts in enumerate(self.color_filter_name):
                f1 = df_age[color_filts[0]]
                f2 = df_age[color_filts[1]]
                # Store color
                theor_tracks[i][j][1 + k] = f1 - f2

                # Individual filters for colors, used for binarity
                cols_dict[color_filts[0]] = f1
                cols_dict[color_filts[1]] = f2

            met_lst.append(cols_dict)

            theor_tracks[i][j][k + 2] = df_age[initial_mass]
        color_filters.append(met_lst)

    return theor_tracks, color_filters


def read(
    model, N_interp, parsec_rm_stage_9, f_paths, met_col, age_col, cols_keep
) -> tuple[list, list]:
    """ """

    if model == "PARSEC":
        isochrones, met_age_arr = read_PARSEC_files(
            f_paths, parsec_rm_stage_9, met_col, age_col, cols_keep, N_interp
        )
    elif model == "MIST":
        isochrones, met_age_arr = read_MIST_files(
            f_paths, met_col, age_col, cols_keep, N_interp
        )
    elif model == "BASTI":
        isochrones, met_age_arr = read_BASTI_files(
            f_paths, met_col, age_col, cols_keep, N_interp
        )

    return isochrones, met_age_arr


def read_PARSEC_files(
    f_paths, parsec_rm_stage_9, met_col, age_col, cols_keep, N_interp, label_col="label"
):
    """
    Each PARSEC file represents a single photometric system.
    """
    isochrones, met_age_arr = [], []
    for isoch_path in f_paths:
        isochs_array, met_age = [], []
        with open(isoch_path, mode="r") as f_iso:
            header = get_header(f_iso)

            # Columns to keep for this photometric system
            cols_keep_ps = list(set(header) & set(cols_keep))

            # Group by metallicity
            met_blocks = get_data_blocks(f_iso, header, met_col)
            met_arr = []
            for _, met_df in met_blocks:
                # Group by age
                a_grouped = met_df.groupby(age_col, sort=False)
                age_arr = []
                for _, df in a_grouped:
                    # Remove post-AGB stage
                    if parsec_rm_stage_9 is True:
                        msk = df["label"] != "9"
                        df = df[msk]

                    # Save (z, a) values
                    met_age.append(
                        [float(df[met_col].values[0]), float(df[age_col].values[0])]
                    )

                    # Interpolate single isochrone
                    age_arr.append(interp_df(N_interp, cols_keep_ps, df))

                # Store single isochrone for this metallicity
                met_arr.append(age_arr)

            # Store all isochrones for this metallicity
            isochs_array += met_arr

        isochrones.append(isochs_array)
        met_age_arr.append(met_age)

    return isochrones, met_age_arr


def read_MIST_files(f_paths, met_col, age_col, cols_keep, N_interp):
    """
    MIST isochrones are organized in one single file per metallicity. Multiple
    files can mean multiple metallicities or multiple photometric systems or both.
    A single file can also mean multiple photometric systems (e.g.: Gaia, UBVRI, etc)
    """

    def get_z_val(f_iso):
        """
        Extract metallicity value for MIST files.
        """
        for i, line in enumerate(f_iso):
            next_line = False
            if met_col in line:
                z_header_cols = line.split()[1:]
                next_line = True
            if next_line:
                z_header_vals = next(f_iso).split()[1:]
                break
        z_df = pd.DataFrame([z_header_vals], columns=z_header_cols)
        zinit = float(z_df[met_col].values)
        return zinit

    # Sort files according to their metallicities BEFORE processing them
    f_paths_zorder = []
    for phot_syst in f_paths:
        zinit_vals = []
        for file_path in phot_syst:
            with open(file_path, mode="r") as f_iso:
                zinit = get_z_val(f_iso)
                zinit_vals.append(zinit)
        z_order = np.argsort(zinit_vals)
        f_paths_zorder.append(list(np.array(phot_syst)[z_order]))
    f_paths = f_paths_zorder

    #
    isochrones, met_age_arr = [], []
    for phot_syst in f_paths:
        isochrones_temp, met_age_temp = [], []
        for file_path in phot_syst:
            isochs_array = []
            # Open the tracks file.
            with open(file_path, mode="r") as f_iso:
                zinit = get_z_val(f_iso)
                header = get_header(f_iso)

                # Columns to keep for this photometric system
                cols_keep_ps = list(set(header) & set(cols_keep))

                age_blocks = get_data_blocks(f_iso, header, age_col)
                met_arr = []
                for name, df in age_blocks:
                    met_age_temp.append([zinit, float(df[age_col].values[0])])
                    # Store single isochrone for this metallicity
                    met_arr.append(interp_df(N_interp, cols_keep_ps, df))

                # Store all isochrones for this metallicity
                isochs_array += met_arr

            isochrones_temp.append(isochs_array)

        isochrones.append(isochrones_temp)
        met_age_arr.append(met_age_temp)

    return isochrones, met_age_arr


def read_BASTI_files(f_paths, met_col, age_col, cols_keep, N_interp):
    """
    BASTI isochrones are organized in one single file per metallicity per age.
    """

    def get_z_a_val(f_iso, met_col, age_col):
        """
        Extract metallicity and age values for BASTI files.
        """
        for i, line in enumerate(f_iso):
            if met_col in line:
                spl_1 = line.split(met_col)[1].replace("=", " ")
                z_header_val = spl_1.split()[0]
                zval = float(z_header_val)

                spl_2 = line.split(age_col)[1].replace("=", " ").strip()
                aval = np.log10(float(spl_2) * 1e6)
                break
        return zval, aval

    # Sort files according to their metallicities and aged BEFORE processing them
    f_paths_zaorder = []
    for phot_syst in f_paths:
        za_vals = []
        for file_path in phot_syst:
            with open(file_path, mode="r") as f_iso:
                zval, aval = get_z_a_val(f_iso, met_col, age_col)
                za_vals.append([file_path, zval, aval])
        # Sort by z then age
        f_paths_zaorder.append(list(sorted(za_vals, key=lambda el: (el[1], el[2]))))

    #
    isochrones, met_age_arr = [], []
    for phot_syst in f_paths_zaorder:
        # Extract initial z value
        z_old = phot_syst[0][1]

        isochs_array, met_arr, met_age_temp = [], [], []
        for file_path in phot_syst:
            zval, aval = file_path[1], file_path[2]
            met_age_temp.append([zval, aval])

            # Open the tracks file.
            with open(file_path[0], mode="r") as f_iso:
                header = get_header(f_iso, True)
                # Columns to keep for this photometric system
                cols_keep_ps = list(set(header) & set(cols_keep))

                df = get_data_blocks(f_iso, header, age_col, True)
                df_interp = interp_df(N_interp, cols_keep_ps, df)

                # Store single isochrone for this metallicity
                if zval == z_old:
                    met_arr.append(df_interp)
                else:
                    isochs_array.append(met_arr)
                    met_arr = [df_interp]

            z_old = zval
        # Store final metallicity
        isochs_array.append(met_arr)

        isochrones.append(isochs_array)
        met_age_arr.append(met_age_temp)

    return isochrones, met_age_arr


def get_header(f_iso, BASTI_flag=False):
    """Iterate through each line in the file to get the header"""
    all_lines = []
    for i, line in enumerate(f_iso):
        all_lines.append(line)
        if not line.startswith("#"):
            break
        header = line

    if BASTI_flag is True:
        header = all_lines[i - 2].replace("#", "").split()
    else:
        header = header.replace("#", "").split()
    return header


def get_data_blocks(f_iso, header, age_col, BASTI_flag=False):
    """ """
    f_iso.seek(0, 0)
    temp_f = []
    for i, line in enumerate(f_iso):
        if line.startswith("#") or line.startswith("\n"):
            pass
        else:
            temp_f.append(line.split())
    df = pd.DataFrame(temp_f, columns=header)
    if BASTI_flag is True:
        return df

    grouped = df.groupby(age_col, sort=False)
    return grouped


def interp_df(N_interp, cols_keep_ps, df):
    """ """
    # Drop columns not in list
    df = df[df.columns.intersection(cols_keep_ps)]
    # Dataframe to floats
    df = df.apply(pd.to_numeric)
    # Interpolate
    xx = np.linspace(0.0, 1.0, N_interp)
    xp = np.linspace(0.0, 1.0, len(df))
    df_new = {}
    for col in cols_keep_ps:
        df_new[col] = np.interp(xx, xp, df[col])
    df_interp = pd.DataFrame(df_new)
    return df_interp
