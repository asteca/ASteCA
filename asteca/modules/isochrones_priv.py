import os

import numpy as np
import pandas as pd

# IMPORTANT
# These are the default column names for the initial mass, metallicity, and age for the
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


def load(
    model: str,
    isochs_path: str,
    magnitude: str,
    color: tuple,
    color2: tuple | None,
    column_names: dict | None,
    N_interp: int,
    parsec_rm_stage_9: bool,
) -> tuple[np.ndarray, list, dict, int]:
    """Load the theoretical isochrones and return processed data.

    :return: Array of isochrones, individual filters for each color defined,
        dictionary with metallicities and ages, and number of files read.
    :rtype: tuple[np.ndarray, list, dict, int]

    :raises ValueError: If there is a shape mismatch in the loaded isochrones
    """
    cols_keep, mass_col, met_col, age_col = get_columns(
        column_names, model, magnitude, color, color2
    )

    f_paths = extract_paths(isochs_path)

    # isochrones.shape = (N_photsyst, N_z, N_a, N_interp, N_cols)
    # met_age_arr.shape = (N_photsyst, N_z*N_a, 2)
    isochrones = read(
        model,
        N_interp,
        parsec_rm_stage_9,
        f_paths,
        met_col,
        age_col,
        cols_keep,
    )

    isochrones, met_age_arr = met_order_merge_massini_check(mass_col, isochrones)
    try:
        np.shape(isochrones)
    except ValueError:
        raise ValueError(
            "Shape mismatch in loaded isochrones. This usually means that an\n"
            + "incorrect number of ages and/or metallicities are stored in the\n"
            + "isochrone files."
        )

    # Create dictionary of mets and ages
    all_m, all_a = np.array(met_age_arr).T
    # Discard duplicate values
    all_m = np.array(list(dict.fromkeys(all_m))).astype(float)
    all_a = np.array(list(dict.fromkeys(all_a))).astype(float)
    met_age_dict = {"met": all_m, "loga": all_a}

    # theor_tracks.shape = (N_z, N_a, N_cols, N_interp)
    theor_tracks, color_filters = shape_isochrones(
        magnitude, color, color2, mass_col, isochrones
    )

    return theor_tracks, color_filters, met_age_dict, len(f_paths)


def get_columns(
    column_names: dict | None,
    model: str,
    magnitude: str,
    color: tuple,
    color2: tuple | None,
) -> tuple[list, str, str, str]:
    """Get the column names for the isochrones.

    :return: List of columns to keep, mass column name, metallicity column name,
        and age column name.
    :rtype: tuple[list, str, str, str]
    """
    if column_names is None:
        mass_col = phot_syst_col_names[model]["mass_col"]
        met_col = phot_syst_col_names[model]["met_col"]
        age_col = phot_syst_col_names[model]["age_col"]
    else:
        mass_col = column_names["mass_col"]
        met_col = column_names["met_col"]
        age_col = column_names["age_col"]

    # Select columns to keep
    all_filters = [magnitude] + list(color)
    if color2 is not None:
        all_filters += list(color2)
    cols_keep = list(dict.fromkeys(all_filters)) + [mass_col]

    return cols_keep, mass_col, met_col, age_col


def extract_paths(isochs_path: str) -> list:
    """Extract isochrone files from `isochs_path`.

    :return: List of isochrone file paths.
    :rtype: list
    :raises FileNotFoundError: If no files are found in the isochrones path.
    """
    # Check if path is to file or folder
    if os.path.isfile(isochs_path):
        f_paths = [isochs_path]
    else:
        f_paths = []
        # Iterate over files in directory
        for path, folders, files in os.walk(isochs_path):
            # Skip hidden folders
            if path.split("/")[-1].startswith("."):
                continue
            for filename in files:
                # Skip hidden files
                if not filename.startswith("."):
                    f_paths.append(os.path.join(path, filename))

        if len(f_paths) == 0:
            raise FileNotFoundError(
                f"No files found in isochrones path '{isochs_path}'"
            )

    return f_paths


def read(
    model: str,
    N_interp: int,
    parsec_rm_stage_9: bool,
    f_paths: list,
    met_col: str,
    age_col: str,
    cols_keep: list,
) -> dict:
    """Read isochrone files and store them in a dictionary.

    :param model: Isochrone model name.
    :type model: str
    :param N_interp: Number of points to interpolate.
    :type N_interp: int
    :param parsec_rm_stage_9: Remove post-AGB stage for PARSEC models.
    :type parsec_rm_stage_9: bool
    :param f_paths: List of isochrone file paths.
    :type f_paths: list
    :param met_col: Metallicity column name.
    :type met_col: str
    :param age_col: Age column name.
    :type age_col: str
    :param cols_keep: List of columns to keep.
    :type cols_keep: list
    :return: Dictionary of isochrones.
    :rtype: dict
    """

    group_col = {"PARSEC": met_col, "MIST": age_col, "BASTI": None}

    isochrones = {}
    for file_path in f_paths:
        # Extract columns names and full header
        col_names, full_header = get_header(file_path)

        # Columns to keep for this photometric system
        cols_keep_ps = list(set(col_names) & set(cols_keep))

        # Group file in blocks as required
        df_blocks = get_data_blocks(file_path, col_names, group_col[model])

        if model == "PARSEC":
            # Process metallicity blocks
            for _, met_df in df_blocks:
                # Group by age
                age_blocks = met_df.groupby(age_col, sort=False)
                for _, df in age_blocks:
                    zinit = df[met_col].values[0]
                    try:
                        isochrones[zinit]
                    except KeyError:
                        isochrones[zinit] = {}

                    # Remove post-AGB stage
                    if parsec_rm_stage_9 is True:
                        msk = df["label"] != "9"
                        df = df[msk]

                    age = df[age_col].values[0]
                    # Interpolated isochrone
                    isochrones = interp_df(
                        N_interp, cols_keep_ps, df, isochrones, zinit, age
                    )

        elif model == "MIST":
            zinit = get_MIST_z_val(met_col, full_header)
            try:
                isochrones[zinit]
            except KeyError:
                isochrones[zinit] = {}

            # Process age blocks
            for _, df in df_blocks:
                age = df[age_col].values[0]
                # Interpolated isochrone
                isochrones = interp_df(
                    N_interp, cols_keep_ps, df, isochrones, zinit, age
                )

        elif model == "BASTI":
            zinit, age = get_BASTI_z_a_val(full_header, met_col, age_col)
            try:
                isochrones[zinit]
            except KeyError:
                isochrones[zinit] = {}

            # Interpolated isochrone
            isochrones = interp_df(
                N_interp, cols_keep_ps, df_blocks, isochrones, zinit, age
            )

    return isochrones


def get_header(file_path: str) -> tuple[list, list]:
    """Iterate through each line in the file to get the header.

    :param file_path: Path to the isochrone file.
    :type file_path: str
    :return: List of column names and full header.
    :rtype: tuple[list, list]
    :raises ValueError: If the header cannot be extracted from the file.
    """
    with open(file_path, mode="r") as f_iso:
        full_header, column_names = [], ""
        for i, line in enumerate(f_iso):
            # Skip BASTI lines that look like this
            if line.startswith("#====="):
                continue
            full_header.append(line)
            if not line.startswith("#"):
                break
            column_names = line

        if column_names == "":
            raise ValueError(f"Could not extract header from file: {file_path}")

        column_names = column_names.replace("#", "").split()

    return column_names, full_header


def get_data_blocks(
    file_path: str, header: list, block_col: str | None = None
) -> pd.DataFrame | pd.core.groupby.generic.DataFrameGroupBy:
    """Extract data in 'block_col' blocks.

    :param file_path: Path to the isochrone file.
    :type file_path: str
    :param header: List of column names.
    :type header: list
    :param block_col: Column to group by, defaults to None
    :type block_col: str | None, optional
    :return: DataFrame or DataFrameGroupBy object.
    :rtype: pd.DataFrame | pd.core.groupby.generic.DataFrameGroupBy
    """

    # IMPORTANT
    # We are assuming that all models use spaces to separate columns. If this ever
    # changes, this line will no longer work
    df = pd.read_csv(file_path, comment="#", header=None, names=header, sep=r"\s+")

    if block_col is None:
        return df

    grouped = df.groupby(block_col, sort=False)
    return grouped


def get_MIST_z_val(met_col: str, full_header: list) -> str:
    """Extract metallicity value for MIST files.

    :param met_col: Metallicity column name.
    :type met_col: str
    :param full_header: Full header of the isochrone file.
    :type full_header: list
    :return: Metallicity value.
    :rtype: str
    :raises ValueError: If the metallicity cannot be read from the header.
    """
    z_header_vals, z_header_cols = None, None
    for i, line in enumerate(full_header):
        if met_col in line:
            z_header_cols = line.split()[1:]
            z_header_vals = full_header[i + 1].split()[1:]

    if z_header_vals is None or z_header_cols is None:
        raise ValueError("Could not read header from MIST isochrone")

    z_df = pd.DataFrame([z_header_vals], columns=z_header_cols)
    zinit = z_df[met_col].values
    return str(zinit[0])


def get_BASTI_z_a_val(full_header: list, met_col: str, age_col: str) -> tuple[str, str]:
    """Extract metallicity and age values for BASTI files.

    :param full_header: Full header of the isochrone file.
    :type full_header: list
    :param met_col: Metallicity column name.
    :type met_col: str
    :param age_col: Age column name.
    :type age_col: str
    :return: Metallicity and age values.
    :rtype: tuple[str, str]
    :raises ValueError: If the metallicity or age cannot be read from the header.
    """
    zval, aval = None, None
    for line in full_header:
        if met_col in line:
            spl_1 = line.split(met_col)[1].replace("=", " ")
            zval = spl_1.split()[0]

            spl_2 = line.split(age_col)[1].replace("=", " ").strip()
            aval = np.log10(float(spl_2) * 1e6)
            aval = str(round(aval, 5))
            break

    if zval is None or aval is None:
        raise ValueError("Could not read header from Basti isochrone")

    return zval, aval


def interp_df(
    N_interp: int,
    cols_keep_ps: list,
    df: pd.DataFrame,
    isochrones: dict,
    zinit: str,
    age: str,
) -> dict:
    """Interpolate the isochrone data.

    :param N_interp: Number of points to interpolate.
    :type N_interp: int
    :param cols_keep_ps: List of columns to keep.
    :type cols_keep_ps: list
    :param df: DataFrame with the isochrone data.
    :type df: pd.DataFrame
    :param isochrones: Dictionary of isochrones.
    :type isochrones: dict
    :param zinit: Metallicity value.
    :type zinit: str
    :param age: Age value.
    :type age: str
    :return: Dictionary of isochrones.
    :rtype: dict
    """

    # Drop columns not in list
    df = df[df.columns.intersection(cols_keep_ps)]
    # Dataframe to floats
    df = df.apply(pd.to_numeric)

    # Only interpolate if there are extra points to add
    if len(df) >= N_interp:
        isoch_interp = df
    else:
        # Interpolate
        xx = np.linspace(0.0, 1.0, N_interp)

        # # Works but the binary sequence is really affected...
        # N_interp = 1000
        # N1, N2, N3, N4 = int(.1 * N_interp), int(.2 * N_interp), int(.3 * N_interp), int(.4 * N_interp)
        # xx = np.array(list(
        #     np.linspace(.0, .25, N1))
        #     + list(np.linspace(.25, .5, N2))
        #     + list(np.linspace(.5, .75, N3))
        #     + list(np.linspace(.75, 1, N4))
        # )

        xp = np.linspace(0.0, 1.0, len(df))

        df_new = {}
        for col in cols_keep_ps:
            df_new[col] = np.interp(xx, xp, df[col])
        isoch_interp = pd.DataFrame(df_new)

    # Add to the dictionary of isochrones
    try:
        isochrones[zinit][age]
    except KeyError:
        isochrones[zinit][age] = []
    isochrones[zinit][age].append(isoch_interp)

    return isochrones


def met_order_merge_massini_check(mass_col: str, isochrones: dict) -> tuple[list, list]:
    """Order by metallicity, combine photometric systems, and check initial masses.

    1. Order by metallicity.
    2. Combine photometric systems if more than one was used.
    3. Check that initial masses are equal across all systems

    The isochrone parameter 'mass_col' is assumed to be equal across
    photometric systems, for a given metallicity and age. We check here that
    this is the case.

    :param mass_col: Name of the mass column.
    :type mass_col: str
    :param isochrones: Dictionary of isochrones.
    :type isochrones: dict
    :return: Ordered isochrones and metallicity-age array.
    :rtype: tuple[list, list]
    :raises ValueError: If initial mass values differ across photometric systems.
    """
    # Metallicities in order
    idx_met_sort = np.argsort(list(map(float, isochrones.keys())))
    all_mets = np.array(list(isochrones.keys()))[idx_met_sort]

    met_age_arr, isochrones_ordered = [], []
    for met_k in all_mets:
        met_vals = []
        for age_k, age_list in isochrones[met_k].items():
            met_age_arr.append([met_k, age_k])

            df_age_0 = age_list[0]
            # Check masses if more than one photometric system was used
            if len(age_list) > 1:
                for df_age_X in age_list[1:]:
                    mass_diff = abs(
                        df_age_0[mass_col].values - df_age_X[mass_col].values
                    ).sum()
                    if mass_diff > 0.001:
                        raise ValueError(
                            "Initial mass values differ across photometric systems"
                        )
                    # Drop 'mass_col' column from this photometric system
                    df_age_X = df_age_X.drop([mass_col], axis=1)
                    # Merge data frames and replace original one
                    df_age_0 = pd.concat([df_age_0, df_age_X], axis=1)

            met_vals.append(df_age_0)
        isochrones_ordered.append(met_vals)

    return isochrones_ordered, met_age_arr


def shape_isochrones(
    magnitude: str,
    color: tuple,
    color2: tuple | None,
    initial_mass: str,
    isochrones: list,
) -> tuple[np.ndarray, list]:
    """Reshape the isochrones array.

    isochrones.shape = Nz, Na, Ni, Nd
    Nz: number of metallicities
    Na: number of log(age)s
    Ni: number of interpolated values
    Nd: number of data columns --> mini + mag + 2 * colors (one mag per color)

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

    :param initial_mass: Name of the initial mass column.
    :type initial_mass: str
    :param isochrones: List of isochrones.
    :type isochrones: list
    :return: Reshaped isochrones array and list of color filters.
    :rtype: tuple[np.ndarray, list]
    """

    Nz, Na, Ni, Nd = np.shape(isochrones)

    all_colors = [color]
    if color2 is not None:
        all_colors.append(color2)
    N_colors = len(all_colors)

    # Array that will store all the interpolated tracks
    theor_tracks = np.zeros([Nz, Na, (1 + 1 + N_colors), Ni])

    # Store the magnitudes to generate the colors separately. Used only by the
    # binarity process
    color_filters = []

    for i, met in enumerate(isochrones):
        met_lst = []
        for j, df_age in enumerate(met):
            # Store magnitude
            theor_tracks[i][j][0] = df_age[magnitude]
            # Store colors
            cols_dict = {}
            for k, color_filts in enumerate(all_colors):
                f1 = df_age[color_filts[0]]
                f2 = df_age[color_filts[1]]
                # Store color
                theor_tracks[i][j][1 + k] = f1 - f2

                # Individual filters for colors, used for binarity
                cols_dict[color_filts[0]] = f1
                cols_dict[color_filts[1]] = f2

            met_lst.append(cols_dict)

            theor_tracks[i][j][k + 2] = df_age[initial_mass]  # TODO: k value?
        color_filters.append(met_lst)

    return theor_tracks, color_filters
