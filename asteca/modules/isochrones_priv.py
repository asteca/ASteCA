import os
import warnings

import numpy as np

# IMPORTANT: this dictionary contains fixed data required to read isochrones from
# each of the models. If this ever changes, this function will no longer work
#
# col_names: default column names for the initial mass, metallicity, and age
# comment_char: character used to indicate comments in the isochrone files
# idx_col_line: index of the line that contains the column names
# myr_to_log10: flag to indicate whether to convert Myr to log10

# DEPRECATED 17/03/25. When I replaced pandas with np.genfromtxt this parameter
# became obsolete because np.genfromtxt does not accept regex. By default it splits
# columns using empty spaces, which works for now at least.
#
# sep_cols: regular expression to separate columns in the isochrone files

phot_systs_data = {
    "PARSEC": {
        "col_names": {
            "eep_col": None,
            "mass_col": "Mini",
            "met_col": "Zini",
            "age_col": "logAge",
        },
        "comment_char": "#",
        # "sep_cols": r"\s+",
        "idx_col_line": -1,
        "parsec_stage_9_col": "label",  # column name for the post-AGB stage
        "parsec_stage_9_id": "9",  # ID for the post-AGB stage
    },
    "MIST": {
        "col_names": {
            "eep_col": "EEP",
            "mass_col": "initial_mass",
            "met_col": "Zinit",
            "age_col": "log10_isochrone_age_yr",
        },
        "comment_char": "#",
        # "sep_cols": r"\s+",
        "idx_col_line": -1,
    },
    "BASTI": {
        "col_names": {
            "eep_col": None,
            "mass_col": "M/Mo(ini)",
            "met_col": "Z =",
            "age_col": "Age (Myr) =",
        },
        "comment_char": "#",
        # "sep_cols": r"\s+",
        "idx_col_line": -2,
        "myr_to_log10": True,
    },
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

    :param model: The model must be one of the supported isochrone services (uppercase).
    :type model: str
    :param isochs_path: Path to the isochrone files.
    :type isochs_path: str
    :param magnitude: Magnitude filter to load.
    :type magnitude: str
    :param color: Tuple with the two filters to generate the first color.
    :type color: tuple
    :param color2: Tuple with the two filters to generate the second color,
        defaults to None
    :type color2: tuple | None
    :param column_names: Dictionary with the column names for the isochrones,
        defaults to None
    :type column_names: dict | None
    :param N_interp: Number of points to interpolate.
    :type N_interp: int
    :param parsec_rm_stage_9: Flag to indicate whether to remove post-AGB stage for
        PARSEC models, defaults to True
    :type parsec_rm_stage_9: bool

    :return: Array of isochrones, individual filters for each color defined,
        dictionary with metallicities and ages, and number of files read.
    :rtype: tuple[np.ndarray, list, dict, int]
    """
    all_filters, eep_col, mass_col, met_col, age_col = get_columns(
        column_names, model, magnitude, color, color2
    )

    f_paths = extract_paths(isochs_path)

    # isochrones.shape = (N_photsyst, N_z, N_a, N_interp, N_cols)
    # met_age_arr.shape = (N_photsyst, N_z*N_a, 2)
    met_age_vals, isoch_arrays = read(
        model,
        parsec_rm_stage_9,
        f_paths,
        all_filters,
        eep_col,
        mass_col,
        met_col,
        age_col,
    )

    met_age_sorted, isochs_sorted = sort_by_met_age(met_age_vals, isoch_arrays)

    # N_interp = 10_000 gives excellent photometric results
    # N_interp = 5000 gives good photometric results
    # N_interp = 2500 gives reasonable photometric results
    # N_interp = 2000 gives borderline photometric results (binary artifacts in low mass region)
    # N_interp = 1000 gives sub-par photometric results
    # N_interp = 500 gives poor photometric results
    # Time performance:
    # 10_000: 1, 5000: 0.74, 2500: 0.63, 2000: 0.59, 1000: 0.53, 500: 0.51
    isochrones = interp_isochrones(
        N_interp, mass_col, eep_col, met_age_sorted, isochs_sorted
    )

    isochrones, met_age_arr = merge_ps_massini_check(mass_col, isochrones)

    # Create dictionary of mets and ages
    all_m, all_a = np.array(met_age_arr).T
    # Discard duplicate values
    all_m = np.array(list(dict.fromkeys(all_m))).astype(float)
    all_a = np.array(list(dict.fromkeys(all_a))).astype(float)
    met_age_dict = {"met": all_m, "loga": all_a}

    if max(met_age_dict["loga"][1:] - met_age_dict["loga"][:-1]) > 0.25:
        warnings.warn(
            "The isochrones grid has very large steps in log(age). This can impact "
            + "the proper generation of synthetic clusters. Consider using a smaller step."
        )

    # theor_tracks.shape = (N_z, N_a, N_cols, N_interp)
    theor_tracks, color_filters = reshape_isochrones(
        magnitude, color, color2, mass_col, isochrones
    )

    return theor_tracks, color_filters, met_age_dict, len(f_paths)


def get_columns(
    column_names: dict | None,
    model: str,
    magnitude: str,
    color: tuple,
    color2: tuple | None,
) -> tuple[list, str | None, str, str, str]:
    """Get the column names for the isochrones.

    :param column_names: Dictionary with the column names for the isochrones,
        defaults to None
    :type column_names: dict | None
    :param model: Isochrone model name.
    :type model: str
    :param magnitude: Magnitude filter to load.
    :type magnitude: str
    :param color: Tuple with the two filters to generate the first color.
    :type color: tuple
    :param color2: Tuple with the two filters to generate the second color,
        defaults to None
    :type color2: tuple | None

    :return: List of filters plus EEP, mass, metallicity, and age column names.
    :rtype: tuple[list, str | None, str, str, str]
    """
    if column_names is None:
        eep_col = phot_systs_data[model]["col_names"]["eep_col"]
        mass_col = phot_systs_data[model]["col_names"]["mass_col"]
        met_col = phot_systs_data[model]["col_names"]["met_col"]
        age_col = phot_systs_data[model]["col_names"]["age_col"]
    else:
        eep_col = column_names["eep_col"]
        mass_col = column_names["mass_col"]
        met_col = column_names["met_col"]
        age_col = column_names["age_col"]

    # Select columns to keep
    all_filters = [magnitude] + list(color)
    if color2 is not None:
        all_filters += list(color2)

    return all_filters, eep_col, mass_col, met_col, age_col


def extract_paths(isochs_path: str) -> list:
    """Extract isochrone files from `isochs_path`.

    :param isochs_path: Path to the isochrone files.
    :type isochs_path: str

    :raises FileNotFoundError: If no files are found in the isochrones path.

    :return: List of isochrone file paths.
    :rtype: list
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
    parsec_rm_stage_9: bool,
    f_paths: list,
    all_filters: list,
    eep_col: str | None,
    mass_col: str,
    met_col: str,
    age_col: str,
) -> tuple[list[str], list[np.ndarray]]:
    """Read isochrone files and store them as numpy arrays along with its associated
    metallicity and age values.

    :param model: Isochrone model name.
    :type model: str
    :param parsec_rm_stage_9: Remove post-AGB stage for PARSEC models.
    :type parsec_rm_stage_9: bool
    :param f_paths: List of isochrone file paths.
    :type f_paths: list
    :param all_filters: List of photometric magnitude  + color(s)
    :type all_filters: list
    :param eep_col: Column name for EEP
    :type eep_col: str | None
    :param mass_col: Column name for initial mass
    :type mass_col: str
    :param met_col: Metallicity column name.
    :type met_col: str
    :param age_col: Age column name.
    :type age_col: str

    :return: First list contains the met and age values (as strings), second list
     contains the associated isochrones as numpy arrays
    :rtype: tuple[list[str], list[np.ndarray]]
    """

    # MIST isochrones require 'eep_col'
    cols_keep = all_filters + [mass_col]
    if eep_col is not None:
        cols_keep += [eep_col]

    met_age_vals, isoch_arrays = [], []
    for file_path in f_paths:
        # Extract columns names and full header
        col_names, full_header = get_header(model, file_path)
        # Columns to keep for this photometric system
        cols_keep_ps = list(set(col_names) & set(cols_keep))

        # Load file
        data = np.genfromtxt(
            file_path,
            comments=phot_systs_data[model]["comment_char"],
            names=col_names,
            # delimiter=phot_systs_data[model]["sep_cols"],
            encoding="utf-8",
            deletechars="",
            autostrip=True,
        )

        if model == "PARSEC":
            # Get unique metallicity values
            unique_mets = np.unique(data[met_col])
            for met_val in unique_mets:
                # Filter by metallicity
                met_mask = data[met_col] == met_val
                met_data = data[met_mask]

                # Get unique age values within this metallicity
                unique_ages = np.unique(met_data[age_col])
                for age_val in unique_ages:
                    # Filter by age
                    age_mask = met_data[age_col] == age_val
                    df = met_data[age_mask]

                    # Remove post-AGB stage
                    if parsec_rm_stage_9 is True:
                        stage_col = phot_systs_data[model]["parsec_stage_9_col"]
                        stage_id = phot_systs_data[model]["parsec_stage_9_id"]
                        msk = df[stage_col] != stage_id
                        df = df[msk]

                    # Extract met, age values
                    met = str(df[met_col][0])
                    age = str(df[age_col][0])

                    # Store data
                    met_age_vals.append([met, age])
                    isoch_arrays.append(df[cols_keep_ps])

        elif model == "MIST":
            met = get_MIST_z_val(met_col, full_header)

            # Get unique age values
            unique_ages = np.unique(data[age_col])
            for age_val in unique_ages:
                # Filter by age
                age_mask = data[age_col] == age_val
                df = data[age_mask]

                age = str(df[age_col][0])

                # Store data
                met_age_vals.append([met, age])
                isoch_arrays.append(df[cols_keep_ps])

        elif model == "BASTI":
            met, age = get_BASTI_z_a_val(full_header, met_col, age_col)

            # Store data
            met_age_vals.append([met, age])
            isoch_arrays.append(data[cols_keep_ps])

    return met_age_vals, isoch_arrays


def get_header(model: str, file_path: str) -> tuple[list, list]:
    """Iterate through each line in the file to get the header.
    Extract the column names from the header.

    :param model: Name of the isochrones model used
    :type model: str
    :param file_path: Path to the isochrone file.
    :type file_path: str

    :return: List of column names and full header.
    :rtype: tuple[list, list]
    """
    # Extract full header
    with open(file_path, mode="r") as f_iso:
        full_header = []
        for line in f_iso:
            if not line.startswith(phot_systs_data[model]["comment_char"]):
                break
            full_header.append(line)

    # Extract column names
    column_line = full_header[phot_systs_data[model]["idx_col_line"]]
    column_names = column_line.replace(
        phot_systs_data[model]["comment_char"], ""
    ).split()

    return column_names, full_header


def get_MIST_z_val(met_col: str, full_header: list) -> str:
    """Extract metallicity value for MIST files.

    :param met_col: Metallicity column name.
    :type met_col: str
    :param full_header: Full header of the isochrone file.
    :type full_header: list

    :raises ValueError: If the metallicity cannot be read from the header.

    :return: Metallicity value.
    :rtype: str
    """
    met = None
    for i, line in enumerate(full_header):
        line = line.replace(phot_systs_data["MIST"]["comment_char"], "")
        if met_col in line:
            # Identify lines that contain the metallicity column
            z_header_cols = line.split()
            # Assume the metallicity value is in the next line
            next_line = full_header[i + 1].replace(
                phot_systs_data["MIST"]["comment_char"], ""
            )
            z_header_vals = next_line.split()
            # Extract metallicity value as string
            met_idx = z_header_cols.index(met_col)
            met = z_header_vals[met_idx]
            break

    if met is None:
        raise ValueError("Could not read header from MIST isochrone")

    return met


def get_BASTI_z_a_val(full_header: list, met_col: str, age_col: str) -> tuple[str, str]:
    """Extract metallicity and age values for BASTI files.

    :param full_header: Full header of the isochrone file.
    :type full_header: list
    :param met_col: Metallicity column name.
    :type met_col: str
    :param age_col: Age column name.
    :type age_col: str

    :raises ValueError: If the metallicity or age cannot be read from the header.

    :return: Metallicity and age values.
    :rtype: tuple[str, str]
    """
    met, age = None, None
    for line in full_header:
        line = line.replace(phot_systs_data["BASTI"]["comment_char"], "").strip()
        if met_col in line:
            # Assume both values exist in the same line
            met = line.split(met_col)[1].split()[0]
            age = line.split(age_col)[1].split()[0]

            # Basti ages are expressed in Myr, convert to log10
            if phot_systs_data["BASTI"]["myr_to_log10"]:
                age = np.log10(float(age) * 1e6)
                # Back to (rounded) string
                age = str(round(age, 5))
            break

    if met is None or age is None:
        raise ValueError("Could not read header from Basti isochrone")

    return met, age


def sort_by_met_age(
    met_age_vals: list,
    isoch_arrays: list[np.ndarray],
) -> tuple[list, list]:
    """
    Sort the isochrones by metallicity and age.

    :param met_age_vals: List of metallicity and ages
    :type met_age_vals: list
    :param isoch_arrays: List of isochrones as numpy arrays
    :type isoch_arrays: list[np.ndarray]

    :return: Lists with sorted met and age values, and sorted isochrones.
    :rtype: tuple[list, list]
    """
    # Convert the string elements to floats for sorting
    met_age_f = [(float(x), float(y)) for x, y in met_age_vals]

    # Get the indexes that sort the list by metallicities first and ages second
    sorted_indexes = sorted(
        range(len(met_age_f)), key=lambda i: (met_age_f[i][0], met_age_f[i][1])
    )
    # Sort metallicity and age values
    met_age_sorted = [met_age_vals[i] for i in sorted_indexes]
    # Sort data frames
    isochs_sorted = [isoch_arrays[i] for i in sorted_indexes]

    return met_age_sorted, isochs_sorted


def interp_isochrones(
    N_interp: int,
    mass_col: str,
    eep_col: str | None,
    met_age_sorted: list,
    isochs_sorted: list[np.ndarray],
) -> dict:
    """
    Interpolate the isochrone data based on either mass or EEP distribution.

    :param N_interp: Number of points to interpolate (for mass distribution).
    :type N_interp: int
    :param mass_col: Column name for the initial mass (for mass distribution).
    :type mass_col: str
    :param eep_col: Column name for the EEP (for EEP distribution).
    :type eep_col: str | None
    :param met_age_sorted: List of sorted metallicity and ages.
    :type met_age_sorted: list
    :param isochs_sorted: List of sorted isochrones as numpy arrays.
    :type isochs_sorted: list[np.ndarray]

    :return: Dictionary of interpolated isochrones.
    :rtype: dict
    """
    if eep_col is None:
        # Define mass distribution parameters
        percs = np.array([0.5, 1, 1.5, 2, 5]) / np.sum([0.5, 1, 1.5, 2, 5])
        b1, b2, b3, b4 = [int(N_interp * p) for p in percs[:4]]
        b5 = N_interp - (b1 + b2 + b3 + b4)
        # Per isochrone mass based distribution
        mass_percs = np.array([0.25, 0.5, 0.75, 0.9])

        # General distribution, isochrone independent. The above works a bit better
        # mp = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 1.])
        # m1 = np.linspace(mp[0], mp[1], b1, endpoint=False)
        # m2 = np.linspace(mp[1], mp[2], b2, endpoint=False)
        # m3 = np.linspace(mp[2], mp[3], b3, endpoint=False)
        # m4 = np.linspace(mp[3], mp[4], b4, endpoint=False)
        # m5 = np.linspace(mp[4], mp[5], b5)
        # mass_dist = np.concatenate([m1, m2, m3, m4, m5])
    else:
        # Define EEP distribution
        eep_min = int(min([min(isoch[eep_col]) for isoch in isochs_sorted]))
        eep_max = int(max([max(isoch[eep_col]) for isoch in isochs_sorted]))
        eep_dist = np.arange(eep_min, eep_max + 1)

    isochrones = {}
    for i, (met, age) in enumerate(met_age_sorted):
        # # Generate entry in dictionary
        # try:
        #     isochrones[met]
        # except KeyError:
        #     isochrones[met] = {}
        # try:
        #     isochrones[met][age]
        # except KeyError:
        #     isochrones[met][age] = []

        # Create nested dictionary structure (equivalent to the block above)
        isochrones.setdefault(met, {}).setdefault(age, [])

        isoch = isochs_sorted[i]

        if eep_col is None:
            # Mass-based interpolation
            mmin, mmax = isoch[mass_col].min(), isoch[mass_col].max()
            mstep1, mstep2, mstep3, mstep4 = mmax * mass_percs  # pyright: ignore
            m1 = np.linspace(mmin, mstep1, b1, endpoint=False)  # pyright: ignore
            m2 = np.linspace(mstep1, mstep2, b2, endpoint=False)  # pyright: ignore
            m3 = np.linspace(mstep2, mstep3, b3, endpoint=False)  # pyright: ignore
            m4 = np.linspace(mstep3, mstep4, b4, endpoint=False)  # pyright: ignore
            m5 = np.linspace(mstep4, mmax, b5)  # pyright: ignore
            mass_dist = np.concatenate([m1, m2, m3, m4, m5])
            interp_dist, interp_col = mass_dist, isoch[mass_col]
            # # General distribution, isochrone independent
            # interp_dist, interp_col = mass_dist, np.linspace(0.0, 1.0, len(isoch))
            # Keep all columns
            cols_keep = isoch.dtype.names
        else:
            # EEP-based interpolation
            interp_dist, interp_col = eep_dist, isoch[eep_col]  # pyright: ignore
            # Remove EEP column
            cols_keep = [_ for _ in isoch.dtype.names if _ != eep_col]

        # Interpolate all columns
        isoch_interp = {}
        for col in cols_keep:
            isoch_interp[col] = np.interp(interp_dist, interp_col, isoch[col])
        isochrones[met][age].append(isoch_interp)

    return isochrones


def merge_ps_massini_check(mass_col: str, isochrones: dict) -> tuple[list, list]:
    """Combine photometric systems, and check initial masses.

    1. Combine photometric systems if more than one was used.
    2. Check that initial masses are equal across all systems

    The isochrone parameter 'mass_col' is assumed to be equal across
    photometric systems, for a given metallicity and age. We check here that
    this is the case.

    :param mass_col: Name of the mass column.
    :type mass_col: str
    :param isochrones: Dictionary of isochrones.
    :type isochrones: dict

    :raises ValueError: If initial mass values differ across photometric systems or
        if there is a shape mismatch in the loaded isochrones

    :return: Ordered isochrones and metallicity-age array.
    :rtype: tuple[list, list]
    """
    met_age_arr, isochrones_ordered = [], []
    for met_k in isochrones.keys():
        met_vals = []
        for age_k, age_list in isochrones[met_k].items():
            met_age_arr.append([met_k, age_k])

            df_age_0 = age_list[0]
            # Check masses if more than one photometric system was used
            if len(age_list) > 1:
                for df_age_X in age_list[1:]:
                    mass_diff = abs(df_age_0[mass_col] - df_age_X[mass_col]).sum()
                    if mass_diff > 0.001:
                        raise ValueError(
                            "Initial mass values differ across photometric systems"
                        )
                    # Drop 'mass_col' column from this photometric system
                    del df_age_X[mass_col]
                    # Merge data frames and replace original one
                    df_age_0.update(df_age_X)

            met_vals.append(df_age_0)
        isochrones_ordered.append(met_vals)

    # Check array sizes and mass order (min to max)
    Nvals = []
    for met in isochrones_ordered:
        for age in met:
            for k, arr in age.items():
                Nvals.append(arr.shape[0])
                if k == mass_col:
                    # Use a mass tolerance instead of the 'sort_flag':
                    #
                    # sort_flag = all(arr[i] <= arr[i + 1] for i in range(len(arr) - 1))
                    #
                    # which checks if all masses are perfectly sorted, because the
                    # PARSEC isochrones contain some rounding errors of the order
                    # 10**-7, 10**-8 for large masses, which would always set the flag
                    # to False
                    max_mass_diff = max(
                        (arr[i] - arr[i + 1]) for i in range(len(arr) - 1)
                    )
                    if max_mass_diff > 0.01:
                        raise ValueError(
                            "Masses are not min-max sorted in loaded isochrones."
                        )

    if not np.all(Nvals):
        raise ValueError(
            "Shape mismatch in loaded isochrones. This usually means that an\n"
            + "incorrect number of ages and/or metallicities are stored in the\n"
            + "isochrone files."
        )

    return isochrones_ordered, met_age_arr


def reshape_isochrones(
    magnitude: str,
    color: tuple,
    color2: tuple | None,
    mass_col: str,
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

    # These columns are added later by the synthetic class:
    fb    : binary magnitude --
    cXb   : binary colors     |
    Minib : binary masses -----

    theor_tracks.shape = (Nz, Na, Nd, Ni)
    Nz: number of metallicities
    Na: number of log(age)s
    Nd: number of data columns
    Ni: number of interpolated values

    :param magnitude: Magnitude filter to load.
    :type magnitude: str
    :param color: Tuple with the two filters to generate the first color.
    :type color: tuple
    :param color2: Tuple with the two filters to generate the second color,
        defaults to None
    :type color2: tuple | None
    :param mass_col: Name of the initial mass column.
    :type mass_col: str
    :param isochrones: List of isochrones.
    :type isochrones: list

    :return: Reshaped isochrones array and list of color filters.
    :rtype: tuple[np.ndarray, list]
    """

    # Notice that Ni and Nd are inverted here compared to the shape of the final
    # 'theor_tracks' array
    Nz = len(isochrones)
    Na = len(isochrones[0])
    Ni = isochrones[0][0][magnitude].shape[0]

    all_colors = [color]
    if color2 is not None:
        all_colors.append(color2)
    N_colors = len(all_colors)

    # Array that will store all the interpolated tracks. The '2' makes room for the
    # magnitude and first color which are required. Adding 'N_colors' to that value
    # makes room for an optional second color, and the initial mass column.
    theor_tracks = np.zeros([Nz, Na, (2 + N_colors), Ni])

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
                # Store color. The '1' makes room for the magnitude that goes first.
                theor_tracks[i][j][1 + k] = f1 - f2

                # Individual filters for colors, used for binarity
                cols_dict[color_filts[0]] = f1
                cols_dict[color_filts[1]] = f2

            met_lst.append(cols_dict)

            # Add initial mass column to the end of the array
            theor_tracks[i][j][2 + N_colors - 1] = df_age[mass_col]
        color_filters.append(met_lst)

    return theor_tracks, color_filters
