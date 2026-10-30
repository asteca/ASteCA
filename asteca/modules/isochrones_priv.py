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
        "parsec_stages_col": "label",  # column name for stages
        "max_mag": 50.0,  # The new files (July 2026) contain bad '999.999' values
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
    parsec_rm_stages: float,
    mass_interpolation: str,
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
    :param parsec_rm_stages: All stages with values equal or larger than this value
        will be removed from the PARSEC models
    :type parsec_rm_stages: float

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
        parsec_rm_stages,
        f_paths,
        all_filters,
        eep_col,
        mass_col,
        met_col,
        age_col,
    )

    met_age_sorted, isochs_sorted = sort_by_met_age(met_age_vals, isoch_arrays)

    # Interpolate the isochrones based on either mass or EEP distribution
    isochrones = interp_isochrones(
        N_interp,
        mass_col,
        eep_col,
        met_age_sorted,
        isochs_sorted,
        mass_interpolation,
        all_filters,
    )

    isochrones, met_age_arr = merge_ps_massini_check(mass_col, isochrones)

    # Create dictionary of mets and ages
    all_m, all_a = np.array(met_age_arr).T
    # Discard duplicate values
    all_m = np.array(list(dict.fromkeys(all_m))).astype(float)
    all_a = np.array(list(dict.fromkeys(all_a))).astype(float)
    met_age_dict = {"met": all_m, "loga": all_a}

    if (
        len(met_age_dict["loga"]) > 1
        and max(met_age_dict["loga"][1:] - met_age_dict["loga"][:-1]) > 0.25
    ):
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
    parsec_rm_stages: float,
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
    :param parsec_rm_stages: Remove stages from PARSEC models.
    :type parsec_rm_stages: float
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

                    # Remove bad values in filters
                    mag_msk = np.logical_and.reduce(
                        [
                            df[_filter] < phot_systs_data[model]["max_mag"]
                            for _filter in all_filters
                        ]
                    )
                    df = df[mag_msk]

                    # Remove all stages equal or larger than this value
                    df = df[
                        df[phot_systs_data[model]["parsec_stages_col"]]
                        < parsec_rm_stages
                    ]

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
    mass_interpolation: str,
    all_filters: list,
    mass_tol: float = 1e-6,
) -> dict:
    """
    Interpolate the isochrone data based on either mass or EEP distribution.

    The variables b1, ..., b5 determine how many of the N_interp interpolation
    points are assigned to five different mass intervals.

    The intention is to sample the high-mass end much more densely, because the
    upper isochrone changes more rapidly and contains crucial evolutionary phases
    where finer interpolation is useful.

    :param N_interp: Number of points to interpolate (for mass distribution).
    :type N_interp: int
    :param mass_col: Column name for the initial mass.
    :type mass_col: str
    :param eep_col: Column name for the EEP.
    :type eep_col: str | None
    :param met_age_sorted: List of sorted metallicity and ages.
    :type met_age_sorted: list
    :param isochs_sorted: List of sorted isochrones as numpy arrays.
    :type isochs_sorted: list[np.ndarray]
    :param mass_interpolation: One of 'mass_uniform', 'mass_skewed', or
        'curve_uniform'. The latter parameterises points by normalised cumulative
        arc-length distance along the (color_col, mag_col) curve instead of mass,
        which better preserves shape in regions of high CMD curvature.
    :type mass_interpolation: str
    :param color_col: Column name for the color, required if
        mass_interpolation == 'curve_uniform'.
    :type color_col: str | None
    :param mag_col: Column name for the magnitude, required if
        mass_interpolation == 'curve_uniform'.
    :type mag_col: str | None

    :raises ValueError: If the isochrone dtype does not contain column names.

    :return: Dictionary of interpolated isochrones.
    :rtype: dict
    """
    # Dummy value to avoid Pyright "possibly unbound" error
    eep_dist = []

    # Mass-based distribution parameters
    mass_fracs = np.array([0.25, 0.5, 0.75, 0.9], dtype=float)
    weights = np.array([0.5, 1, 1.5, 2, 5], dtype=float)
    percs = weights / weights.sum()
    b1, b2, b3, b4 = [int(N_interp * p) for p in percs[:4]]
    b5 = N_interp - (b1 + b2 + b3 + b4)

    if eep_col is not None:
        # EEP-based global distribution
        if not isochs_sorted:
            raise ValueError("Empty 'isochs_sorted'")

        # Determine the global min and max EEP values across all isochrones
        eep_min = int(min(np.min(isoch[eep_col]) for isoch in isochs_sorted))
        eep_max = int(max(np.max(isoch[eep_col]) for isoch in isochs_sorted))
        eep_dist = np.arange(eep_min, eep_max + 1)

    isochrones = {}
    for (met, age), isoch in zip(met_age_sorted, isochs_sorted):
        if isoch.dtype.names is None:
            raise ValueError(
                f"No 'names' found for isoch.dtype: (met, age)=({met}, {age})"
            )

        # Create nested dictionary structure for metallicity and age
        isochrones.setdefault(met, {}).setdefault(age, [])

        if mass_interpolation == "adaptive":
            # ---- Curve-arclength-based interpolation, uniform ----
            isoch_interp = adaptive_interpolation(
                age, isoch, all_filters, mass_col, N_interp
            )

        else:
            if eep_col is None:
                # interp_col = isoch[mass_col]
                mass = isoch[mass_col]
                if np.any(np.diff(mass) < -mass_tol):
                    raise ValueError(
                        f"Masses not monotonic in raw isochrone: (met, age)=({met}, {age})"
                    )
                # Numerical noise below tol: clamp without reordering
                interp_col = np.maximum.accumulate(mass)

                if mass_interpolation == "mass_uniform":
                    # ---- Mass-based interpolation, uniform ----
                    mmin = np.min(isoch[mass_col])
                    mmax = np.max(isoch[mass_col])
                    interp_dist = np.linspace(mmin, mmax, N_interp)

                elif mass_interpolation == "mass_skewed":
                    # ---- Mass-based interpolation, skewed ----
                    mmin = np.min(isoch[mass_col])
                    mmax = np.max(isoch[mass_col])
                    mstep1, mstep2, mstep3, mstep4 = mmax * mass_fracs
                    m1 = np.linspace(mmin, mstep1, b1, endpoint=False)
                    m2 = np.linspace(mstep1, mstep2, b2, endpoint=False)
                    m3 = np.linspace(mstep2, mstep3, b3, endpoint=False)
                    m4 = np.linspace(mstep3, mstep4, b4, endpoint=False)
                    m5 = np.linspace(mstep4, mmax, b5)
                    interp_dist = np.concatenate((m1, m2, m3, m4, m5))
                else:
                    raise ValueError(
                        f"Invalid mass_interpolation value: {mass_interpolation}"
                    )

                cols_keep = isoch.dtype.names

            else:
                # ---- EEP-based interpolation ----
                interp_dist = eep_dist
                interp_col = isoch[eep_col]
                cols_keep = tuple(c for c in isoch.dtype.names if c != eep_col)

            # ---- Interpolate columns ----
            isoch_interp: dict = {
                col: np.interp(interp_dist, interp_col, isoch[col]) for col in cols_keep
            }

        # Store the interpolated isochrone in the dictionary
        isochrones[met][age].append(isoch_interp)

    return isochrones


def adaptive_interpolation(
    age: float, isoch: np.ndarray, all_filters, mass_col, N_interp, min_dist=0.001
) -> dict:
    """
    Interpolate the isochrone data based on normalized cumulative arc-length distance
    """
    filter_1, filter_2 = isoch[all_filters[1]], isoch[all_filters[2]]
    color = filter_1 - filter_2
    mag = isoch[all_filters[0]]
    mass_ini = isoch[mass_col]

    # Normalized color, magnitude in [0, 1]
    color_norm = (color - np.min(color)) / (np.max(color) - np.min(color))
    mag_norm = (mag - np.min(mag)) / (np.max(mag) - np.min(mag))

    # Iterate through each point in the CMD
    prev_cmd_point = np.array([color[0], mag[0]])
    prev_cmd_point_norm = np.array([color_norm[0], mag_norm[0]])
    prev_filter_1, prev_filter_2 = filter_1[0], filter_2[0]
    prev_mass_ini = mass_ini[0]

    mag_interp, filter_1_interp, filter_2_interp, mass_ini_interp = [], [], [], []
    for point in zip(
        color_norm[1:],
        mag_norm[1:],
        color[1:],
        mag[1:],
        filter_1[1:],
        filter_2[1:],
        mass_ini[1:],
    ):
        (
            next_color_norm,
            next_mag_norm,
            next_color,
            next_mag,
            next_f1,
            next_f2,
            next_mini,
        ) = point

        next_cmd_point = np.array([next_color, next_mag])
        next_cmd_point_norm = np.array([next_color_norm, next_mag_norm])

        dist = np.linalg.norm(next_cmd_point_norm - prev_cmd_point_norm)

        if dist <= min_dist:
            mag_interp.append(next_mag)
            filter_1_interp.append(next_f1)
            filter_2_interp.append(next_f2)
            mass_ini_interp.append(next_mini)
            prev_cmd_point = next_cmd_point
            prev_cmd_point_norm = next_cmd_point_norm
            prev_filter_1, prev_filter_2 = next_f1, next_f2
            prev_mass_ini = next_mini
            continue

        N_interp_pts = int(np.ceil(dist / min_dist))
        if N_interp_pts < 2:
            mag_interp.append(next_mag)
            filter_1_interp.append(next_f1)
            filter_2_interp.append(next_f2)
            mass_ini_interp.append(next_mini)
            prev_cmd_point = next_cmd_point
            prev_cmd_point_norm = next_cmd_point_norm
            prev_filter_1, prev_filter_2 = next_f1, next_f2
            prev_mass_ini = next_mini
            continue

        interp_points = np.linspace(prev_cmd_point, next_cmd_point, N_interp_pts)
        mag_interp.extend(interp_points[:, 1])

        color_interp = interp_points[:, 0]
        f1_interp = np.linspace(prev_filter_1, next_f1, N_interp_pts)
        f2_interp = np.linspace(prev_filter_2, next_f2, N_interp_pts)

        filter_1_interp.extend(color_interp + f2_interp)
        filter_2_interp.extend(f1_interp - color_interp)

        mass_ini_interp.extend(np.linspace(prev_mass_ini, next_mini, N_interp_pts))

        prev_cmd_point = next_cmd_point
        prev_cmd_point_norm = next_cmd_point_norm
        prev_filter_1, prev_filter_2 = next_f1, next_f2
        prev_mass_ini = next_mini

    # Remove points non-uniformly from interpolated arrays until the number of
    # total points == N_interp, preferentially keeping points around the
    # index closest to mass=1 Msun (typical MS+turn-off region)
    mag_interp = np.array(mag_interp)
    filter_1_interp = np.array(filter_1_interp)
    filter_2_interp = np.array(filter_2_interp)
    mass_ini_interp = np.array(mass_ini_interp)
    N_old, m_ini_min, m_ini_max = (
        len(mag_interp),
        mass_ini_interp[0],
        mass_ini_interp[-1],
    )
    mass_mid_point = 1
    if N_old > N_interp:
        # mass_ini_interp is monotonic increasing; find the first index
        # where mass crosses 1 Msun, starting from the beginning
        mass_mid_point = min(mass_mid_point, (m_ini_min + m_ini_max) / 2)
        i_ref = np.searchsorted(mass_ini_interp, mass_mid_point)
        i_ref = np.clip(i_ref, 0, N_old - 1)

        # Weight by distance in index space from i_ref (peaked there)
        idx_arr = np.arange(N_old)
        scale = N_old / 10  # controls how sharply weight falls off
        weights = np.exp(-0.5 * ((idx_arr - i_ref) / scale) ** 2)
        weights += 1e-3  # floor so no region has exactly zero probability
        weights /= weights.sum()

        idx = np.random.choice(N_old, size=N_interp - 2, replace=False, p=weights)
        idx = np.concatenate(([0], idx, [N_old - 1]))  # force-keep endpoints
        idx = np.unique(idx)
        if len(idx) < N_interp:
            remaining = np.setdiff1d(idx_arr, idx)
            extra = np.random.choice(
                remaining,
                size=N_interp - len(idx),
                replace=False,
                p=weights[remaining] / weights[remaining].sum(),
            )
            idx = np.concatenate((idx, extra))
        idx = np.sort(idx)

        mag_interp = mag_interp[idx]
        filter_1_interp = filter_1_interp[idx]
        filter_2_interp = filter_2_interp[idx]
        mass_ini_interp = mass_ini_interp[idx]

    print(
        f"logage={age}, mass_mid_point:{mass_mid_point:.3f}, N: {len(mag)} --> {N_old} --> {len(mag_interp)} ; Mini: [{m_ini_min:.3f}, {m_ini_max:.2f}] --> [{mass_ini_interp[0]:.3f}, {mass_ini_interp[-1]:.2f}]"
    )

    isoch_interp = {
        all_filters[0]: mag_interp,
        all_filters[1]: filter_1_interp,
        all_filters[2]: filter_2_interp,
        mass_col: mass_ini_interp,
    }

    # if mass_ini_interp[-1]-mass_ini_interp[0]<0.2:
    #     import matplotlib.pyplot as plt

    #     plt.subplot(121)
    #     plt.scatter(np.array(filter_1_interp) - np.array(filter_2_interp), mag_interp, s=5, alpha=.5)
    #     plt.scatter(isoch[all_filters[1]] - isoch[all_filters[2]], isoch[all_filters[0]], c='k', s=10)
    #     plt.gca().invert_yaxis()
    #     plt.subplot(122)
    #     plt.hist(mass_ini_interp)
    #     plt.show()
    #     breakpoint()

    return isoch_interp


def merge_ps_massini_check(mass_col: str, isochrones: dict) -> tuple[list, list]:
    """Combine photometric systems, and check initial masses.

    1. Combine photometric systems if more than one was used.
    2. Check that initial masses are equal across all systems
    3. Check that all isochrones have the same number of values

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
    # Process all isochrones, merging photometric systems if needed
    met_age_arr, isochrones_ordered = [], []
    for met_k, age_dict in isochrones.items():
        met_vals = []
        for age_k, age_list in age_dict.items():
            met_age_arr.append([met_k, age_k])

            df_age_0 = age_list[0]
            # Check masses if more than one photometric system was used
            if len(age_list) > 1:
                for df_age_X in age_list[1:]:
                    mass_diff = abs(df_age_0[mass_col] - df_age_X[mass_col]).sum()
                    if mass_diff > 0.001:
                        raise ValueError(
                            "\nMore than one photometric system was detected at "
                            f"(met={met_k}, age{age_k}) and the initial mass values "
                            "differ across them."
                        )
                    # Drop 'mass_col' column from this photometric system
                    del df_age_X[mass_col]
                    # Merge data frames and replace original one
                    df_age_0.update(df_age_X)

            met_vals.append(df_age_0)
        isochrones_ordered.append(met_vals)

    # Check whether the (metallicity, age) grid is regular
    met_age = np.asarray(met_age_arr, dtype=float)
    met_vals = np.unique(met_age[:, 0])
    age_vals = np.unique(met_age[:, 1])

    issues = []

    def _find_irregular_spacing(values: np.ndarray, name: str) -> list[str]:
        if len(values) <= 2:
            return []

        diffs = np.diff(values)
        step = np.median(diffs)
        bad_idx = np.where(~np.isclose(diffs, step))[0]

        if not len(bad_idx):
            return []

        lines = [f"Irregular {name} spacing (expected step ~{step:g}):"]
        for i in bad_idx:
            i0 = max(0, i - 1)
            i1 = min(len(values), i + 3)
            context = values[i0:i1]
            lines.append(
                f"  around {values[i]:g} -> {values[i + 1]:g} "
                f"(step={diffs[i]:g}): {context}"
            )

        return lines

    issues.extend(_find_irregular_spacing(met_vals, "metallicity"))
    issues.extend(_find_irregular_spacing(age_vals, "age"))

    # Check for missing (metallicity, age) combinations
    grid = {(met, age) for met, age in met_age}
    for met in met_vals:
        ages_present = np.sort(np.array([age for m, age in grid if np.isclose(m, met)]))
        missing = [age for age in age_vals if (met, age) not in grid]
        for age in missing:
            idx = np.searchsorted(ages_present, age)
            adjacent = ages_present[max(0, idx - 1) : min(len(ages_present), idx + 1)]
            issues.append(
                f"Missing grid point (met={met:g}, age={age:g}); "
                f"adjacent ages at this metallicity: {adjacent}"
            )

    if issues:
        warnings.warn(
            "The (metallicity, age) isochrone grid is irregular:\n" + "\n".join(issues)
        )

    # Check array sizes and mass order (min to max)
    Nvals = []
    for met in isochrones_ordered:
        for age in met:
            for arr in age.values():
                Nvals.append(arr.shape[0])

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
    ageX = [f, c1, c2,.., Mini, fb, c1b, c2b,.., Minib]

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
