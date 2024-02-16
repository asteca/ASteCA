import numpy as np
import pandas as pd
from typing import Tuple, List


def get(model, N_interp, f_paths, met_col, age_col, cols_keep) -> Tuple[List, List]:
    """ """
    if model == "PARSEC":
        isochrones, met_age_arr = read_PARSEC_files(
            f_paths, met_col, age_col, cols_keep, N_interp
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


def read_PARSEC_files(f_paths, met_col, age_col, cols_keep, N_interp):
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
