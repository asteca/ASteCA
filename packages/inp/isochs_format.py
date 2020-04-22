
import sys


def age_format(evol_track):
    """
    Define reg expression to isolate the age of an isochrone from its
    commented title line.
    """
    if evol_track[:3] == 'PAR':
        age_format = r"Age = (.+?) yr"
    else:
        # TODO in place for #275
        pass

    return age_format


def line_start_format(evol_track):
    """
    Return the format of the line where column names are read from, for each
    set of evolutionary tracks.
    """
    if evol_track[:3] == 'PAR':
        # String that identifies the beginning of a new isochrone.
        line_start = "# Zini"
    else:
        # TODO in place for #275
        pass

    return line_start


def read_line_start(met_f, evol_track, line_start):
    """
    Read the line where the filter names are positioned, for each set of
    evolutionary tracks.
    """
    if evol_track[:3] == 'PAR':
        # Open the metallicity file.
        with open(met_f, mode="r") as f_iso:
            # Iterate through each line in the file.
            for line in f_iso:
                # When this line is found, extract columns names from the
                # following line.
                if line.startswith(line_start):
                    # Split into column names.
                    ls = line.split()
                    # Remove comment character.
                    if ls[0] == '#':
                        del ls[0]
                    else:
                        sys.exit("\nERROR: theoretical isochrones are not\n"
                                 "formatted as expected!\n")
                    break
    else:
        # TODO in place for #275
        pass

    return ls


def common_ids(evol_track, col_names, line):
    """
    Find the indexes of the filter / extra pars to read in each metallicity
    file.

    # TODO #275: does this apply to other tracks?
    """
    if evol_track[:3] == 'PAR':
        col_ids = []
        for col in col_names:
            col_ids.append(line.index(col))

        return col_ids

    else:
        # TODO in place for #275
        pass
