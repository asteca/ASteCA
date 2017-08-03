
import os


def effWavel(line):
    """
    Properly read effective wavelengths from line.
    """
    effW = []
    for l0 in line.split():
        effW.append(float(l0.split("(")[0]))

    return tuple(effW)


def main():
    '''
    Dictionary that stores the names and column names for each
    filter defined in each photometric system, as presented in the CMD
    Girardi et al. service: http://stev.oapd.inaf.it/cgi-bin/cmd

    Data stored in the 'CMD_systs.dat' file.

    Capitalization of the filter names matters!

    cmd_systs = {photSyst0, photSyst1, ..., photSystN}
    photSyst['X'] = ('nameX', filtersX, effWavelengthsX)
    filtersX = (filter1, filter2, ..., filterQ)
    effWavelengthsX = (float1, float2, ..., floatQ)

    Returns
    -------
    cmd_systs : dictionary

    '''
    fn = os.path.join(os.path.dirname(__file__), 'CMD_systs.dat')
    cmd_systs = {}
    with open(fn) as f:
        j = 0
        f_lines = f.readlines()
        for i, li in enumerate(f_lines):
            if not li.startswith("#"):
                # Read odd lines
                if i % 2 != 0:
                    # Photometric system's name, and filters' names.
                    ls = li.split()
                    # Filters' effective wavelengths.
                    effW = effWavel(f_lines[i + 1])
                    cmd_systs[str(j)] = (ls[0], tuple(ls[1:]), effW)
                    j += 1

    return cmd_systs


if __name__ == '__main__':
    main()
