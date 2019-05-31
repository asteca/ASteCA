
import os


def main():
    '''
    Dictionary that stores the names and column names for each
    filter defined in each photometric system, as presented in the CMD
    Girardi et al. service: http://stev.oapd.inaf.it/cgi-bin/cmd

    Data stored in the 'CMD_systs.dat' file.

    Capitalization of the filter names matters!

    cmd_systs = {0, 1, ..., N}
    X = ('nameX', filtersX, effWavelengthsX)
    filtersX = (filter1, filter2, ..., filterQ)
    effWavelengthsX = (float1, float2, ..., floatQ)

    Returns
    -------
    cmd_systs : dictionary
    '''
    fn = os.path.join(os.path.dirname(__file__), 'CMD_systs.dat')
    cmd_systs = {}
    with open(fn) as f:
        f_lines = f.readlines()
        for i, li in enumerate(f_lines):
            if not li.startswith("#"):
                # Photometric system's name, and filters' names.
                ls = li.split()
                Nf = int((len(ls) - 1) / 3.)
                cmd_systs[str(ls[0])] = (
                    ls[1], tuple(ls[2:Nf + 2]),
                    tuple(map(float, ls[Nf + 2:2 * Nf + 2])))

    return cmd_systs


if __name__ == '__main__':
    main()
