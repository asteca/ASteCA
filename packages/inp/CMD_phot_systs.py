
from pathlib import Path


def main():
    """
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
    """
    fn = Path.cwd() / 'packages' / 'defvals' / 'CMD_systs.dat'

    cmd_systs = {}
    with open(fn) as f:
        f_lines = f.readlines()
        for i, li in enumerate(f_lines):
            if not li.startswith("#"):
                # Photometric system's name, and filters' names.
                _id = li.split()[0]
                ls = li[105:].split()
                phot_syst = ls[0]
                # Number of filters/lambdas/omegas
                Nf = int((len(ls) - 1) / 3.)
                # Filters
                filters = tuple(ls[1:Nf + 1])
                # Lambdas
                lambdas = tuple(map(float, ls[Nf + 1:2 * Nf + 1]))
                cmd_systs[_id] = (phot_syst, filters, lambdas)

    return cmd_systs


if __name__ == '__main__':
    main()
