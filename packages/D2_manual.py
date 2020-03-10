
import pickle
from .out import make_D2_plot


def main(path):
    """
    This file is used to manually re-generate the D2 plot through a
    '_D2.pickle' file.

    To do this, create a file at the top level of the package (next to
    'asteca.py') with the following content:

    from packages import D2_manual
    path = './output/sixteen_clusts/cluster_folder/cluster_file'
    D2_manual.main(path)

    where 'path' is the path to the pickle file.
    """

    fname = path + '_D2.pickle'
    with open(fname, 'rb') as f:
        data = pickle.load(f)
    make_D2_plot.main(*data)
