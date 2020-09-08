
from ..inp import data_IO
from astropy.table import Table


def main(clp, npd, **kwargs):
    """
    Create output data file with stars inside the cluster radius along with
    their membership probabilities and 'clean region' selection identifier.
    """

    # Add ID associated to the use of the each star in the fundamental
    # parameters estimation process (ie: after cleaning the cluster region).
    data, idx = [], ['1', '0']
    for i, reg in enumerate([clp['cl_reg_fit'], clp['cl_reg_no_fit']]):
        for st in reg:
            # Identify stars selected by the removal function.
            data.append([st[0], st[9], idx[i]])

    # TODO: this block gets really slow for large clusters
    # Add "incomplete" data in cluster region to file.
    ids_data = list(zip(*data))[0]
    for i, st in enumerate(clp['cl_region_i']):
        if st[0] not in ids_data:
            # Identify stars selected by the removal function.
            data.append([
                st[0], round(clp['memb_probs_cl_region_i'][i], 2), '-1'])

    t = Table(list(zip(*data)), names=['ID', 'MP', 'sel'])
    data_IO.dataSave(t, npd['memb_file_out'])
    print("Cluster region and MPs saved to file")
