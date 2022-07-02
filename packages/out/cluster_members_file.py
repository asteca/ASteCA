
from ..inp import data_IO


def main(clp, npd, id_ids, **kwargs):
    """
    Create output data file with membership probabilities and selection
    identifiers:

    -1 means outside of the cluster region
    0 means inside but not selected as member
    1 means selected as member
    """

    all_data = data_IO.dataRead(None, npd['data_file'], 'r')
    if 'sel' in all_data.keys():
        all_data.remove_column('sel')
    if 'memb_probs' in all_data.keys():
        all_data.remove_column('memb_probs')
    all_data.add_column(-1, name='sel')
    all_data.add_column(0., name='memb_probs')
    IDs = list(all_data[id_ids].astype(str))

    # Add ID associated to the use of the each star in the fundamental
    # parameters estimation process and MPs.
    for i, reg in enumerate([clp['cl_reg_no_fit'], clp['cl_reg_fit']]):
        for st in reg:
            # st[9] == membership probability
            try:
                j = IDs.index(st[0])
                all_data['sel'][j] = i
                all_data['memb_probs'][j] = st[9]
            except ValueError:
                pass

    data_IO.dataSave(all_data, npd['memb_file_out'], 'w')
    print("Cluster region and MPs saved to file")
