
import json


def main(pd, npd):
    """
    Save 'pd' dictionary to .json file.
    """
    with open(npd['json_out'], "w") as jsonFile:
        pd2 = dict(pd)
        # Don't store these entries, they are large and unnecessary.
        for key in [
                'cmd_systs', 'plot_isoch_data', 'theor_tracks',
                'fundam_params']:
            try:
                pd2.pop(key)
            except KeyError:
                pass
        json.dump(pd2, jsonFile)
