
# DEPRECATED 31/12/18


def fundParams(clp):
    '''
    Round fundamental fitted parameters and their errors.
    '''
    cp_e, lenght = [], [5, 4, 3, 3, 0, 2]
    for i, e in enumerate(clp['isoch_fit_errors']):
        cp_e.append(round(e, lenght[i]))
    clp['fit_errors_r'] = ['{:g}'.format(_) for _ in cp_e]

    # Round mean solutions.
    clp['fit_params_mean_r'] = [
        '{:g}'.format(_) for _ in clp['isoch_fit_params']['mean_sol']]

    # Round median solutions.
    clp['fit_params_median_r'] = [
        '{:g}'.format(_) for _ in clp['isoch_fit_params']['median_sol']]

    # Round MAP/ML solutions.
    clp['fit_params_map_r'] = [
        '{:g}'.format(_) for _ in clp['isoch_fit_params']['map_sol']]

    return clp
