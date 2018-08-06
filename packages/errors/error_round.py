

def fundParams(clp):
    '''
    Round fundamental fitted parameters and their errors.
    '''
    cp_e, lenght = [], [5, 3, 3, 3, 0, 2]
    for i, e in enumerate(clp['isoch_fit_errors']):
        cp_e.append(round(e, lenght[i]))
    clp['fit_errors_r'] = ['{:g}'.format(_) for _ in cp_e]
    clp['fit_params_r'] = [
        '{:g}'.format(_) for _ in clp['isoch_fit_params']['best_sol']]

    return clp
