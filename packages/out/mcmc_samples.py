
from astropy.table import Table, vstack
from astropy.io import ascii


def main(clp, pd, mcmc_file_out, **kwargs):
    '''
    Create output data file with the 5 best chains of the MCMC sampler, for
    each parameter.
    Store only the latest 1000 steps or less.
    '''
    if pd['best_fit_algor'] == 'emcee':
        varIdxs = clp['isoch_fit_params']['varIdxs']
        chains_nruns = clp['isoch_fit_params']['pars_chains']
        min_at_5c = clp['isoch_fit_params']['min_at_5c']

        tt = Table()
        params = ['metal', 'log(age)', 'E_BV', 'dm', 'mass', 'bf']
        # TODO better column names
        for i, par in enumerate(params):
            if i in varIdxs:
                c_model = varIdxs.index(i)
                bc = chains_nruns[c_model][min_at_5c[c_model]][:, -1000:]
                tt = vstack([tt, Table(bc.T)])

        ascii.write(
            tt, mcmc_file_out, overwrite=True,
            formats={'col0': '%.5f', 'col1': '%.5f', 'col2': '%.5f',
                     'col3': '%.5f', 'col4': '%.5f'})

        print('MCMC samples saved to file.')
