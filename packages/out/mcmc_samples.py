
from astropy.table import Table
from astropy.io import ascii


def main(clp, pd, mcmc_file_out, **kwargs):
    '''
    Create output data file with the "best" chain of the MCMC sampler (smallest
    acorr time), for each parameter.
    Store only the latest 10000 steps or less.
    '''
    if pd['best_fit_algor'] in ('emcee', 'ptemcee'):
        varIdxs = clp['isoch_fit_params']['varIdxs']
        chains_nruns = clp['isoch_fit_params']['pars_chains']
        min_at_c = clp['isoch_fit_params']['min_at_c']

        tt, fmt = Table(), {}
        params = ['metal', 'log(age)', 'E_BV', 'dm', 'mass', 'bf']
        # TODO better column names
        for i, par in enumerate(params):
            if i in varIdxs:
                c_model = varIdxs.index(i)
                bc = chains_nruns[c_model][min_at_c[c_model]][-10000:]
                tt[par], fmt[par] = bc, '%.5f'

        ascii.write(tt, mcmc_file_out, overwrite=True, formats=fmt)

        print('MCMC samples saved to file.')
