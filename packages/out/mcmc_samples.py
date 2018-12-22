
from astropy.table import Table
from astropy.io import ascii


def main(clp, pd, mcmc_file_out, **kwargs):
    '''
    Create output data file with the MCMC sampler values for each parameter,
    for all chains/walkers. Maximum size of the output file is 10 Mbs.
    '''
    if pd['best_fit_algor'] in ('emcee', 'ptemcee'):
        varIdxs = clp['isoch_fit_params']['varIdxs']
        chains_nruns = clp['isoch_fit_params']['pars_chains']

        # Assume that each value stored occupies 10 bytes. Use a maximum file
        # size of 20 Mbs.
        max_sz = 20. * 1024. * 1024.
        ndim, nwalkers, nsteps = chains_nruns.shape
        # Maximum number of steps to store.
        max_n = int(max_sz / (ndim * nwalkers * 10.))

        tt, fmt = Table(), {}
        params = ['metal', 'log(age)', 'E_BV', 'dm', 'mass', 'binar_f']
        for i, par in enumerate(params):
            if i in varIdxs:
                c_model = varIdxs.index(i)
                for w in range(nwalkers):
                    bc = chains_nruns[c_model][w][-max_n:]
                    col_n = par + '_' + str(w)
                    tt[col_n], fmt[col_n] = bc, '%.5f'

        ascii.write(tt, mcmc_file_out, overwrite=True, formats=fmt)

        print('MCMC samples saved to file.')
