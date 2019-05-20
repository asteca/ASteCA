
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import uniform_filter1d


def xxx():
    """
    """
    ax = plt.subplot(gs[0:2, 2:8])


def pl_MAP_lkl(dummy, gs, prob_mean, map_lkl, map_lkl_final):
    '''
    Evolution of MAP likelihood values.
    '''
    ax = plt.subplot(gs[2:4, 4:6])
    x, y = list(zip(*map_lkl))
    ax.plot(x, y, label=r"$L_{{min}}={:.1f}$".format(map_lkl_final))
    x, y = list(zip(*prob_mean))
    ax.plot(x, y, label="Mean LP")
    plt.xlabel("steps", fontsize=14)
    plt.ylabel("Lkl (MAP)", fontsize=14)
    ax.legend(fontsize='small', loc=0)  # , handlelength=0.)


def pl_param_chain(
    par_name, gs, best_fit_algor, mean_sol, min_max_p, nwalkers, nburn, nsteps,
        model_done, varIdxs, pre_bi, post_bi, acorr_t, med_at_c,
        mcmc_ess):
    '''
    Parameter sampler chain.
    '''
    plot_dict = {
        'metal': [8, 12, 0, 1, 0], 'age': [8, 12, 1, 2, 1],
        'ext': [8, 12, 2, 3, 2], 'dist': [8, 12, 3, 4, 3],
        'mass': [8, 12, 4, 5, 4], 'binar': [8, 12, 5, 6, 5]
    }

    labels = [r'$z$', r'$\log(age)$', r'$E_{{(B-V)}}$', r'$(m-M)_o$',
              r'$M\,(M_{{\odot}})$', r'$b_{{frac}}$']

    gs_x1, gs_x2, gs_y1, gs_y2, cp = plot_dict[par_name]
    ax = plt.subplot(gs[gs_y1:gs_y2, gs_x1:gs_x2])
    if cp == 5:
        plt.xlabel("Steps")
    else:
        ax.tick_params(labelbottom=False)
    plt.ylabel(labels[cp])
    if best_fit_algor in ('ptemcee', 'emcee'):
        N_bi, N_tot = nburn, nburn + nsteps
    elif best_fit_algor == 'abc':
        N_bi, N_tot = nburn, nsteps
    ax.set_xlim(0, N_tot)

    if cp in varIdxs:
        c_model = varIdxs.index(cp)

        # Chain with median Tau
        # Burn-in stage
        pre_bi_max_at = pre_bi[c_model][med_at_c[c_model]]
        plt.plot(range(N_bi), pre_bi_max_at, c='grey', lw=.5, alpha=0.5)
        # Post burn-in.
        post_bi_max_at = post_bi[c_model][med_at_c[c_model]]
        plt.plot(np.arange(N_bi, N_tot), post_bi_max_at, c='k', lw=.8,
                 ls='-', alpha=0.5)
        # # Best chain
        # # Burn-in stage
        # pre_bi_min_at = pre_bi[c_model][min_at_c[c_model]]
        # plt.plot(range(N_bi), pre_bi_min_at, c='grey', lw=.5, alpha=0.5)
        # # Post burn-in.
        # post_bi_min_at = post_bi[c_model][min_at_c[c_model]]
        # plt.plot(
        #     np.arange(N_bi, N_tot), post_bi_min_at, c='k', lw=.8, alpha=0.5)

        # Filtered mean of all chains.
        N = post_bi.shape[-1]
        xavr = uniform_filter1d(
            np.mean(post_bi[c_model], axis=0), int(.02 * N))
        plt.plot(np.arange(N_bi, N_tot), xavr, c='g')

        # Mean
        plt.axhline(
            y=float(mean_sol[cp]), linestyle='--', color='blue', zorder=4)
        #  16th and 84th percentiles (1 sigma) around median.
        ph = np.percentile(model_done[c_model], 84)
        pl = np.percentile(model_done[c_model], 16)
        plt.axhline(y=ph, linestyle=':', color='orange', zorder=4)
        plt.axhline(y=pl, linestyle=':', color='orange', zorder=4)
        # plt.axhline(
        #     y=float(cp_r[cp]), color='k', ls='--', lw=1.2, zorder=4,
        #     label=r"$\tau={:.0f}\;(\hat{{n}}_{{eff}}={:.0f})$".format(
        #         acorr_t[c_model], mcmc_ess[c_model]))

        ax.set_title(
            r"$\hat{{\tau}}_{{c}}={:.0f}\;(\hat{{n}}_{{eff}}={:.0f})$".format(
                acorr_t[c_model], mcmc_ess[c_model]))
        ax.set_ylim(min_max_p[cp][0], min_max_p[cp][1])
        # ax.legend(fontsize='small', loc=0, handlelength=0.)


def pl_betas(dummy, gs, best_fit_algor, betas_pt):
    '''
    Evolution of Temp swaps AFs.
    '''
    ax = plt.subplot(gs[0:2, 2:4])
    x, betas = betas_pt
    Nt = len(betas) - 1
    ax.set_title(r"$N_{{temps}}={}$".format(Nt + 1), fontsize=10)
    for i, y in enumerate(betas):
        if i == 0:
            lbl_ls = ("Cold", '--', 1.5, 4)
        elif i == Nt:
            lbl_ls = ("Hot", ':', 1.5, 4)
        else:
            lbl_ls = (None, '-', .5, 1)
        ax.plot(
            x, y[:len(x)], label=lbl_ls[0], ls=lbl_ls[1], lw=lbl_ls[2],
            zorder=lbl_ls[3])
    plt.xlabel("steps", fontsize=14)
    plt.ylabel(r"$\beta$", fontsize=14)
    ax.legend(fontsize='small', loc=0)


def pl_Tswaps(dummy, gs, best_fit_algor, tswaps_afs):
    '''
    Evolution of Temp swaps AFs.
    '''
    ax = plt.subplot(gs[0:2, 4:6])
    # ax.set_title(
    #     r"$MAF_{{[T=1]}}={:.3f}$".format(maf_steps[1][0][-1]), fontsize=10)
    x, y_replicas = tswaps_afs
    Nt = len(y_replicas) - 1
    for i, y in enumerate(y_replicas):
        if i == 0:
            lbl_ls = ("Cold", '--', 1.5, 4)
        elif i == Nt:
            lbl_ls = ("Hot", ':', 1.5, 4)
        else:
            lbl_ls = (None, '-', .5, 1)
        ax.plot(
            x, y, label=lbl_ls[0], ls=lbl_ls[1], lw=lbl_ls[2],
            zorder=lbl_ls[3])
    plt.xlabel("steps", fontsize=14)
    plt.ylabel("Tswaps AF", fontsize=14)
    ax.legend(fontsize='small', loc=0)


def pl_MAF(dummy, gs, best_fit_algor, maf_steps):
    '''
    Evolution of MAF values.
    '''
    ax = plt.subplot(gs[0:2, 6:8])
    ax.set_title(
        r"$MAF_{{[T=1]}}={:.3f}$".format(maf_steps[1][0][-1]), fontsize=10)
    # x, y = list(zip(*maf_steps))
    x, y_replicas = maf_steps
    Nt = len(y_replicas) - 1
    for i, y in enumerate(y_replicas):
        if i == 0:
            lbl_ls = ("Cold", '--', 1.5, 4)
        elif i == Nt:
            lbl_ls = ("Hot", ':', 1.5, 4)
        else:
            lbl_ls = (None, '-', .5, 1)
        ax.plot(
            x, y, label=lbl_ls[0], ls=lbl_ls[1], lw=lbl_ls[2],
            zorder=lbl_ls[3])
    plt.xlabel("steps", fontsize=14)
    plt.ylabel("MAF", fontsize=14)
    if best_fit_algor in ('ptemcee', 'emcee'):
        plt.axhline(y=.25, color='grey', ls=':', lw=1.2, zorder=4)
        plt.axhline(y=.5, color='grey', ls=':', lw=1.2, zorder=4)
    ax.legend(fontsize='small', loc=0)  # , handlelength=0.)


def pl_mESS(dummy, gs, mESS, minESS, minESS_epsilon):
    '''
    mESS plot.
    '''
    ax = plt.subplot(gs[4:6, 6:8])
    plt.xlabel(r"$CI\;(=1-\alpha)$", fontsize=14)
    plt.ylabel(r"$\epsilon$", fontsize=16)
    ax.set_xlim(.01, 1.01)

    plt.plot(
        1. - np.array(minESS_epsilon[0]), minESS_epsilon[1],
        label="minESS ({:.0f})".format(minESS))
    plt.plot(
        1. - np.array(minESS_epsilon[0]), minESS_epsilon[2],
        label="mESS ({:.0f})".format(mESS))
    ax.legend(fontsize='small', loc=0)


def pl_tau_histo(dummy, gs, all_taus):
    """
    """
    ax = plt.subplot(gs[4:6, 6:8])
    plt.title(r"$\tau$ for all chains and parameters (N={})".format(
        len(all_taus)), fontsize=10)
    plt.hist(np.nan_to_num(all_taus), bins=20, color='#64B1D7')
    plt.axvline(
        np.nanmean(all_taus), color="r",
        label=r"$\hat{{\tau}}_{{c\,,\,p}}$={:.0f}".format(
            np.nanmean(all_taus)))
    plt.xlabel(r"$\tau$", fontsize=14)
    plt.ylabel(r"$N$", fontsize=14)
    ax.legend(fontsize='small', loc=0)


def pl_tau(dummy, gs, N_steps_conv, N_conv, tol_conv, tau_index, tau_autocorr):
    '''
    Tau vs steps plot.
    '''
    ax = plt.subplot(gs[6:8, 10:12])
    plt.title(r"$N_{{conv}}={:.0f}, tol_{{conv}}={:.2f}$".format(
        N_conv, tol_conv), fontsize=10)
    plt.xlabel("steps", fontsize=14)
    plt.ylabel(r"$\hat{\tau}_{{c\,|\,p}}$", fontsize=14)

    n = N_steps_conv * np.arange(1, tau_index + 1)
    plt.plot(n, n / 50., "--g", label="N/50")
    plt.plot(n, n / 100., "--b", label="N/100")
    plt.plot(n, n / 500., "--r", label="N/500")
    plt.plot(n, tau_autocorr)
    plt.xlim(0, n.max())
    try:
        plt.ylim(
            0, np.nanmax(tau_autocorr) + 0.1 *
            (np.nanmax(tau_autocorr) - np.nanmin(tau_autocorr)))
    except ValueError:
        print("  WARNING: no mean autocorrelation values to plot.")
    ax.legend(fontsize='small', loc=0)


def pl_lags(dummy, gs, varIdxs, lag_zero, acorr_function):
    '''
    lags plot.
    '''
    ax = plt.subplot(gs[6:8, 8:10])
    plt.xlabel("Lag", fontsize=14)
    plt.ylabel("ACF", fontsize=14)

    plot_dict = ['metal', 'age', 'ext', 'dist', 'mass', 'binar']
    for i, par_name in enumerate(plot_dict):
        if i in varIdxs:
            c_model = varIdxs.index(i)
            p = acorr_function[c_model]
            plt.plot(
                range(len(p)), p, lw=.8, alpha=0.5,
                label="{}".format(par_name))
    ax.legend(fontsize='small', loc=0)
    plt.plot((lag_zero, lag_zero), (0, .5), ls='--', c='k')
    ax.text(lag_zero, .52, "{}".format(lag_zero))


def pl_GW(dummy, gs, varIdxs, geweke_z):
    '''
    Geweke plot.
    '''
    ax = plt.subplot(gs[8:10, 10:12])
    ax.set_title("Geweke", fontsize=10)
    plt.xlabel("First iteration in segment", fontsize=14)
    plt.ylabel("z-score", fontsize=14)
    plt.axhline(y=2., color='grey', ls=':', lw=1.2, zorder=4)
    plt.axhline(y=-2., color='grey', ls=':', lw=1.2, zorder=4)

    plot_dict = ['metal', 'age', 'ext', 'dist', 'mass', 'binar']
    ymin, ymax = [], []
    for i, par_name in enumerate(plot_dict):
        if i in varIdxs:
            c_model = varIdxs.index(i)
            p = geweke_z[c_model]
            idx, zscore = list(zip(*p))
            plt.plot(
                idx, zscore, ls="-.", linewidth=.7,
                label="{}".format(par_name))
            ymin.append(np.nanmin(zscore))
            ymax.append(np.nanmax(zscore))

    if ymin and ymax:
        ax.set_ylim(max(-2.1, min(ymin)), min(2.1, max(ymax)))
    ax.legend(fontsize='small', loc=0)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''
    plt.style.use('seaborn-darkgrid')
    plt_map = {
        0: [pl_MAP_lkl, ' MAP likelihood values'],
        1: [pl_MAF, ' MAF vs steps'],
        2: [pl_betas, ' Betas vs steps'],
        3: [pl_Tswaps, ' Tswaps AFs vs steps'],
        4: [pl_param_chain, args[0] + ' sampler chain'],
        5: [pl_tau, args[0]],
        6: [pl_mESS, args[0]],
        7: [pl_lags, args[0]],
        8: [pl_GW, args[0]],
        9: [pl_tau_histo, args[0]]
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except Exception:
        import traceback
        print(traceback.format_exc())
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))
    plt.style.use('default')
