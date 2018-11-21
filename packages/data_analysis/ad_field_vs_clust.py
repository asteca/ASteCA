
import numpy as np
from scipy.stats import anderson_ksamp, gaussian_kde


def main(pd, clp, cld_c):
    """

    AD test for k-samples: "tests the null hypothesis that k-samples are drawn
    from the same population"

    We want to REJECT the null hypothesis to state that the cluster is
    reasonably different from the field region.

    Significance level: probability of rejecting the null hypothesis when it is
    true (Type I Error). This is: the probability of being wrong saying that
    the samples are drawn from different populations. Equivalent: the
    probability that both samples come from the same population.

    If AD test value > 4 then I can say that the samples come from different
    populations (reject null hypothesis). Ideally I want:

    * large AD values for cluster vs field regions
    * small AD values for field vs field regions
    """

    # TODO incorporate to params_input.dat
    pd['flag_kde_test'] = False

    # Skip test if < 10 members are found within the cluster's radius.
    flag_few_members = False if len(clp['cl_region_c']) > 10 else True

    # Check if test is to be applied or skipped. Check if field regions
    # where found.
    flag_kde_test, kde_test_params = False, []
    if not pd['flag_kde_test']:
        print('Skipping field vs cluster KDE test for cluster.')

    elif clp['flag_no_fl_regs_c']:
        print('No field regions. Skipping field vs cluster KDE test.')

    elif flag_few_members:
        print('  WARNING: < 10 stars in cluster region.'
              '  Skipping field vs cluster KDE test.')

    # Run process.
    else:

        data_cl = dataExtract(clp['cl_region_c'])
        data_fr = []
        for fr in clp['field_regions_c']:
            data_fr.append(dataExtract(fr))

        ad_vals_cl, ad_vals_f = [], []
        for f_idx, data_fl in enumerate(data_fr):

            ad_vals_cl.append(ADtest(data_cl, data_fl))

            # Compare the field region used above with all the remaining
            # field regions. This results in [N*(N+1)/2] combinations of
            # field vs field comparisons.
            # grid_vals = gridVals(data_fl)
            for data_fl2 in data_fr[(f_idx + 1):]:
                ad_vals_f.append(ADtest(data_fl, data_fl2))

        # Extract AD values and capped pvalues
        a2vals_cl = np.array(ad_vals_cl)[:, :, 0].ravel()
        a2vals_f = np.array(ad_vals_f)[:, :, 0].ravel()
        # Cap p_values
        pvals_cl = pvalFix(np.array(ad_vals_cl).reshape(a2vals_cl.shape[0], 2))
        pvals_f = pvalFix(np.array(ad_vals_f).reshape(a2vals_f.shape[0], 2))

        # import matplotlib.pyplot as plt
        # plt.hist(a2vals_cl, bins=25, alpha=.5, density=True, label='cl')
        # plt.hist(a2vals_f, bins=25, alpha=.5, density=True, label='fr')
        # plt.legend()
        # plt.show()

        kdeplot(pvals_cl, pvals_f)

        import pdb; pdb.set_trace()  # breakpoint f79d6346 //


    clp.update({
        'flag_kde_test': flag_kde_test, 'kde_test_params': kde_test_params})
    return clp


def dataExtract(region):
    """
    """
    # Main magnitude. Must have shape (1, N)
    mags = np.array(list(zip(*list(zip(*region))[3])))
    # One or two colors
    cols = np.array(list(zip(*list(zip(*region))[5])))

    # Check if Plx and/or PM data exist. TODO
    # Plx + pm_ra + pm_dec
    kins = np.array(list(zip(*list(zip(*region))[7])))[:3]

    data_all = np.concatenate([mags, cols, kins])

    return data_all


def ADtest(data_x, data_y):
    """
    Obtain Anderson-Darling test for each data dimension.
    """
    ad_vals = []
    # For each dimension
    for i, dd in enumerate(data_x):
        ad_stts = list(anderson_ksamp([dd, data_y[i]]))
        # Store A-D value and p-value.
        ad_vals.append([ad_stts[0], ad_stts[2]])

    return np.array(ad_vals)


def pvalFix(ad_vals):
    """
    TODO Fix taken from Scipy v1.2.0:
    https://github.com/scipy/scipy/blob/
    dfc9a9c73ced00e2588dd8d3ee03f9e106e139bf/scipy/stats/morestats.py#L2032
    """

    critical = np.array([0.325, 1.226, 1.961, 2.718, 3.752, 4.592, 6.546])
    sig = np.array([0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001])

    p_vals = []
    for A2, pval in ad_vals:

        if A2 < critical.min():
            # p-value capped
            p = sig.max()
        elif A2 > critical.max():
            # p-value floored
            p = sig.min()
        else:
            p = pval
        p_vals.append(p)

    return np.array(p_vals)


def kdeplot(p_vals_cl, p_vals_f):
    from scipy.integrate import quad

    # Define KDE limits.
    xmin, xmax = -1., 2.
    x_kde = np.mgrid[xmin:xmax:1000j]

    # Obtain the 1D KDE for the cluster region (stars inside cluster's
    # radius) vs field regions.
    kernel_cl = gaussian_kde(p_vals_cl)
    # KDE for plotting.
    kde_cl_1d = np.reshape(kernel_cl(x_kde).T, x_kde.shape)

    # Check if field regions were compared among each other.
    if p_vals_f.any():
        # Obtain the 1D KDE for the field regions vs field regions.
        kernel_f = gaussian_kde(p_vals_f)

        # KDE for plotting.
        kde_f_1d = np.reshape(kernel_f(x_kde).T, x_kde.shape)

        # Calculate overlap between the two KDEs.
        def y_pts(pt):
            y_pt = min(kernel_cl(pt), kernel_f(pt))
            return y_pt

        overlap = quad(y_pts, -1., 2.)
        # Store y values for plotting the overlap filled.
        y_over = [float(y_pts(x_pt)) for x_pt in x_kde]

        # Probability value for the cluster.
        prob_cl_kde = 1 - overlap[0]
    else:
        # If not, assign probability as 1 minus the average of the cluster
        # vs the single field region used.
        prob_cl_kde = 1. - np.mean(p_vals_cl)
        # Pass empty lists for plotting.
        kde_f_1d, y_over = np.asarray([]), []

    pl_p_vals(
        True, p_vals_cl, p_vals_f, prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde,
        y_over)


def pl_p_vals(
    flag_kde_test, p_vals_cl, p_vals_f, prob_cl_kde, kde_cl_1d, kde_f_1d,
        x_kde, y_over):
    '''
    Distribution of KDE p_values.
    '''
    import matplotlib.pyplot as plt
    if flag_kde_test:
        ax = plt.subplot(111)  # gs[4:6, 2:4])
        # plt.xlim(-0.1, 1.05)
        plt.xlim(-0.09, .4)
        plt.ylim(0, 1.02)
        plt.xlabel('p-values', fontsize=12)
        plt.ylabel('Density (normalized)', fontsize=12)
        ax.minorticks_on()
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
                zorder=1)
        # Grid to background.
        ax.set_axisbelow(True)
        # Plot field vs field KDE.
        if kde_f_1d.any():
            max_kde = max(max(kde_f_1d), max(kde_cl_1d))
            plt.plot(x_kde, kde_f_1d / max_kde, color='b', ls='-', lw=1.,
                     label=r'$KDE_{{fl}}\,({})$'.format(len(p_vals_f)),
                     zorder=2)
        else:
            max_kde = max(kde_cl_1d)
        # Plot cluster vs field KDE.
        plt.plot(x_kde, kde_cl_1d / max_kde, color='r', ls='-', lw=1.,
                 label=r'$KDE_{{cl}}\,({})$'.format(len(p_vals_cl)), zorder=2)
        # Fill overlap.
        if y_over:
            plt.fill_between(x_kde, np.asarray(y_over) / max_kde, 0,
                             color='grey', alpha='0.5')
        text = '$P_{cl}^{KDE} = %0.2f$' % round(prob_cl_kde, 2)
        plt.text(0.05, 0.92, text, transform=ax.transAxes,
                 bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
        # Legend.
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, labels, loc='upper right', numpoints=1,
                        fontsize=12)
        leg.get_frame().set_alpha(0.6)

        plt.show()
