
from .. import update_progress
import numpy as np
from scipy.stats import anderson_ksamp, gaussian_kde
from scipy.integrate import quad
import warnings


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
    pd['flag_kde_test'] = True

    # Skip test if < 10 members are found within the cluster's radius.
    flag_few_members = False if len(clp['cl_region_c']) > 10 else True

    flag_kde_test = False
    ad_cl, ad_fr, pv_cl, pv_fr = [[[], []] for _ in range(4)]
    ad_cl_fr_p, ad_cl_fr_pk = [], []

    # Check if test is to be applied or skipped. Check if field regions
    # where found.
    if not pd['flag_kde_test']:
        print('Skipping field vs cluster A-D test for cluster.')

    elif clp['flag_no_fl_regs_c']:
        print('No field regions. Skipping field vs cluster A-D test.')

    elif flag_few_members:
        print('  WARNING: < 10 stars in cluster region.'
              '  Skipping field vs cluster A-D test.')

    else:
        print("        A-D test")
        flag_kde_test = True

        error_runs = 100
        run_total, runs = 2. * int(error_runs * len(clp['field_regions_c'])), 0
        # Run first only for photometric data, and then for all data (if more
        # data exists)
        for i, kflag in enumerate((False, True)):
            for run_num in range(error_runs):

                data_cl = dataExtract(clp['cl_region_c'], kflag)
                # Field regions
                data_fr = []
                for fr in clp['field_regions_c']:
                    data_fr.append(dataExtract(fr, kflag))

                # Compare to each defined field region.
                for f_idx, data_fl in enumerate(data_fr):

                    ad_pv = ADtest(data_cl, data_fl)
                    ad_cl[i] += list(ad_pv[0])
                    pv_cl[i] += list(ad_pv[1])

                    # Compare the field region used above with all the
                    # remaining field regions. This results in [N*(N+1)/2]
                    # combinations of field vs field comparisons.
                    for data_fl2 in data_fr[(f_idx + 1):]:
                        ad_pv = ADtest(data_fl, data_fl2)
                        ad_fr[i] += list(ad_pv[0])
                        pv_fr[i] += list(ad_pv[1])

                    runs += 1
                update_progress.updt(run_total, runs)

        # Cap p_values
        pvals_cl = [pvalFix(ad_cl[0], pv_cl[0]), pvalFix(ad_cl[1], pv_cl[1])]
        pvals_fr = [pvalFix(ad_fr[0], pv_fr[0]), pvalFix(ad_fr[1], pv_fr[1])]

        ad_cl_fr_p = kdeplot(pvals_cl[0], pvals_fr[0])
        ad_cl_fr_pk = kdeplot(pvals_cl[1], pvals_fr[1])

    clp.update({
        'flag_kde_test': flag_kde_test, 'ad_cl': ad_cl, 'ad_fr': ad_fr,
        'ad_cl_fr_p': ad_cl_fr_p, 'ad_cl_fr_pk': ad_cl_fr_pk})
    return clp


def dataExtract(region, kin_flag):
    """
    """
    # Main magnitude. Must have shape (1, N)
    mags = np.array(list(zip(*list(zip(*region))[3])))
    e_mag = np.array(list(zip(*list(zip(*region))[4])))
    mags = normErr(mags, e_mag)

    # One or two colors
    cols = np.array(list(zip(*list(zip(*region))[5])))
    e_col = np.array(list(zip(*list(zip(*region))[6])))
    c_err = []
    for i, c in enumerate(cols):
        c_err.append(normErr(c, e_col[i]))
    cols = np.array(c_err)

    data_all = np.concatenate([mags, cols])

    if kin_flag:
        # Plx + pm_ra + pm_dec
        kins = np.array(list(zip(*list(zip(*region))[7])))[:3]
        e_kin = np.array(list(zip(*list(zip(*region))[8])))[:3]
        k_err = []
        for i, k in enumerate(kins):
            # Only process if any star contains at least one not 'nan' data.
            if np.any(~np.isnan(k)):
                k_err.append(normErr(k, e_kin[i]))
        # If any of the parallax+PMs dimensions was processed, add it.
        if k_err:
            data_all = np.concatenate([data_all, np.array(k_err)])

    return data_all


def normErr(x, e_x):
    # Randomly move mag and color through a Gaussian function.
    return x + np.random.normal(0, 1, len(x)) * e_x


def ADtest(data_x, data_y):
    """
    Obtain Anderson-Darling test for each data dimension.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        ad_vals = []
        # For each dimension
        for i, dd in enumerate(data_x):
            ad_stts = list(anderson_ksamp([dd, data_y[i]]))
            # Store A-D value and p-value.
            ad_vals.append([ad_stts[0], ad_stts[2]])

    return np.array(ad_vals).T


def pvalFix(ad_vals, p_vals):
    """
    TODO Fix taken from Scipy v1.2.0:
    https://github.com/scipy/scipy/blob/
    dfc9a9c73ced00e2588dd8d3ee03f9e106e139bf/scipy/stats/morestats.py#L2032
    """
    critical = np.array([0.325, 1.226, 1.961, 2.718, 3.752, 4.592, 6.546])
    sig = np.array([0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001])

    p_vals_c = []
    for i, A2 in enumerate(ad_vals):

        if A2 < critical.min():
            # p-value capped
            # if p_vals[i] > 1.:
            #     p = 1.
            # else:
            #     p = p_vals[i]
            p = sig.max()
        elif A2 > critical.max():
            # p-value floored
            p = sig.min()
        else:
            p = p_vals[i]

        p_vals_c.append(p)

    return np.array(p_vals_c)


def kdeplot(p_vals_cl, p_vals_f):
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

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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

    return p_vals_cl, p_vals_f, prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde,\
        y_over
