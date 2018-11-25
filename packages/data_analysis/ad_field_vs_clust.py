
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

    # Skip test if < 10 members are found within the cluster's radius.
    flag_few_members = False if len(clp['cl_region_c']) > 10 else True

    flag_ad_test = False
    ad_cl, ad_fr, pv_cl, pv_fr = [[[], []] for _ in range(4)]
    ad_cl_fr_p, ad_cl_fr_pk = [], []

    # Check if test is to be applied or skipped. Check if field regions
    # where found.
    if pd['ad_runs'] <= 0:
        print('Skipping field vs cluster A-D test.')

    elif clp['flag_no_fl_regs_c']:
        print('No field regions. Skipping field vs cluster A-D test.')

    elif flag_few_members:
        print('  WARNING: < 10 stars in cluster region.'
              '  Skipping field vs cluster A-D test.')

    else:
        print("    A-D test ({})".format(pd['ad_runs']))
        flag_ad_test = True

        run_total = 2. * int(pd['ad_runs'] * len(clp['field_regions_c']))
        runs = 0
        # Run first only for photometric data, and then for all data (if more
        # data exists)
        for i in range(2):
            for run_num in range(pd['ad_runs']):

                data_cl = dataExtract(clp['cl_region_c'], pd['ad_k_comb'], i)
                # Field regions
                data_fr = []
                for fr in clp['field_regions_c']:
                    data_fr.append(dataExtract(fr, pd['ad_k_comb'], i))

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

        ad_cl_fr_p = kdeplot(pvals_cl[0], pvals_fr[0], 'phot')
        data_id = 'phot+Plx+PM' if pd['ad_k_comb'] else 'Plx+PM'
        ad_cl_fr_pk = kdeplot(pvals_cl[1], pvals_fr[1], data_id)

    clp.update({
        'flag_ad_test': flag_ad_test, 'ad_cl': ad_cl, 'ad_fr': ad_fr,
        'ad_cl_fr_p': ad_cl_fr_p, 'ad_cl_fr_pk': ad_cl_fr_pk})
    return clp


def dataExtract(region, kin_flag, idx):
    """
    """
    def photData():
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
        return np.concatenate([mags, cols])

    def kinData():
        # Plx + pm_ra + pm_dec
        kins = np.array(list(zip(*list(zip(*region))[7])))[:3]
        e_kin = np.array(list(zip(*list(zip(*region))[8])))[:3]
        k_err = []
        for i, k in enumerate(kins):
            # Only process if any star contains at least one not 'nan'
            # data.
            if np.any(~np.isnan(k)):
                k_err.append(normErr(k, e_kin[i]))
        return k_err

    if idx == 0:
        data_all = photData()

    elif idx == 1:
        k_err = kinData()

        if kin_flag is True:
            phot_data = photData()
            # If any of the parallax+PMs dimensions was processed, add it.
            if k_err:
                data_all = np.concatenate([phot_data, np.array(k_err)])
            else:
                data_all = phot_data
                print("  WARNING: no valid Plx and/or PM data found.")

        elif kin_flag is False:
            if k_err:
                data_all = np.array(k_err)
            else:
                data_all = np.array([np.nan])
                print("  WARNING: no valid Plx and/or PM data found.")

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

            # anderson_darling_k([dd, data_y[i]])

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


def kdeplot(p_vals_cl, p_vals_f, data_id):
    # Define KDE limits.
    pmin = np.min(np.concatenate([p_vals_cl, p_vals_f]))
    pmax = np.max(np.concatenate([p_vals_cl, p_vals_f]))

    if pmin < pmax:
        xmin = max(-1., pmin)
        xmax = min(2., pmax)
        xrng = (xmax - xmin) * .3
        xmin = xmin - xrng
        xmax = xmax + xrng
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
                overlap = quad(y_pts, xmin, xmax)
            # Store y values for plotting the overlap filled.
            y_over = [float(y_pts(x_pt)) for x_pt in x_kde]

            # Probability value for the cluster.
            prob_cl_kde = 1 - overlap[0]
        else:
            # If not, assign probability as 1 minus the average of the cluster
            # vs the single field region used.
            prob_cl_kde = 1. - np.mean(p_vals_cl)
            # Pass empty lists for plotting.
            kde_f_1d, y_over = np.array([]), []
    else:
        print("  WARNING could not obtain p-value distribution" +
              " for '{}' data".format(data_id))
        prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde, y_over =\
            np.nan, np.array([]), np.array([]), [], []

    return p_vals_cl, p_vals_f, prob_cl_kde, kde_cl_1d, kde_f_1d, x_kde,\
        y_over


def anderson_darling_k(samples):
    """Apply the Anderson-Darling k-sample test.
    This test evaluates whether it is plausible that all the samples are drawn
    from the same distribution, based on Scholz and Stephens 1987. The
    statistic computed is their A_kn (rather than A_akn, which differs in
    how it handles ties). The significance of the result is computed by
    producing a scaled and standardized result T_kN, which is returned and
    compared against a list of standard significance levels. The next-larger
    p-value is returned.
    """

    _ps = np.array([0.25, 0.10, 0.05, 0.025, 0.01])
    # allow interpolation to get above _tm
    _b0s = np.array([0.675, 1.281, 1.645, 1.960, 2.326])
    _b1s = np.array([-0.245, 0.250, 0.678, 1.149, 1.822])
    _b2s = np.array([-0.105, -0.305, -0.362, -0.391, -0.396])

    samples = [np.array(sorted(s)) for s in samples]
    all = np.concatenate(samples + [[np.inf]])

    values = np.unique(all)
    L = len(values) - 1
    fij = np.zeros((len(samples), L), dtype=np.int)
    H = 0
    for (i, s) in enumerate(samples):
        c, be = np.histogram(s, bins=values) #, new=True)

        assert np.sum(c) == len(s)

        fij[i, :] = c
        H += 1. / len(s)

    ni = np.sum(fij, axis=1)[:, np.newaxis]
    N = np.sum(ni)
    k = len(samples)
    lj = np.sum(fij, axis=0)
    Mij = np.cumsum(fij, axis=1)
    Bj = np.cumsum(lj)

    A2 = np.sum(
        ((1. / ni) * lj / float(N) * (N * Mij - ni * Bj)**2 /
            (Bj * (N - Bj)))[:, :-1])

    h = np.sum(1. / np.arange(1, N))

    i = np.arange(1, N, dtype=np.float)[:, np.newaxis]
    j = np.arange(1, N, dtype=np.float)
    g = np.sum(np.sum((i < j) / ((N - i) * j)))

    a = (4 * g - 6) * (k - 1) + (10 - 6 * g) * H
    b = (2 * g - 4) * k**2 + 8 * h * k + \
        (2 * g - 14 * h - 4) * H - 8 * h + 4 * g - 6
    c = (6 * h + 2 * g - 2) * k**2 + \
        (4 * h - 4 * g + 6) * k + (2 * h - 6) * H + 4 * h
    d = (2 * h + 6) * k**2 - 4 * h * k

    sigmaN2 = (a * N**3 + b * N**2 + c * N + d) / ((N - 1) * (N - 2) * (N - 3))

    sigmaN = np.sqrt(sigmaN2)

    TkN = (A2 - (k - 1)) / sigmaN

    tkm1 = _b0s + _b1s / np.sqrt(k - 1) + _b2s / (k - 1)

    ix = np.searchsorted(tkm1, TkN)
    if ix > 0:
        p = _ps[ix - 1]
    else:
        p = 1.
    return A2, TkN, (tkm1, _ps.copy()), p
