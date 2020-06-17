
from . import BayesPlots


def pl_KP_Bys(
    gs, coord, kp_nburn, KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS,
    KP_samples, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta,
        KP_Bayes_kde):
    """
    """

    # Convert from deg to arcmin if (ra,dec) were used.
    if coord == 'deg':
        KP_samples[:, :, :2] = KP_samples[:, :, :2] * 60.
        KP_Bys_rc, KP_Bys_rt = KP_Bys_rc * 60., KP_Bys_rt * 60.
        KP_Bayes_kde[0][0], KP_Bayes_kde[1][0] =\
            KP_Bayes_kde[0][0] * 60., KP_Bayes_kde[1][0] * 60.
        coord2, dec_places = 'arcmin', "{:.2f}"
    else:
        coord2, dec_places = 'px', "{:.0f}"

    gsy, gsx = (0, 2), (0, 2)
    BayesPlots.autocorr(
        gs, gsx, gsy, KP_steps, KP_tau_autocorr, KP_ESS)

    gsy, gsx = (2, 4), (0, 2)
    BayesPlots.meanAF(gs, gsx, gsy, KP_steps, KP_mean_afs)

    #
    gsy, gsx = (0, 1), (2, 6)
    xylabel = r"$r_{{c}}$ [{}]".format(coord2)
    BayesPlots.traceplot(
        gs, gsx, gsy, KP_samples[:, :, 0], KP_Bys_rc, kp_nburn, xylabel,
        False)
    gsy, gsx = (1, 2), (2, 6)
    xylabel = r"$r_{{t}}$ [{}]".format(coord2)
    BayesPlots.traceplot(
        gs, gsx, gsy, KP_samples[:, :, 1], KP_Bys_rt, kp_nburn, xylabel)

    # Core vs tidal radii
    gsy, gsx = (4, 6), (0, 2)
    xylabel = r"$r_{{c}}$ [{}]".format(coord2)
    BayesPlots.histogram(
        gs, gsx, gsy, KP_samples[:, :, 0], KP_Bys_rc, KP_Bayes_kde[0],
        xylabel, dec_places)
    #
    gsy, gsx = (4, 6), (2, 4)
    xylabel = r"$r_{{t}}$ [{}]".format(coord2)
    BayesPlots.histogram(
        gs, gsx, gsy, KP_samples[:, :, 1], KP_Bys_rt, KP_Bayes_kde[1],
        xylabel, dec_places)
    #
    gsy, gsx = (4, 6), (4, 6)
    xylabel = (
        r"$r_{{c}}$ [{}]".format(coord2), r"$r_{{t}}$ [{}]".format(coord2))
    x_samples, y_samples = KP_samples[:, :, 0], KP_samples[:, :, 1]
    BayesPlots.twoParDens(
        gs, gsx, gsy, x_samples, y_samples, KP_Bys_rc, KP_Bys_rt, xylabel)

    try:
        gsy, gsx = (2, 3), (2, 6)
        xylabel = r"$ecc$"
        BayesPlots.traceplot(
            gs, gsx, gsy, KP_samples[:, :, 2], KP_Bys_ecc, kp_nburn, xylabel)
        gsy, gsx = (3, 4), (2, 6)
        xylabel = r"$\theta$ [rad]"
        BayesPlots.traceplot(
            gs, gsx, gsy, KP_samples[:, :, 3], KP_Bys_theta, kp_nburn, xylabel)

        # Eccentricity vs theta
        gsy, gsx = (6, 8), (0, 2)
        xylabel = r"$ecc$"
        BayesPlots.histogram(
            gs, gsx, gsy, KP_samples[:, :, 2], KP_Bys_ecc, KP_Bayes_kde[2],
            xylabel, dec_places)
        #
        gsy, gsx = (6, 8), (2, 4)
        xylabel = r"$\theta$ [rad]"
        BayesPlots.histogram(
            gs, gsx, gsy, KP_samples[:, :, 3], KP_Bys_theta, KP_Bayes_kde[3],
            xylabel, dec_places)
        #
        gsy, gsx = (6, 8), (4, 6)
        xylabel = (r"$ecc$", r"$\theta$ [rad]")
        x_samples, y_samples = KP_samples[:, :, 2], KP_samples[:, :, 3]
        BayesPlots.twoParDens(
            gs, gsx, gsy, x_samples, y_samples, KP_Bys_ecc, KP_Bys_theta, xylabel)
    except:
        pass


def plot(N, *args):
    """
    Handle each plot separately.
    """

    plt_map = {
        0: [pl_KP_Bys, 'King profile Bayes plots']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except Exception:
        print("  WARNING: error when plotting {}".format(plt_map.get(N)[1]))
        import traceback
        print(traceback.format_exc())
