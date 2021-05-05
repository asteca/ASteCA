
from . import BayesPlots


def pl_KP_Bys(
    gs, coord, kp_nburn, KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc,
        KP_Bys_theta):
    """
    """

    # Convert from deg to arcmin if (ra,dec) were used.
    KP_samples, KP_kde = KP_plot['KP_samples'], KP_plot['KP_kde']
    KP_samples[:, :, :2] = KP_samples[:, :, :2] * 60.
    KP_Bys_rc, KP_Bys_rt = KP_Bys_rc * 60., KP_Bys_rt * 60.
    KP_kde[0][0], KP_kde[1][0] = KP_kde[0][0] * 60., KP_kde[1][0] * 60.
    coord2, dec_places = 'arcmin', "{:.2f}"

    # Decide how to accommodate the plots
    if KP_samples.shape[-1] in (2, 3):
        gsy_AC, gsy_AF = (0, 1), (1, 2)
        gsy_Hrc, gsy_Hrt, gsy_Hrcrt = (2, 4), (2, 4), (2, 4)
    elif KP_samples.shape[-1] == 4:
        gsy_AC, gsy_AF = (0, 2), (2, 4)
        gsy_Hrc, gsy_Hrt, gsy_Hrcrt = (4, 6), (4, 6), (4, 6)

    gsx = (0, 2)
    BayesPlots.autocorr(
        gs, gsx, gsy_AC, KP_plot['KP_steps'], KP_plot['KP_tau_autocorr'],
        KP_plot['KP_ESS'])
    gsx = (0, 2)
    BayesPlots.meanAF(
        gs, gsx, gsy_AF, KP_plot['KP_steps'], KP_plot['KP_mean_afs'])

    # Number of burn-in samples
    nburn = int(KP_samples.shape[0] * kp_nburn)

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
    gsx = (0, 2)
    xylabel = r"$r_{{c}}$ [{}]".format(coord2)
    BayesPlots.histogram(
        gs, gsx, gsy_Hrc, KP_samples[nburn:, :, 0], KP_Bys_rc, KP_kde[0],
        xylabel, dec_places)
    #
    gsx = (2, 4)
    xylabel = r"$r_{{t}}$ [{}]".format(coord2)
    BayesPlots.histogram(
        gs, gsx, gsy_Hrt, KP_samples[nburn:, :, 1], KP_Bys_rt, KP_kde[1],
        xylabel, dec_places)
    #
    gsx = (4, 6)
    xylabel = (
        r"$r_{{c}}$ [{}]".format(coord2), r"$r_{{t}}$ [{}]".format(coord2))
    x_samples, y_samples = KP_samples[nburn:, :, 0], KP_samples[nburn:, :, 1]
    BayesPlots.twoParDens(
        gs, gsx, gsy_Hrcrt, x_samples, y_samples, KP_Bys_rc, KP_Bys_rt,
        xylabel)

    # Eccentricity, theta
    if KP_samples.shape[-1] == 4:
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
            gs, gsx, gsy, KP_samples[nburn:, :, 2], KP_Bys_ecc,
            KP_kde[2], xylabel, "{:.2f}")
        #
        gsy, gsx = (6, 8), (2, 4)
        xylabel = r"$\theta$ [rad]"
        BayesPlots.histogram(
            gs, gsx, gsy, KP_samples[nburn:, :, 3], KP_Bys_theta,
            KP_kde[3], xylabel, "{:.2f}")
        #
        gsy, gsx = (6, 8), (4, 6)
        xylabel = (r"$ecc$", r"$\theta$ [rad]")
        x_samples, y_samples = KP_samples[nburn:, :, 2],\
            KP_samples[nburn:, :, 3]
        BayesPlots.twoParDens(
            gs, gsx, gsy, x_samples, y_samples, KP_Bys_ecc, KP_Bys_theta,
            xylabel)


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
