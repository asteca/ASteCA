
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt


def main(clp):

    # Synthetic photometric and mass data.
    # shape: (N: number of stars, M: number of mags & colors)
    phot_synth = np.array(clp['synth_clst'][0][0]).T
    m_ini = np.array(clp['synth_clst'][1][2])
    m_act = np.array(clp['synth_clst'][1][3])

    # Observed photometric data
    phot_obs = []
    for st in clp['cl_max_mag']:
        phot_obs.append(list(st[3]) + list(st[5]))
    phot_obs = np.array(phot_obs)

    # Find synthetic stars closest to all the observed stars.
    tree = spatial.cKDTree(phot_synth)
    mindist, minid = tree.query(phot_obs)
    # Mass associated to the observed stars.
    m_ini_obs, m_act_obs = m_ini[minid], m_act[minid]

    # Obtain IMF factor.
    h1, x = np.histogram(m_ini_obs)
    # Mask -inf values that can appear when a bin contains 0 elements.
    h1_log = np.log10(h1)
    msk = h1_log == -np.inf
    alpha_i, _ = np.polyfit(np.log10(x[:-1][~msk]), h1_log[~msk], 1)
    m, b = np.polyfit(x[:-1][~msk], h1_log[~msk], 1)
    plt.plot(
        x, 10**np.poly1d((m, b))(x), c='#7BC5D4',
        label=r"$\alpha_{{i}}={:.3f}$".format(alpha_i))
    plt.scatter(x[1:], h1, c='#7BC5D4', label=r'$M_{i}^{obs}$')

    # Obtain PDMF factor.
    h2, x = np.histogram(m_act_obs)
    h2_log = np.log10(h2)
    msk = h2_log == -np.inf
    alpha_p, _ = np.polyfit(np.log10(x[:-1][~msk]), h2_log[~msk], 1)
    m, b = np.polyfit(x[1:], h2_log[~msk], 1)
    plt.plot(
        x, 10**np.poly1d((m, b))(x), c='orange',
        label=r"$\alpha_{{p}}={:.3f}$".format(alpha_p))
    plt.scatter(x[1:], h2, c='orange', label=r'$M_{p}^{obs}$')

    y_max = max(h1.max(), h2.max())
    plt.xlabel(r"$m\,(M_{\odot})$")
    plt.ylabel(r"$\log[\xi(m)]$")
    plt.yscale('log')
    plt.ylim(0., y_max + .35 * y_max)
    plt.legend()
    plt.show()

    import pdb; pdb.set_trace()  # breakpoint b2c3ec89 //

