
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt


def main(clp):  # m_ini_obs
    """
    TODO
    1. use the LF maximum magnitude cut, otherwise the incompleteness of low
    mass stars impacts on the slope estimation.
    2. perhaps fit two distinct lines: m<1 and m>1? This is where most IMFs
    seem to make a break in the slope.
    3. discard outliers with distances too large
    4. Use the Bayesian linear fit instead of polyfit()
    """

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
    def func(maxm, bins):
        m_ini_obs, m_act_obs = m_ini[minid], m_act[minid]

        # Filter out stars with large photometric distances or below the LF
        # 'maximum magnitude'.
        d_med, d_std = np.median(mindist), np.std(mindist)
        msk = (mindist < d_med * 3. * d_std) & (phot_obs.T[0] < maxm)
        print(len(mindist), sum(msk))
        m_ini_obs, m_act_obs = m_ini_obs[msk], m_act_obs[msk]

        # Obtain IMF factor.
        h1, x = np.histogram(m_ini_obs, bins=bins)
        # Mask -inf values that can appear when a bin contains 0 elements.
        h1_log = np.log10(h1)
        msk = h1_log == -np.inf
        x_log = np.log10(x[:-1][~msk])
        alpha_i, b = np.polyfit(x_log, h1_log[~msk], 1)
        # m, b = np.polyfit(x[:-1][~msk], h1_log[~msk], 1)
        plt.plot(
            x, 10**(b) * x**(alpha_i),
            c='#7BC5D4', label=r"$\alpha_{{i}}={:.3f}$".format(-alpha_i))
        plt.scatter(x[1:], h1, c='#7BC5D4', label=r'$M_{i}^{obs}$')

        # Obtain PDMF factor.
        h2, x = np.histogram(m_act_obs, bins=bins)
        h2_log = np.log10(h2)
        msk = h2_log == -np.inf
        x_log = np.log10(x[:-1][~msk])
        alpha_p, b = np.polyfit(x_log, h2_log[~msk], 1)
        # m, b = np.polyfit(x[1:], h2_log[~msk], 1)
        plt.plot(
            x, 10**(b) * x**(alpha_p),
            c='orange', label=r"$\alpha_{{p}}={:.3f}$".format(-alpha_p))
        plt.scatter(x[1:], h2, c='orange', label=r'$M_{p}^{obs}$')

        # y_max = max(h1.max(), h2.max())
        # y_max = h1_log[~msk].max()
        plt.xlabel(r"$m\,[M_{\odot}]$")
        plt.ylabel(r"$\xi(m) \Delta m$")
        plt.xscale('log')
        plt.yscale('log')
        plt.minorticks_on()
        # plt.ylim(0., y_max + .35 * y_max)
        plt.legend()
        plt.show()
    func(18.5, 20)

    import pdb; pdb.set_trace()  # breakpoint b2c3ec89 //


if __name__ == '__main__':
    m_ini_obs = np.array([1.28296444, 1.28296444, 1.2789754 , 1.27510215, 1.2789754 ,
           1.26973795, 1.264633  , 1.26258727, 1.25953579, 1.25953579,
           1.25953579, 1.25423097, 1.25423097, 1.25423097, 1.25423097,
           1.25423097, 1.25030026, 1.25030026, 1.25030026, 1.25030026,
           1.25030026, 1.25030026, 1.25030026, 1.25030026, 1.25030026,
           1.25030026, 1.25030026, 1.22025002, 1.25030026, 1.24707987,
           1.24707987, 1.24707987, 1.24707987, 1.22025002, 1.22025002,
           1.19865899, 1.22025002, 1.22025002, 1.22426759, 1.22025002,
           1.22025002, 1.22025002, 1.22025002, 1.22025002, 1.20524902,
           1.19865899, 1.22426759, 1.22426759, 1.22426759, 1.20524902,
           1.20524902, 1.1881938 , 1.1881938 , 1.1881938 , 1.17673556,
           1.17795167, 1.1881938 , 1.23558953, 1.1881938 , 1.1881938 ,
           1.1881938 , 1.1881938 , 1.17184569, 1.23558953, 1.17184569,
           1.17184569, 1.17184569, 1.23596523, 1.16632464, 1.16632464,
           1.13775778, 1.15363711, 1.14244811, 1.13775778, 1.13775778,
           1.12121604, 1.12121604, 1.12121604, 1.12121604, 1.10536567,
           1.0649941 , 1.08368668, 1.0909959 , 1.0909959 , 1.0909959 ,
           1.08368668, 1.0909959 , 1.0909959 , 1.08368668, 1.08368668,
           1.08003206, 1.08003206, 1.0649941 , 1.06878804, 1.0649941 ,
           1.06116414, 1.06116414, 1.06116414, 1.06116414, 1.0649941 ,
           1.06116414, 1.06116414, 1.05383514, 1.0649941 , 1.05383514,
           1.05267715, 1.0649941 , 1.0649941 , 1.05267715, 1.0479548 ,
           1.0479548 , 1.0479548 , 1.03877277, 1.03877277, 1.02305627,
           1.01895533, 1.01895533, 1.00419911, 1.00082745, 1.00251328,
           1.00082745, 1.00082745, 0.99333663, 0.99333663, 0.99333663,
           0.99333663, 0.9841773 , 0.9841773 , 0.9841773 , 0.96638566,
           0.96638566, 0.96638566, 0.96638566, 0.96638566, 0.96638566,
           0.95748984, 0.95099106, 0.94846753, 0.95099106, 0.93328222,
           0.94846753, 0.94087487, 0.94087487, 0.93328222, 0.94087487,
           0.94087487, 0.94087487, 0.92568956, 0.93328222, 0.93328222,
           0.91162141, 0.91162141, 0.91162141, 0.91162141, 0.89726957,
           0.89726957, 0.8807354 , 0.8807354 , 0.89726957, 0.89726957,
           0.86597362, 0.8807354 , 0.8807354 , 0.86597362, 0.86597362,
           0.86597362, 0.86597362, 0.86597362, 0.86597362, 0.83081342,
           0.83081342, 0.83081342, 0.83081342, 0.83081342, 0.83081342,
           0.83081342, 0.83081342, 0.80278639, 0.80278639, 0.83081342,
           0.77904254, 0.77904254, 0.77904254])
    main(m_ini_obs)
