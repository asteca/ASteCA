import numpy as np
from scipy.spatial import KDTree

from functools import wraps
from time import time
def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print('func:%r took: %2.4f sec' % (f.__name__, te-ts))
        return result
    return wrap


@timing
def newfastMP(
    ra_c, dec_c, pmra_c, pmde_c, plx_c, ra_v, dec_v, pmra_v, pmde_v, plx_v,
    e_pmra_v, e_pmde_v, e_plx_v, N_cluster, rng, N_runs
):
    """
    """

    xy = np.array([ra_v, dec_v]).T
    d_pm = np.sqrt((pmra_c - pmra_v)**2 + (pmde_c - pmde_v)**2)
    e_d_pm = np.sqrt(
        ((pmra_v*e_pmra_v)**2 + (pmde_v*e_pmde_v)**2) / (pmra_v**2 + pmde_v**2)
    )
    d_plx = abs(plx_c - plx_v)

    N_tot = len(ra_v)
    N_break = 50
    probs_all = np.zeros(N_tot)
    prob_old = np.zeros(N_tot)
    for r in range(N_runs):

        # Reset step
        N_step = 25

        d_idxs_s = [np.argsort(rng.normal(d_pm, e_d_pm))]
        d_idxs_s += [np.argsort(rng.normal(d_plx, e_plx_v))]
        d_idxs_s = np.array(d_idxs_s)

        while True:
            i_select = list(set.intersection(*map(set, d_idxs_s[:, :N_step])))
            if len(i_select) < N_cluster:
                N_step += N_step
            else:
                break
        i_select = np.array(i_select)

        if len(i_select) > N_cluster:
            mdists = mean_dists_kdtree(xy[i_select])
            idxs = np.argsort(mdists)[:N_cluster]
            i_select = i_select[idxs]

        probs_all[i_select] += 1
        probs = probs_all / (r+1)

        msk = probs > 0.5
        # Check that all P>0.5 probabilities converged to 1%
        if (abs(prob_old[msk] - probs[msk]) < 0.01).all() and r > N_break:
            break
        else:
            prob_old = np.array(probs)

    if r < N_runs:
        print(f"Convergence reached at {r} runs")
    else:
        print(f"Maximum number of runs reached: {N_runs}")

    return probs


def mean_dists_kdtree(coordinates, N_neigh=5):
    tree = KDTree(coordinates)
    # Query the tree for the N_neigh nearest neighbors (including the point itself)
    distances, closest_indices = tree.query(coordinates, k=N_neigh + 1)
    # Mean distance to the NN_dd neighbors. Remove the distance to itself
    return distances[:, 1:].mean(1)
