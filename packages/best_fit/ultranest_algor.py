
import ultranest
from .bf_common import getSynthClust
from . import likelihood


def main(lkl_method, obs_clust, varIdxs, ranges, synthcl_args):
    """
    Did not get good results at least using the default options
    """

    def getLikelihood(params):
        model = params
        synth_clust = getSynthClust(model, True, synthcl_args)[0]
        lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

        return lkl

    def prior_transform(cube, varIdxs=varIdxs, ranges=ranges):
        # the argument, cube, consists of values from 0 to 1
        # we have to convert them to physical scales
        params = cube.copy()

        lo, hi = ranges[varIdxs].T
        # transform location parameter: uniform prior
        params = cube * (hi - lo) + lo

        return params

    param_names = ['z', 'a', 'Av', 'dr', 'dm']

    sampler = ultranest.ReactiveNestedSampler(
        param_names, getLikelihood, prior_transform)
    result = sampler.run(max_iters=3000)
    breakpoint()

    sampler.print_results()
    import matplotlib.pyplot as plt
    sampler.plot_trace()
    plt.show()
    sampler.plot_corner()
    plt.show()

    # from ultranest.plot import cornerplot
    # cornerplot(result)
    # from ultranest.plot import traceplot
    # traceplot(result)  # KeyError: 'logvol' ??
    # plt.show()
