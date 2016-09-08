
import time
import gc  # Garbage collector.
#
from inp import names_paths
from inp import get_data_semi
from inp import get_data
from structure import trim_frame
from structure import histo_2d
from structure import center
from structure import radial_dens_prof
from structure import field_density
from structure import radius
from structure import king_profile
from errors import err_accpt_rejct
from structure import stars_in_out_cl_reg
from structure import cluster_area
from phot_analysis import members_number
from phot_analysis import contamination_index
from structure import field_regions
from phot_analysis import luminosity_func
from phot_analysis import integrated_mag
from phot_analysis import kde_pvalue
from decont_algors import bayesian_da
from phot_analysis import members_N_compare
from decont_algors import cl_region_clean
from out import cluster_members_file
from best_fit import synth_cl_err
from best_fit import best_fit_synth_cl
from out import synth_cl_file
from out import create_out_data_file
from out import add_data_output
from out import top_tiers
from out import make_A_plot
from out import make_B_plot
from out import make_C_plot
from out import make_plots
from out import done_move


def main(cl_file, pd):
    '''
    Container that holds the calls to all the moduli and functions.
    '''

    # Start timing this loop.
    start = time.time()

    # Get file names (n) and paths (p) dictionary (d).
    npd = names_paths.main(cl_file, **pd)

    # Get data from semi-data input file. Add to dictionary.
    pd = get_data_semi.main(pd, **npd)

    # Get cluster's data from file, as dictionary.
    cld = get_data.main(npd, **pd)

    # If Manual mode is set, display frame and ask if it should be trimmed.
    cld = trim_frame.main(cld, **pd)

    # Obtain 2D coordinates histogram for the observed frame.
    # Return cluster's parameters dictionary 'clp'.
    clp = histo_2d.main(pd, **cld)

    # Get cluster's center coordinates and errors.
    clp = center.main(cld, clp, **pd)

    # Get density profile
    clp = radial_dens_prof.main(clp)

    # Get field density value in stars/<area unit>.
    clp = field_density.main(clp, **pd)

    # Get cluster radius
    clp = radius.main(cld, clp, **pd)
    # clust_rad = radius_params[0]

    # Get King profiles based on the density profiles.
    clp = king_profile.main(clp, **pd)

    # Accept and reject stars based on their errors.
    clp, pd = err_accpt_rejct.main(cld, clp, pd)

    # Get stars in and out of cluster's radius.
    clp = stars_in_out_cl_reg.main(clp)

    # Get cluster's area.
    clp = cluster_area.main(clp, **cld)

    # Get approximate number of cluster's members.
    clp = members_number.main(clp)

    # Get contamination index.
    clp = contamination_index.main(clp)

    # Field regions around the cluster's center.
    clp = field_regions.main(clp, **pd)

    # Make A block plots.
    # make_A_plot.main(npd, cld, pd, **clp)

    # Get luminosity function and completeness level for each magnitude bin.
    clp = luminosity_func.main(clp, **cld)

    # Calculate integrated magnitude.
    clp = integrated_mag.main(clp, **pd)

    # Get physical cluster probability based on p_values distribution.
    clp = kde_pvalue.main(clp, **pd)

    # Make B block plots.
    # make_B_plot.main(npd, cld, pd, **clp)

    # Apply decontamination algorithm.
    clp = bayesian_da.main(clp, npd, **pd)

    # Obtain members parameter.
    clp = members_N_compare.main(clp)

    # Reduce number of stars in cluster according to a lower membership
    # probability or magnitude limit.
    clp = cl_region_clean.main(clp, **pd)

    # Make C block plots.
    make_C_plot.main(npd, cld, pd, **clp)

    # Create data file with membership probabilities.
    cluster_members_file.main(clp, **npd)

    # Obtain exponential error function parameters to use by the
    # synthetic cluster creation function.
    clp = synth_cl_err.main(cld, clp, **pd)

    # Obtain best fitting parameters for cluster.
    clp = best_fit_synth_cl.main(clp, **pd)

    # Create output synthetic cluster file if one was found
    synth_cl_file.main(clp, npd, **pd)

    # Create output data file in /output dir if it doesn't exist.
    create_out_data_file.main(npd)

    # Add cluster data and flags to output file
    add_data_output.main(npd, pd, **clp)

    # Output top tiers models if best fit parameters were obtained.
    top_tiers.main(npd, cld, pd, **clp)

    # Make plots
    make_plots.main(npd, cld, pd, **clp)

    # Move file to 'done' dir (if flag is set).
    done_move.main(pd, **npd)

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print('End of analysis for {} in {:.0f}m {:.0f}s.\n'.format(
        npd['clust_name'], m, s))

    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()
