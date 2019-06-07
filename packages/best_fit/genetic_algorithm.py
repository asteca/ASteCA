
import time as t
import textwrap
import random
import numpy as np
from .bf_common import varPars, rangeCheck, synthClust, random_population
# from ..synth_clust import synth_cluster
from . import likelihood

#############################################################
# # Timer function: http://stackoverflow.com/a/21860100/1391441
# from contextlib import contextmanager
# import time


# @contextmanager
# def timeblock(label):
#     start = time.clock()
#     try:
#         yield
#     finally:
#         end = time.clock()
#         print ('{} elapsed: {}'.format(label, end - start))
#############################################################


def main(
    available_secs, ran_pop, flag_print_perc, N_popl, N_gener, max_mag_syn,
    obs_clust, ext_coefs, st_dist_mass, N_fc, cmpl_rnd, err_rnd, e_max,
    err_lst, completeness, lkl_method, fundam_params, theor_tracks, R_V,
        fit_diff, cross_prob, cross_sel, mut_prob, N_el, N_ei, N_es, **kwargs):
    '''
    Genetic algorithm. Finds the best fit model-observation by minimizing the
    likelihood function.
    '''

    # Check if N_popl is odd. If it is sum 1 to avoid conflict if cross_sel
    # '2P' was selected.
    N_popl += N_popl % 2

    # Get number of binary digits to use.
    n_bin, p_delta, p_mins = num_binary_digits(fundam_params)

    # Fitness.
    # Rank-based breeding probability. Independent of the fitness values,
    # only depends on the total number of chromosomes N_popl and the fitness
    # differential fit_diff.
    fitness = [1. / N_popl + fit_diff * (N_popl + 1. - 2. * (i + 1.)) /
               (N_popl * (N_popl + 1.)) for i in range(N_popl)]

    # For plotting purposes.
    # Stores parameters of the solutions already processed and the likelihoods
    # obtained.
    models_GA = np.zeros((N_gener * N_popl, 6))
    lkls_GA = np.zeros(N_gener * N_popl)
    lkl_best = np.full(N_gener, np.inf)
    lkl_mean = np.zeros(N_gener)

    varIdxs, ndim, ranges = varPars(fundam_params)

    # Evaluate initial random solutions in the objective function.
    generation, lkl = evaluation(
        lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
        obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
        err_rnd, varIdxs, ranges, ran_pop)

    # Store best solution(s) for passing along in the 'Elitism' block.
    best_sol = generation[:N_el]

    # Stores indexes where a new best solution was found.
    new_bs_indx = []
    # Initiate counters.
    best_sol_count, ext_imm_count = 0, 0
    # Minimum likelihood found after each application of the
    # Extinction/Immigration operator.
    lkl_best_old, lkl_ei = lkl[0], np.inf

    # Print percentage done.
    step, elapsed_in, start_in = 0, 0., t.time()
    milestones = list(range(10, 101, 10))
    # Begin processing the populations up to N_gen generations.
    while elapsed_in < available_secs:

        # *** Selection/Reproduction ***
        # with timeblock("Selec/Repo"):
        # Select chromosomes for breeding from the current generation of
        # solutions according to breed_prob to generate the intermediate
        # population.
        int_popul = selection(generation, fitness)

        # Encode intermediate population's solutions into binary
        # chromosomes.
        chromosomes = encode(n_bin, p_delta, p_mins, int_popul)

        # *** Breeding ***
        # with timeblock("Breeding"):
        # Pair chromosomes by randomly shuffling them.
        random.shuffle(chromosomes)

        # Apply crossover operation on each subsequent pair of chromosomes
        # with a cross_prob probability (crossover probability)
        cross_chrom = crossover(chromosomes, cross_prob, cross_sel)

        # Apply mutation operation on random genes for every chromosome.
        mut_chrom = mutation(cross_chrom, mut_prob)

        # *** Evaluation ***
        # Decode the chromosomes into solutions to form the new generation.
        # with timeblock("Decode"):
        p_lst_d = decode(fundam_params, n_bin, p_delta, p_mins, mut_chrom)

        # Elitism: make sure that the best N_el solutions from the previous
        # generation are passed unchanged into this next generation.
        # with timeblock("Elitism"):
        p_lst_e = elitism(best_sol, p_lst_d)

        # Evaluate each new solution in the objective function and sort
        # according to the best solutions found.
        # with timeblock("Evaluation"):
        generation, lkl = evaluation(
            lkl_method, e_max, err_lst, completeness, max_mag_syn,
            fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
            st_dist_mass, N_fc, cmpl_rnd, err_rnd, varIdxs, ranges, p_lst_e)

        # *** Extinction/Immigration ***
        # If the best solution has remained unchanged for N_ei
        # generations, remove all chromosomes but the best ones (extinction)
        # and fill with random new solutions (immigration).

        # Check if new best solution is better than the previous one.
        if lkl_best_old <= lkl[0]:
            # Increase counter.
            best_sol_count += 1

            # Check how many times the best_sol has remained unchanged.
            # If the number equals N_ei, apply Extinction/Immigration operator.
            if best_sol_count == N_ei:

                # *** Exit switch ***
                # If N_es runs of the Ext/Imm operator have been applied with
                # no changes to the best solution, apply the exit switch.
                if np.allclose(lkl_ei, lkl_best_old, .01):
                    # Increase Ext/Imm operator counter.
                    ext_imm_count += 1
                    if ext_imm_count == N_es:
                        if flag_print_perc:
                            print("  GA exit switch applied.")
                        elapsed_in = t.time() - start_in
                        break
                else:
                    # Update best solution.
                    lkl_ei = lkl_best_old
                    # Reset counter.
                    ext_imm_count = 0

                # Apply Extinction/Immigration operator.
                generation = ext_imm(best_sol, fundam_params, N_popl)
                # Reset best solution counter.
                best_sol_count = 0

        else:
            lkl_best_old = lkl[0]
            # For plotting purposes. Save index where a new best solution
            # was found.
            new_bs_indx.append(step)
            # Update best solution for passing along in the 'Elitism' block.
            best_sol = generation[:N_el]
            # Reset counter.
            best_sol_count = 0

        if flag_print_perc:

            # For plotting purposes.
            models_GA[N_popl * step:N_popl * (step + 1)] = generation
            lkls_GA[N_popl * step:N_popl * (step + 1)] = lkl
            lkl_best[step] = lkl[0]
            # Discard large values associated with empty arrays from mean.
            lkl_mean[step] = np.mean(
                np.asarray(lkl)[np.asarray(lkl) < 9.9e08])

            # Time used to check how fast the sampler is advancing.
            elapsed_in = t.time() - start_in
            percentage_complete = (100. * elapsed_in / available_secs)
            if len(milestones) > 0 and percentage_complete >= milestones[0]:
                m, s = divmod(max(1., available_secs - elapsed_in), 60)
                h, m = divmod(m, 60)
                # m += s / 60.
                print("{:>3}% LP={:.1f} ({:.5f}, {:.3f}, {:.3f}, "
                      "{:.2f}, {:g}, {:.2f})".format(
                          milestones[0], lkl[0], *generation[0]) +
                      " [{:.0f} m/s | {:.0f}h{:.0f}m]".format(
                          (N_popl * step) / elapsed_in, h, m))
                # Remove that milestone from the list.
                milestones = milestones[1:]

            # If the hardcoded maximum number of generations is reached.
            if step == N_gener:
                print(" WARNING: maximum allowed number of " +
                      "generations ({}) reached.".format(N_gener))
                break

        else:
            if step == N_gener:
                break
        step += 1

    # If this is a bootstrap run, return the best model found only.
    if not flag_print_perc:
        return generation[0]

    # For plotting
    # In case not all the generations were processed.
    models_GA = models_GA[:(step - 1) * N_popl]
    lkls_GA = lkls_GA[:(step - 1) * N_popl]
    lkl_best = lkl_best[:(step - 1)]
    lkl_mean = lkl_mean[:(step - 1)]

    # Remove possible large Lkl values.
    msk = lkls_GA < 9.9e08
    models_GA, lkls_GA = models_GA[msk], lkls_GA[msk]
    # Sort and trim to N_max models to plot
    idx_sort = np.argsort(lkls_GA)
    N_max = 10000
    models_GA = models_GA[idx_sort][:N_max].T
    lkls_GA = lkls_GA[idx_sort][:N_max]

    # Prevent very small 'binary fraction' values from disrupting the synthetic
    # cluster plotting.
    map_sol = generation[0]
    map_sol[5] = 0. if map_sol[5] < 0.01 else map_sol[5]

    # Store the ML solution as 'map_sol' for consistency.
    isoch_fit_params = {
        'OF_final_generation': generation, 'OF_elapsed': elapsed_in,
        'map_sol': map_sol, 'lkl_best': lkl_best, 'lkl_mean': lkl_mean,
        'OF_steps': step + 1, 'OF_models': (step + 1) * N_popl,
        'new_bs_indx': new_bs_indx, 'models_GA': models_GA, 'lkls_GA': lkls_GA}

    return isoch_fit_params


def num_binary_digits(fundam_params):
    '''
    Store parameters ranges and calculate the minimum number of binary digits
    needed to encode the solutions.
    '''

    p_mins, p_delta = [], []
    for param in fundam_params:
        p_mins.append(min(param))  # Used by the encode operator.
        # If delta is zero it means a single value is being used. Set delta
        # to a small value to avoid issues with Elitism operator.
        p_delta.append(max(max(param) - min(param), 1e-10))

    p_interv = np.array([len(_) for _ in fundam_params])

    # Number of binary digits used to create the chromosomes.
    # The max function prevents an error when all parameter ranges are set to
    # a unique value in which case the np.log is a negative float.
    # We add 10 since it's very cheap and adds a lot of accuracy.
    n_bin = max(int(np.log(max(p_interv)) / np.log(2)) + 1, 1) + 10

    return n_bin, p_delta, p_mins


def encode(n_bin, p_delta, p_mins, int_popul):
    '''
    Encode the solutions into binary string chromosomes to be bred.
    '''

    chromosomes = []
    for sol in int_popul:
        # Convert floats to binary strings.
        p_binar = []
        for i, p_del in enumerate(p_delta):
            p_binar.append(str(bin(int(
                ((sol[i] - p_mins[i]) / p_del) * (2 ** n_bin - 1))))[2:].zfill(
                    n_bin))

        # Combine binary strings to generate a single chromosome.
        chromosomes.append(''.join(p_binar))

    return chromosomes


def decode(fundam_params, n_bin, p_delta, p_mins, mut_chrom):
    '''
    Decode the chromosomes into its real values to be evaluated by the
    objective function.
    '''

    p_lst = []
    for chrom in mut_chrom:
        # Split chromosome string.
        chrom_split = textwrap.wrap(chrom, n_bin)

        # Map integers to the real parameter values.
        sol = []
        for i, p_del in enumerate(p_delta):
            # Convert binary into an integer.
            b2i = int(chrom_split[i], 2)
            # Map integer to the real parameter value.
            sol.append(p_mins[i] + (b2i * p_del / (2 ** n_bin - 1)))
        p_lst.append(sol)

    return p_lst


def evaluation(
    lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
    obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
        err_rnd, varIdxs, ranges, p_lst):
    '''
    Evaluate each model in the objective function to obtain the fitness of
    each one.
    '''

    # lkl_lst, generation_list = [], []
    lkl_arr = np.full(len(p_lst), np.inf)
    # Process each model selected.
    for i, model in enumerate(p_lst):

        rangeFlag = rangeCheck(np.array(model)[varIdxs], ranges, varIdxs)
        if rangeFlag:

            # Generate synthetic cluster.
            synth_clust = synthClust(
                fundam_params, varIdxs, model, (
                    theor_tracks, e_max, err_lst, completeness, max_mag_syn,
                    st_dist_mass, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd))

            # Call likelihood function for this model.
            # with timeblock(" Likelihood"):
            lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

            # DEPRECATED May 2019, there's no need now that the IMF is
            # sampled once outside of the best fit process.
            # # Check if this model was already processed. Without this check here,
            # # the extinction/immigration operator becomes useless, since the best
            # # solution's likelihood value can vary slightly due to the re-sampling
            # # of the mass distribution.
            # try:
            #     # Get old likelihood value for this model.
            #     lkl_old = model_done[1][model_done[0].index(model)]
            #     # Compare with new value. If old one is smaller, pass old one.
            #     # This keeps the likelihood always trending downwards.
            #     if lkl_old < lkl:
            #         lkl = lkl_old
            # except ValueError:
            #     # Model was not processed before.
            #     pass

            lkl_arr[i] = lkl

    # Sort according to the likelihood list. This puts the best model (ie:
    # the one with the minimum likelihood value) first.
    # with timeblock(" sort"):
    generation = [x for y, x in sorted(list(zip(lkl_arr, p_lst)))]

    # Sort list in place putting the likelihood minimum value first.
    lkl_arr.sort()

    return generation, lkl_arr


def selection(generation, breed_prob):
    '''
    Select random chromosomes from the chromosome list passed according to
    the breeding probability given by their fitness.
    '''
    select_chrom = []
    # Draw len(generation) random numbers uniformly distributed between [0,1)
    ran_lst = np.random.uniform(0, sum(breed_prob), len(generation))
    # For each of these numbers, obtain the corresponding likelihood from the
    # breed_prob CDF.
    gen_breed = list(zip(*[generation, breed_prob]))
    for r in ran_lst:
        s = 0.
        for sol, num in gen_breed:
            s += num
            if s >= r:
                select_chrom.append(sol)
                break

    return select_chrom


def chunker(seq, size):
    '''
    Helper function for 'crossover'. Returns elements in list in "chunks"
    rather than one at a time.
    '''
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def crossover(chromosomes, cross_prob, cross_sel):
    '''
    Applies the crossover operator over each chromosome.
    '''
    cross_chrom = []
    # Take two chromosomes at a time.
    for chrom_pair in chunker(chromosomes, 2):
        r = random.random()
        if r <= cross_prob:

            if cross_sel == '1P':
                # Select one random crossover point.
                cp = random.randint(0, len(chrom_pair[0]))
                # Apply crossover on these two chromosomes.
                cross_chrom.append(chrom_pair[0][:cp] + chrom_pair[1][cp:])
                cross_chrom.append(chrom_pair[1][:cp] + chrom_pair[0][cp:])
            elif cross_sel == '2P':
                # Select two random crossover points.
                cp0, cp1 = np.sort(random.sample(range(0, len(chrom_pair[0])),
                                                 2))
                # Apply crossover on these two chromosomes.
                cross_chrom.append(chrom_pair[0][:cp0] +
                                   chrom_pair[1][cp0:cp1] +
                                   chrom_pair[0][cp1:])
                cross_chrom.append(chrom_pair[1][:cp0] +
                                   chrom_pair[0][cp0:cp1] +
                                   chrom_pair[1][cp1:])
        else:
            # Skip crossover operation.
            cross_chrom.append(chrom_pair[0])
            cross_chrom.append(chrom_pair[1])

    return cross_chrom


def mutation(cross_chrom, mut_prob):
    '''
    Applies the mutation operator over random genes in each chromosome.
    '''
    # For each chromosome flip their genes according to the probability
    # mut_prob.
    for i, elem in enumerate(cross_chrom):
        cross_chrom[i] = ''.join(char if random.random() > mut_prob else
                                 str(1 - int(char)) for char in elem)

    return cross_chrom


def elitism(best_sol, p_lst):
    '''
    Pass the best N_el solutions unchanged to the next generation.
    '''

    # Append the N_el best solutions to the beginning of the list, pushing
    # out N_el solutions located on the end of the list.
    p_lst_r = best_sol + p_lst[:-len(best_sol)]

    return p_lst_r


def ext_imm(best_sol, fundam_params, N_popl):
    '''
    Append a new random population to the best solution so far.
    '''

    # Generate (N_popl-N_el) random solutions.
    n_ran = N_popl - len(best_sol)
    p_lst_r = random_population(fundam_params, (0, 1, 2, 3, 4, 5), n_ran)

    # Append immigrant random population to the best solution.
    generation_ei = best_sol + p_lst_r.tolist()

    return generation_ei
