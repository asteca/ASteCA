
import random
import numpy as np
from ..synth_clust import synth_cluster
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


def encode(n_bin, p_delta, p_mins, int_popul):
    '''
    Encode the solutions into binary string chromosomes to be bred.
    '''

    chromosomes = []
    for sol in int_popul:
        # Convert floats to binary strings.
        p_binar = []
        for i, p_del in enumerate(p_delta):
            p_binar.append(str(bin(int(((sol[i] - p_mins[i]) / p_del) *
                           (2 ** n_bin - 1))))[2:].zfill(n_bin))

        # Combine binary strings to generate chromosome.
        chrom = ''.join(p_binar)
        chromosomes.append(chrom)

    return chromosomes


def chunker(seq, size):
    '''
    Helper function for 'crossover'. Returns elements in list in "chunks"
    rather than one at a time.
    '''
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


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


def decode(fundam_params, n_bin, p_delta, p_mins, mut_chrom):
    '''
    Decode the chromosomes into its real values to be evaluated by the
    objective function.
    '''

    # Initialize empty list with as many sub-lists as parameters.
    p_lst = [[] for _ in range(len(fundam_params))]

    for chrom in mut_chrom:
        # Split chromosome string.
        chrom_split = [chrom[i:i + n_bin] for i in xrange(0, len(chrom),
                       n_bin)]

        # Convert each binary into an integer for all parameters.
        b2i = [int(i, 2) for i in chrom_split]

        # Map integers to the real parameter values.
        for i, p_del in enumerate(p_delta):
            # Map integer to the real parameter value.
            p_r = p_mins[i] + (b2i[i] * p_del / (2 ** n_bin - 1))
            # Find the closest value in the parameters list.
            p = min(fundam_params[i], key=lambda x: abs(x - p_r))
            # Store this last value.
            p_lst[i].append(p)

    return p_lst


def elitism(best_sol, p_lst):
    '''
    Pass the best N_el solutions unchanged to the next generation.
    '''

    # Append the N_el best solutions to the beginning of the list, pushing
    # out N_el solutions located on the end (right) of the list.
    p_lst_r = []
    for i, pars in enumerate(p_lst):
        p_lst_r.append(list(zip(*best_sol)[i]) + pars[:-len(best_sol)])

    return p_lst_r


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
    gen_breed = zip(*[generation, breed_prob])
    for r in ran_lst:
        s = 0.
        for sol, num in gen_breed:
            s += num
            if s >= r:
                select_chrom.append(sol)
                break

    return select_chrom


def evaluation(lkl_method, e_max, err_lst, completeness, max_mag_syn,
               fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
               st_dist_mass, N_fc, cmpl_rnd, err_rnd, p_lst, model_done):
    '''
    Evaluate each model in the objective function to obtain the fitness of
    each one.
    '''

    lkl_lst, generation_list = [], []
    # Process each model selected.
    for model in zip(*p_lst):

        # Metallicity and age indexes to identify isochrone.
        m_i = fundam_params[0].index(model[0])
        a_i = fundam_params[1].index(model[1])
        isochrone = theor_tracks[m_i][a_i]

        # Generate synthetic cluster.
        synth_clust = synth_cluster.main(
            e_max, err_lst, completeness, max_mag_syn, st_dist_mass, isochrone,
            R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, model)

        # Call likelihood function for this model.
        # with timeblock(" Likelihood"):
        lkl = likelihood.main(lkl_method, synth_clust, obs_clust)

        # Check if this model was already processed. Without this check here,
        # the extinction/immigration operator becomes useless, since the best
        # solution's likelihood value can vary slightly due to the re-sampling
        # of the mass distribution.
        # if model in model_done[0]:
        # with timeblock(" Compare"):
        try:
            # Get old likelihood value for this model.
            lkl_old = model_done[1][model_done[0].index(model)]
            # Compare with new value. If old one is smaller, pass old one.
            # This keeps the likelihood always trending downwards.
            if lkl_old < lkl:
                lkl = lkl_old
        except ValueError:
            # Model was not processed before.
            pass

        # Append data to the lists that will be erased with each call
        # to this function.
        generation_list.append(model)
        lkl_lst.append(lkl)

    # Sort according to the likelihood list. This puts the best model (ie:
    # the one with the minimum likelihood value) first.
    # with timeblock(" sort"):
    generation = [x for y, x in sorted(zip(lkl_lst, generation_list))]
    # Sort list in place putting the likelihood minimum value first.
    lkl_lst.sort()

    # Append data identifying the isochrone and the obtained
    # likelihood value to this *persistent* list.
    # with timeblock(" Append"):
    model_done[0] = generation + model_done[0]
    model_done[1] = lkl_lst + model_done[1]

    return generation, lkl_lst, model_done


def random_population(fundam_params, n_ran):
    '''
    Generate a random set of parameter values to use as a random population.
    '''
    # Pick n_ran initial random solutions from each list storing all the
    # possible parameters values. These lists store real values.
    p_lst = []
    for param in fundam_params:
        p_lst.append([random.choice(param) for _ in range(n_ran)])

    return p_lst


def ext_imm(best_sol, fundam_params, N_pop):
    '''
    Append a new random population to the best solution so far.
    '''

    # Generate (N_pop-N_el) random solutions.
    n_ran = N_pop - len(best_sol)
    p_lst_r = random_population(fundam_params, n_ran)

    # Append immigrant random population to the best solution.
    generation_ei = best_sol + zip(*p_lst_r)

    return generation_ei


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


def main(lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
         obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
         err_rnd, N_pop, N_gen, fit_diff, cross_prob, cross_sel, mut_prob,
         N_el, N_ei, N_es, flag_print_perc):
    '''
    Genetic algorithm. Finds the best fit model-observation.
    '''

    # Check if N_pop is odd. If it is sum 1 to avoid conflict if cross_sel
    # '2P' was selected.
    N_pop += N_pop % 2

    # Get number of binary digits to use.
    n_bin, p_delta, p_mins = num_binary_digits(fundam_params)

    # Fitness.
    # Rank-based breeding probability. Independent of the fitness values,
    # only depends on the total number of chromosomes N_pop and the fitness
    # differential fit_diff.
    fitness = [1. / N_pop + fit_diff * (N_pop + 1. - 2. * (i + 1.)) /
               (N_pop * (N_pop + 1.)) for i in range(N_pop)]

    # *** Initial random population evaluation ***
    p_lst_r = random_population(fundam_params, N_pop)

    # Stores parameters of the solutions already processed and the likelihoods
    # obtained.
    model_done = [[], []]

    # Evaluate initial random solutions in the objective function.
    generation, lkl, model_done = evaluation(
        lkl_method, e_max, err_lst, completeness, max_mag_syn, fundam_params,
        obs_clust, theor_tracks, R_V, ext_coefs, st_dist_mass, N_fc, cmpl_rnd,
        err_rnd, p_lst_r, model_done)

    # Store best solution for passing along in the 'Elitism' block.
    best_sol = generation[:N_el]

    # For plotting purposes.
    lkl_old = [[], []]
    # Stores indexes where a new best solution was found.
    new_bs_indx = []

    # Initiate counters.
    best_sol_count, ext_imm_count = 0, 0
    # Initiate empty list. Stores the best solution found after each
    # application of the Extinction/Immigration operator.
    best_sol_ei = []

    # Print percentage done.
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # Begin processing the populations up to N_gen generations.
    for i in range(N_gen):

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
        generation, lkl, model_done = evaluation(
            lkl_method, e_max, err_lst, completeness, max_mag_syn,
            fundam_params, obs_clust, theor_tracks, R_V, ext_coefs,
            st_dist_mass, N_fc, cmpl_rnd, err_rnd, p_lst_e, model_done)

        # *** Extinction/Immigration ***
        # If the best solution has remained unchanged for N_ei
        # generations, remove all chromosomes but the best ones (extinction)
        # and fill with random new solutions (immigration).

        # Check if new best solution is equal to the previous one.
        if generation[0] == best_sol[0]:
            # Increase counter.
            best_sol_count += 1

            # Check how many times the best_sol has remained unchanged.
            # If the number equals N_ei, apply Extinction/Immigration operator.
            if best_sol_count == N_ei:

                # *** Exit switch ***
                # If N_es runs of the Ext/Imm operator have been applied with
                # no changes to the best solution, apply the exit switch,
                # ie: exit the GA.
                if best_sol[0] == best_sol_ei:
                    # Increase Ext/Imm operator counter.
                    ext_imm_count += 1
                    if ext_imm_count == N_es:
                        # Exit generations loop.
                        if flag_print_perc:
                            print("  GA exit switch applied.")
                        break
                else:
                    # Update best solution.
                    best_sol_ei = best_sol[0]
                    # Reset counter.
                    ext_imm_count = 0

                # Apply Extinction/Immigration operator.
                generation = ext_imm(best_sol, fundam_params, N_pop)

                # Reset best solution counter.
                best_sol_count = 0

        else:
            # For plotting purposes. Save index where a new best solution
            # was found.
            new_bs_indx.append([i])
            # Update best solution for passing along in the 'Elitism' block.
            best_sol = generation[:N_el]
            # Reset counter.
            best_sol_count = 0

        # For plotting purposes.
        lkl_old[0].append(lkl[0])
        # Discard large values associated with empty arrays from mean.
        lkl_old[1].append(np.mean(np.asarray(lkl)[np.asarray(lkl) < 9.9e08]))

        if flag_print_perc:
            percentage_complete = (100.0 * (i + 1) / N_gen)
            while len(milestones) > 0 and percentage_complete >= milestones[0]:
                print (" {:>3}%  L={:.1f} ({:g}, {:g}, {:g}, {:g}, {:g},"
                       " {:g})".format(milestones[0], lkl[0], *generation[0]))
                # Remove that milestone from the list.
                milestones = milestones[1:]

        # print i, generation[0], lkl[0], len(model_done[0])

    isoch_fit_params = {
        'best_sol': generation[0], 'lkl_old': lkl_old,
        'new_bs_indx': new_bs_indx, 'model_done': model_done}

    return isoch_fit_params
