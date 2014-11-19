# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:43:14 2014

@author: gabriel
"""
from get_likelihood import isoch_likelihood as i_l
import random
import numpy as np


def encode(n_bin, p_delta, p_mins, int_popul):
    '''
    Encode the solutions into binary string chromosomes to be bred.
    '''

    chromosomes = []
    for sol in int_popul:
        # Convert floats to binary strings.
        p_binar = []
        for i, p_del in enumerate(p_delta):
            p_binar.append(str(bin(int(((sol[i] - p_mins[i]) / p_delta[i]) *
                (2 ** n_bin))))[2:].zfill(n_bin))

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


def crossover(chromosomes, p_cross, cr_sel):
    '''
    Applies the crosssover operator over each chromosome.
    '''
    cross_chrom = []
    # Take two chromosomes at a time.
    for chrom_pair in chunker(chromosomes, 2):
        r = random.random()
        if r <= p_cross:

            if cr_sel == '1P':
                # Select one random crossover point.
                cp = random.randint(0, len(chrom_pair[0]))
                # Apply crossover on these two chromosomes.
                cross_chrom.append(chrom_pair[0][:cp] + chrom_pair[1][cp:])
                cross_chrom.append(chrom_pair[1][:cp] + chrom_pair[0][cp:])
            elif cr_sel == '2P':
                # Select two random crossover points.
                cp0, cp1 = np.sort(random.sample(range(0, len(chrom_pair[0])),
                                                 2))
                # Apply crossover on these two chromosomes.
                cross_chrom.append(chrom_pair[0][:cp0] +
                chrom_pair[1][cp0:cp1] + chrom_pair[0][cp1:])
                cross_chrom.append(chrom_pair[1][:cp0] +
                chrom_pair[0][cp0:cp1] + chrom_pair[1][cp1:])
        else:
            # Skip crossover operation.
            cross_chrom.append(chrom_pair[0])
            cross_chrom.append(chrom_pair[1])

    return cross_chrom


def mutation(cross_chrom, p_mut):
    '''
    Applies the mutation operator over random genes in each chromosome.
    '''
    # For each chromosome flip their genes according to the probability p_mut.
    for i, elem in enumerate(cross_chrom):
        cross_chrom[i] = ''.join(char if random.random() > p_mut else
        str(1 - int(char)) for char in elem)

    return cross_chrom


def decode(param_values, n_bin, p_delta, p_mins, mut_chrom):
    '''
    Decode the chromosomes into its real values to be evaluated by the
    objective function.
    '''

    # Initialize empty list with as many sub-lists as parameters.
    p_lst = [[] for _ in range(len(param_values))]

    for chrom in mut_chrom:
        # Split chromosome string.
        chrom_split = [chrom[i:i + n_bin] for i in xrange(0, len(chrom), n_bin)]

        # Convert each binary into an integer for all parameters.
        b2i = [int(i, 2) for i in chrom_split]

        # Map integers to the real parameter values.
        for i, p_del in enumerate(p_delta):
            # Map integer to the real parameter value.
            p_r = p_mins[i] + (b2i[i] * p_del / (2 ** n_bin))
            # Find the closest value in the parameters list.
            p = min(param_values[i], key=lambda x: abs(x - p_r))
            # Store this last value.
            p_lst[i].append(p)

    return p_lst


def elitism(best_sol, p_lst):
    '''
    Pass the best n_el solutions unchanged to the next generation.
    '''

    # Append the n_el best solutions to the beginning of the list, pushing
    # out n_el solutions located on the end (right) of the list.
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


def evaluation(err_lst, obs_clust, completeness, isoch_list, param_values,
                 p_lst, st_d_bin_mr, model_done, cmd_sel):
    '''
    Evaluate each model in the objective function to obtain the fitness of
    each one.
    '''

    likel_lst, generation_list = [], []
    # Process each model selected.
    for model in zip(*p_lst):

        # Check if this model was already processed.
        if model in model_done[0]:
            # Get likel_val value for this isochrone.
            likelihood = model_done[1][model_done[0].index(model)]
        else:
            # Metallicity and age indexes to identify isochrone.
            m_i = param_values[0].index(model[0])
            a_i = param_values[1].index(model[1])
            isochrone = isoch_list[m_i][a_i]

            # Call likelihood function for this model.
            likelihood = i_l(err_lst, obs_clust, completeness, st_d_bin_mr,
                            isochrone, model, cmd_sel)
            # Append data identifying the isochrone and the obtained
            # likelihood value to this *persistent* list.
            model_done[0].append(model)
            model_done[1].append(likelihood)

        # Append data to the lists that will be erased with each call
        # to this function.
        generation_list.append(model)
        likel_lst.append(likelihood)

    # Sort according to the likelihood list. This puts the best model (ie:
    # the one with the minimum likelihood value) first.
    generation = [x for y, x in sorted(zip(likel_lst, generation_list))]
    # Sort list in place putting the likelihood minimum value first.
    likel_lst.sort()

    return generation, likel_lst, model_done


def random_population(param_values, n_ran):
    '''
    Generate a random set of parameter values to use as a random population.
    '''
    # Pick n_ran initial random solutions from each list storing all the
    # possible parameters values. These lists store real values.
    p_lst = []
    for param in param_values:
        p_lst.append([random.choice(param) for _ in range(n_ran)])

    return p_lst


def ext_imm(best_sol, param_values, n_pop):
    '''
    Append a new random population to the best solution so far.
    '''

    # Generate (n_pop-n_el) random solutions.
    n_ran = n_pop - len(best_sol)
    p_lst_r = random_population(param_values, n_ran)

    # Append immigrant random population to the best solution.
    generation_ei = best_sol + zip(*p_lst_r)

    return generation_ei


def num_binary_digits(param_rs):
    '''
    Store parameters ranges and calculate the minimum number of binary digits
    needed to encode the solutions.
    '''

    p_mins, p_delta, p_step = [], [], []
    for param in param_rs:
        p_mins.append(param[0])  # Used by the encode operator.
        p_delta.append(param[1] - param[0])  # max - min value
        p_step.append(param[2])  # step

    p_interv = np.array(p_delta) / np.array(p_step)

    # Number of binary digits used to create the chromosomes.
    # The max function prevents an error when all parameter ranges are set to
    # a unique value in which case the np.log is a negative float.
    n_bin = max(int(np.log(max(p_interv)) / np.log(2)) + 1, 1)

    return n_bin, p_delta, p_mins


def gen_algor(flag_print_perc, err_lst, obs_clust, completeness, ip_list,
    st_d_bin_mr, ga_params, cmd_sel):
    '''
    Genetic algorithm adapted to find the best fit model-obervation.
    '''

    # Unpack.
    isoch_list, param_values, param_rs = ip_list
    n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es = ga_params
    # Check if n_pop is odd. If it is sum 1 to avoid conflict if cr_sel
    # '2P' was selected.
    n_pop += n_pop % 2

    # Get number of binary digits to use.
    n_bin, p_delta, p_mins = num_binary_digits(param_rs)

    # Fitness.
    # Rank-based breeding probability. Independent of the fitness values,
    # only depends on the total number of chromosomes n_pop and the fitness
    # differential fdif.
    fitness = [1. / n_pop + fdif * (n_pop + 1. - 2. * (i + 1.)) /
        (n_pop * (n_pop + 1.)) for i in range(n_pop)]

    ### Initial random population evaluation. ###
    p_lst_r = random_population(param_values, n_pop)

    # Stores parameters of the solutions already processed and the likelihhods
    # obtained.
    model_done = [[], []]

    # Evaluate initial random solutions in the objective function.
    generation, lkl, model_done = evaluation(err_lst, obs_clust,
        completeness, isoch_list, param_values, p_lst_r,
        st_d_bin_mr, model_done, cmd_sel)

    # Store best solution for passing along in the 'Elitism' block.
    best_sol = generation[:n_el]

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
    milestones = [20, 40, 60, 80, 100]
    # Begin processing the populations up to n_gen generations.
    for i in range(n_gen):

        #### Selection/Reproduction ###
        # Select chromosomes for breeding from the current generation of
        # solutions according to breed_prob to generate the intermediate
        # population.
        int_popul = selection(generation, fitness)

        # Encode intermediate population's solutions into binary chromosomes.
        chromosomes = encode(n_bin, p_delta, p_mins, int_popul)

        #### Breeding ###
        # Pair chromosomes by randomly shuffling them.
        random.shuffle(chromosomes)

        # Apply crossover operation on each subsequent pair of chromosomes
        # with a p_cross probability (crossover probability)
        cross_chrom = crossover(chromosomes, p_cross, cr_sel)

        # Apply mutation operation on random genes for every chromosome.
        mut_chrom = mutation(cross_chrom, p_mut)

        ### Evaluation ###
        # Decode the chromosomes into solutions to form the new generation.
        p_lst_d = decode(param_values, n_bin, p_delta, p_mins, mut_chrom)

        # Elitism: make sure that the best n_el solutions from the previous
        # generation are passed unchanged into this next generation.
        p_lst_e = elitism(best_sol, p_lst_d)

        # Evaluate each new solution in the objective function and sort
        # according to the best solutions found.
        generation, lkl, model_done = evaluation(err_lst, obs_clust,
            completeness, isoch_list, param_values, p_lst_e, st_d_bin_mr,
            model_done, cmd_sel)

        ### Extinction/Immigration ###
        # If the best solution has remained unchanged for n_ei
        # generations, remove all chromosomes but the best ones (extinction)
        # and fill with random new solutions (immigration).

        # Check if new best solution is equal to the previous one.
        if generation[0] == best_sol[0]:
            # Increase counter.
            best_sol_count += 1

            # Check how many times the best_sol has remained unchanged.
            # If the number equals n_ei, apply Extinction/Immigration operator.
            if best_sol_count == n_ei:

                ### Exit switch. ##
                # If n_es runs of the Ext/Imm operator have been applied with
                # no changes to the best solution, apply the exit switch,
                # ie: exit the GA.
                if best_sol[0] == best_sol_ei:
                    # Increase Ext/Imm operator counter.
                    ext_imm_count += 1
                    if ext_imm_count == n_es:
                        # Exit generations loop.
                        break
                else:
                    # Update best solution.
                    best_sol_ei = best_sol[0]
                    # Reset counter.
                    ext_imm_count = 0

                # Apply Extinction/Immigration operator.
                generation = ext_imm(best_sol, param_values, n_pop)

                # Reset best solution counter.
                best_sol_count = 0

        else:
            # For plotting purposes. Save index where a new best solution
            # was found.
            new_bs_indx.append([i])
            # Update best solution for passing along in the 'Elitism' block.
            best_sol = generation[:n_el]
            # Reset counter.
            best_sol_count = 0

        if flag_print_perc:
            percentage_complete = (100.0 * (i + 1) / n_gen)
            while len(milestones) > 0 and percentage_complete >= milestones[0]:
                print ("  {}% done \t ({:g}, {:g}, {:g}, {:g}, {:g},"
                " {:g})".format(milestones[0], *generation[0]))
                # Remove that milestone from the list.
                milestones = milestones[1:]

        # For plotting purposes.
        lkl_old[0].append(lkl[0])
        lkl_old[1].append(np.mean(lkl))

        #print i, generation[0], lkl[0], len(model_done[0])

    isoch_fit_params = [generation[0], lkl_old, new_bs_indx, model_done]

    return isoch_fit_params
