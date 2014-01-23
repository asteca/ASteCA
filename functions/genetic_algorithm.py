# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:43:14 2014

@author: gabriel
"""
from isoch_likelihood import isoch_likelihood as i_l

import random
import numpy as np
import time



def encode(mm_m, mm_a, mm_e, mm_d, n, int_popul):
    '''
    Encode the solutions into binary string chromosomes to be bred.
    '''
    
    delta_m, delta_a, delta_e, delta_d = (mm_m[1]-mm_m[0]), (mm_a[1]-mm_a[0]),\
    (mm_e[1]-mm_e[0]), (mm_d[1]-mm_d[0])
    chromosomes = []
    for sol in int_popul:
        # Convert floats to binary strings.
        m_binar = str(bin(int(((sol[0]-mm_m[0])/delta_m)*(2**n))))[2:].zfill(n)
        a_binar = str(bin(int(((sol[1]-mm_a[0])/delta_a)*(2**n))))[2:].zfill(n)
        e_binar = str(bin(int(((sol[2]-mm_e[0])/delta_e)*(2**n))))[2:].zfill(n)
        d_binar = str(bin(int(((sol[3]-mm_d[0])/delta_d)*(2**n))))[2:].zfill(n)
        chrom = m_binar + a_binar + e_binar + d_binar
#        print len(chrom), m_binar, a_binar, e_binar, d_binar, sol
        chromosomes.append(chrom)
    
    return chromosomes
    
    
def chunker(seq, size):
    '''
    Helper function for 'crossover'. Returns elemenst in list in "chunks"
    rather than one at a time.
    '''
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))
    
    
def crossover(chromosomes, p_cross):
    '''
    Applies the crosssover operator over each chromosome.
    '''
    cross_chrom= []
    # Take two chromosomes at a time.
    for chrom_pair in chunker(chromosomes, 2):
        r = random.random()
        if r<= p_cross:
            # Select random crossover point.
            cp = random.randint(0, len(chrom_pair[0]))
            cross_chrom.append(chrom_pair[0][:cp]+chrom_pair[1][cp:])
            cross_chrom.append(chrom_pair[1][:cp]+chrom_pair[0][cp:])
        else:
            cross_chrom.append(chrom_pair[0])
            cross_chrom.append(chrom_pair[1])
        
    return cross_chrom
    
    
def mutation(cross_chrom, p_mut):
    '''
    Applies the mutation operator over random genes in each chromosome.
    '''
    #                 
    for i,elem in enumerate(cross_chrom):
        cross_chrom[i] = ''.join(char if random.random()>p_mut else str(1-int(char)) for char in elem)
    
    return cross_chrom
    
    
def decode(mm_m, mm_a, mm_e, mm_d, n_bin, isoch_ma, isoch_ed, mut_chrom):
    '''
    Decode the chromosomes into its real values to be evaluated by the
    objective function.
    '''
    import itertools
#    import bisect
    
    delta_m, delta_a, delta_e, delta_d = (mm_m[1]-mm_m[0]), (mm_a[1]-mm_a[0]),\
    (mm_e[1]-mm_e[0]), (mm_d[1]-mm_d[0])    
    
    # Flat out metallicity and ages list.
    flat_ma = zip(*list(itertools.chain(*isoch_ma)))
#    met_lst = [i[0][0] for i in isoch_ma]
    
    ma_lst, e_lst, d_ls = [], [], []
    for chrom in mut_chrom:
        # Split chromosome string.
        chrom_split = [chrom[i:i+n_bin] for i in xrange(0, len(chrom), n_bin)]
        
        # Convert binary to integer for each parameter.
        km, ka, ke, kd = (int(i, 2) for i in chrom_split)
        
        # Map integers to the real parameter values.
        xm = mm_m[0]+(km*delta_m/(2**n_bin))
        xa = mm_a[0]+(ka*delta_a/(2**n_bin))
        xe = mm_e[0]+(ke*delta_e/(2**n_bin))
        xd = mm_d[0]+(kd*delta_d/(2**n_bin))
        
#        print xm, xa, xe, xd
#        tik=time.time()
        # Find the closest values in the parameters list and store its index
        # in the case of metallicity and age and real values for extinction
        # and distance modulus.
        m = min(flat_ma[0], key=lambda x:abs(x-xm))
        a = min(flat_ma[1], key=lambda x:abs(x-xa))
        e = min(isoch_ed[0], key=lambda x:abs(x-xe))
        d = min(isoch_ed[1], key=lambda x:abs(x-xd))
#        print ' ', m,a
        # Find the indexes for these metallicity and age values.
        [m, a] = next(((i,j) for i,x in enumerate(isoch_ma) for j,y in enumerate(x) if y == [m, a]), None)
        
#        print ' ', m,a,e,d,time.time()-tik
#        
#        tik=time.time()
#        m = bisect.bisect(met_lst, xm)
#        a = bisect.bisect(flat_ma[1][:len(isoch_ma[0])], xa)
#        e = bisect.bisect(isoch_ed[0], xe)
#        d = bisect.bisect(isoch_ed[1], xd)
#        print ' ', m,a
#        print ' ', isoch_ma[m][a][0], isoch_ma[m][a][1]
#        print ' ', m,a,isoch_ed[0][e],isoch_ed[1][d],time.time()-tik, '\n'
        
        # Append indexes (for m,a) and real values (for e,d) to lists.
        ma_lst.append([m, a])
        e_lst.append(e)
        d_ls.append(d)
        
    return ma_lst, e_lst, d_ls
    
    

def selection(generation, breed_prob):
    '''
    Select random chromosomes from the chromose list passed according to
    the breeding probability given by their fitness.
    '''

    select_chrom = []
    # Draw n_pop random numbers uniformelly distributed between [0,1)
    ran_lst = np.random.uniform(0, sum(breed_prob), len(generation))
    # For each of these numbers, obtain the corresponding likelihood from the
    # breed_prob CDF.
    for r in ran_lst:
        s = 0.
        for sol,num in zip(*[generation, breed_prob]):
            s += num
            if s >= r:
                select_chrom.append(sol)
                break

    return select_chrom



def fitness_eval(sys_select, isoch_list, obs_clust, mass_dist, isoch_ma,
                 ma_lst, e_lst, d_lst):
    '''
    Evaluate each n_pop random solution in the objective function to obtain
    the fitness of each solution.
    '''
    likelihood, generation_list = [], []
    for indx, e in enumerate(e_lst):
        # Get isochrone with m,a values.
        m, a = ma_lst[indx]
        isochrone = isoch_list[m][a]
        d = d_lst[indx]
        # Call likelihood function with m,a,e,d values.
        likel_val = i_l(sys_select, isochrone, e, d, obs_clust, mass_dist)
        likelihood.append(likel_val)
        generation_list.append([isoch_ma[m][a][0], isoch_ma[m][a][1], e,d])
    # Sort according to the likelihood list.
    generation = [x for y, x in sorted(zip(likelihood, generation_list))]
    # Sort list in place putting the likelihood minimum value first.
    likelihood.sort()
    
    return generation, likelihood
#    return generation
    


def gen_algor(sys_select, obs_clust, isoch_list, isoch_ma, isoch_ed, mass_dist,
              ranges_steps, n_pop, n_gen, fdif, p_cross, p_mut):
    '''
    Main function.
    
    Genetic algorithm adapted to find the best fit model-obervation.
    
    n_pop: number of chromosomes in the population.
    n_gen: number of generations to process.
    fdif: Fitness differential. Establishes the 'selection pressure' for the
    algorithm.
    '''
    
    # Store parameters ranges and calculate the minimum number of binary digits
    # needed to encode the solutions.
    mm_m, step_m = ranges_steps[0], (isoch_ma[1][0][0]-isoch_ma[0][0][0])
    mm_a, step_a = ranges_steps[1], (isoch_ma[0][1][1]-isoch_ma[0][0][1])
    mm_e, step_e = [ranges_steps[2][0], ranges_steps[2][1]], ranges_steps[2][2]
    mm_d, step_d = [ranges_steps[3][0], ranges_steps[3][1]], ranges_steps[3][2]
    delta_m, delta_a, delta_e, delta_d = (mm_m[1]-mm_m[0]), (mm_a[1]-mm_a[0]),\
    (mm_e[1]-mm_e[0]), (mm_d[1]-mm_d[0])
    n_bin = int(np.log(max((delta_m/step_m), (delta_a/step_a), (delta_e/step_e), (delta_d/step_d)))/np.log(2))+1
    
    
    ### Initial random population evaluation. ###
    
    # Pick n_pop initial random solutions from each list storing all the possible
    # parameters values. These lists store real values.
    e_lst = [random.choice(isoch_ed[0]) for _ in range(n_pop)]
    d_lst = [random.choice(isoch_ed[1]) for _ in range(n_pop)]
    # Flat array so every [metal,age] combination has the same probability
    # of being picked. This list stores indexes.
    ma_flat = [(i, j) for i in range(len(isoch_ma)) for j in range(len(isoch_ma[i]))]
    ma_lst = [random.choice(ma_flat) for _ in range(n_pop)]
    
    
    # Evaluate each random solution in the objective function.
    generation, lkl = fitness_eval(sys_select, isoch_list, obs_clust, mass_dist,
                              isoch_ma, ma_lst, e_lst, d_lst)

    # Rank-based breeding probability. Independent of the fitness values,
    # only depends on the total number of chromosomes and the fitness
    # differential fdif.
    breed_prob = [1./n_pop + fdif*(n_pop+1.-2.*(i+1.))/(n_pop*(n_pop+1.)) \
    for i in range(n_pop)]    
    
    
    # Begin processing the populations up to n_gen generations.
    for i in range(n_gen):
        
        tik0 = time.time()
        #### Selection/Reproduction ###
        
        # Store best solution for passing along in the 'Elitism' block.
        best_sol = [generation[0]]
        
        # Select chromosomes for breeding from the current  generation of
        # solutions according to breed_prob to generate the intermediate
        # population.
        int_popul = selection(generation, breed_prob)

        # Encode intermediate population's solutions into binary chromosomes.
        chromosomes = encode(mm_m, mm_a, mm_e, mm_d, n_bin, int_popul)

        
        #### Breeding ###
        
        # Pair chromosomes by randomly shuffling them.
        random.shuffle(chromosomes)

        # Apply crossover operation on each subsequent pair of chromosomes.
        # Select only 100*p_cross% of pairs to apply the crossover to,
        # where p_cross is the crossover probability.
        cross_chrom = crossover(chromosomes, p_cross)

        # Apply mutation operation on random genes for every chromosome.
        mut_chrom = mutation(cross_chrom, p_mut)
        
        # Elitism: make sure that the best solution from the previous generation
        # is passed unchanged into this next generation.
        best_chrom = encode(mm_m, mm_a, mm_e, mm_d, n_bin, best_sol)
        mut_chrom[0] = best_chrom[0]
        
        
        ### Evaluation/fitness ###
        
        # Decode the chromosomes into solutions to form the new generation.
        tik6 = time.time()
        ma_lst, e_lst, d_lst = decode(mm_m, mm_a, mm_e, mm_d, n_bin, isoch_ma, isoch_ed, mut_chrom)
        print 'decod', time.time()-tik6
        
        # Evaluate each new solution in the objective function and sort
        # according to the best solutions found.
        tik7 = time.time()
        generation, lkl = fitness_eval(sys_select, isoch_list, obs_clust, mass_dist,
                                  isoch_ma, ma_lst, e_lst, d_lst)
        print 'fitne', time.time()-tik7
                                  
        print i, lkl[0], generation[0], time.time()-tik0

    return generation[0]
