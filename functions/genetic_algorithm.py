# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:43:14 2014

@author: gabriel
"""
from isoch_likelihood import isoch_likelihood as i_l

import random
import numpy as np


def encode(mm_m, mm_a, mm_e, mm_d, n, int_popul):
    '''
    Encode the solutions into binary string chromosomes to be bred.
    '''
    
    delta_m, delta_a, delta_e, delta_d = (mm_m[1]-mm_m[0]), (mm_a[1]-mm_a[0]),\
    (mm_e[1]-mm_e[0]), (mm_d[1]-mm_d[0])
    chromosomes = []
    for sol in int_popul:
        m_binar = str(bin(int(((sol[0]-mm_m[0])/delta_m)*(2**n))))[2:].zfill(n)
        a_binar = str(bin(int(((sol[1]-mm_a[0])/delta_a)*(2**n))))[2:].zfill(n)
        e_binar = str(bin(int(((sol[2]-mm_e[0])/delta_e)*(2**n))))[2:].zfill(n)
        d_binar = str(bin(int(((sol[3]-mm_d[0])/delta_d)*(2**n))))[2:].zfill(n)
        chrom = m_binar + a_binar + e_binar + d_binar
        print len(chrom), m_binar, a_binar, e_binar, d_binar, sol
        chromosomes.append(chrom)
    
    return chromosomes
    
    
def decode():
    '''
    Decode the chromosomes into its real values to be evaluated by the
    objective function.
    '''
    return
    
    

def selection(generation, breed_prob):
    '''
    Select np random chromosomes from the chromose list passed according to
    the breeding probability given by their fitness.
    '''
    select_chrom = []
    # Draw n_pop random numbers uniformelly distributed between [0,1)
    ran_lst = np.random.uniform(0, 1, len(generation))
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
        m, a = ma_lst[indx]
        isochrone = isoch_list[m][a]
        d = d_lst[indx]
        likel_val = i_l(sys_select, isochrone, e, d, obs_clust, mass_dist)
        likelihood.append(likel_val)
        generation_list.append([isoch_ma[m][a][0], isoch_ma[m][a][1], e,d])
    # Sort according to the likelihood list.
    generation = [x for y, x in sorted(zip(likelihood, generation_list))]
    # Sort list in place putting the likelihood minimum value first.
#    likelihood.sort()
    
#    return likelihood, generation
    return generation
    


def gen_algor(sys_select, obs_clust, isoch_list, isoch_ma, isoch_ed, mass_dist,
              ranges_steps, n_pop, n_gen, fdif):
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
    print delta_m,step_m, delta_a,step_a, delta_e,step_e, delta_d,step_d
    n_bin = int(np.log(max((delta_m/step_m), (delta_a/step_a), (delta_e/step_e), (delta_d/step_d)))/np.log(2))+1
    print n_bin
    
    
    ### Initial random population evaluation. ###
    
    # Pick n_pop initial random solutions from each list storing all the possible
    # parameters values. The lists store indexes.
    e_lst = [random.choice(isoch_ed[0]) for _ in range(n_pop)]
    d_lst = [random.choice(isoch_ed[1]) for _ in range(n_pop)]
    # Flat array so every [metal,age] combination has the same probability
    # of being picked.
    ma_flat = [(i, j) for i in range(len(isoch_ma)) for j in range(len(isoch_ma[i]))]
    ma_lst = [random.choice(ma_flat) for _ in range(n_pop)]
    
    
    # Evaluate each random solution in the objective function.
    generation = fitness_eval(sys_select, isoch_list, obs_clust, mass_dist,
                              isoch_ma, ma_lst, e_lst, d_lst)

    # Rank-based breeding probability. Independent of the fitness values,
    # only depends on the total number of chromosomes and the fitness
    # differential fdif.
    breed_prob = [1./n_pop + fdif*(n_pop+1.-2.*(i+1.))/(n_pop*(n_pop+1.)) \
    for i in range(n_pop)]    
    
    # Begin processing the populations up to n_gen generations.
    for i in range(n_gen):
    
       
        #### Selection/Reproduction ###
        
        # Elitism: pass the best solution to ensure it breeds.
        int_popul0 = [generation[0]]
        
        # Select (n_pop-n_elit) chromosomes for breeding from the above
        # generation of solutions according to breed_prob to generate the
        # intermediate population.
        int_popul = int_popul0 + selection(generation[1:], breed_prob[1:])
        
        # Encode intermediate population's solutions into binary chromosomes.
        chromosomes = encode(mm_m, mm_a, mm_e, mm_d, n_bin, int_popul)
        print chromosomes
        raw_input()
        
        #### Breeding ###
        
        # Randomly pair chromosomes.
        
        # Apply crossover operation on random pair of chromosomes. Select only
        # 100*p_cross% of pairs to apply the crossover to, where p_cross is the
        # crossover probability.
        
        # Apply mutation operation on random chromosomes.
        
        
        ### Evaluation/fitness ###
        # Evaluate the N random chromosomes in the evaluation/objective function
        # (likelihood)
        
        # Decode the chromosomes into solutions first.
        decode()
        
        # Evaluate each random solution in the objective function.



