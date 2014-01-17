# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:43:14 2014

@author: gabriel
"""

#import random
import numpy as np


def encode():
    '''
    Encode the solutions into binary string chromosomes to be bred.
    '''
    return 
    
    
def decode():
    '''
    Decode the chromosomes into its real values to be evaluated by the
    objective function.
    '''
    return
    
    

def selection(likel_list, breed_prob):
    '''
    Select np random chromosomes from the chromose list passed according to
    the breeding probability given by their fitness.
    '''
    select_chrom = []
    # Draw np random numbers uniformelly distributed between [0,1)
    ran_lst = np.random.uniform(0, 1, len(likel_list))
    # For each of these numbers, obtain the corresponding element from the CDF.
    for r in ran_lst:
        s = 0.
        for prob,num in zip(*[likel_list, breed_prob]):
            s += num
            if s >= r:
                select_chrom.append(prob)
                break

    return select_chrom



def gen_algor(obs_clust, isochrones, isoch_params, dist_mass, n_pop, n_gen, fdif):
    '''
    Main function.
    
    Genetic algorithm adapted to find the best fit model-obervation.
    
    n_pop: number of chromosomes in the population.
    n_gen: number of generations to process.
    fdif: Fitness differential. Establishes the 'selection pressure' for the
    algorithm.
    '''
    
    # Initial random population evaluation.
    
    # Pick N initial solutions at random from each
    # list storing all the possible parameters values.
    
    # Evaluate each solution in the objective function.
    
    
    # Rank-based breeding probability. Independent of the fitness values,
    # only depends on the toal number of chromosomes and the fitness
    # differential fdif.
    breed_prob = [1./n_pop + fdif*(n_pop+1.-2.*(i+1.))/(n_pop*(n_pop+1.)) for i in range(n_pop)]    
    
    # Begin processing the populations up to n_gen generations.
    for i in range(n_gen):
    
        #### Breeding ###
        
        # Encode solutions into chromosomes.
        chromosomes = encode(generation)
        
        # Apply crossover operation on random chromosomes.
        
        # Apply mutation operation on random chromosomes.
        
        
        ### Evaluation/fitness ###
        
        # Evaluate the N random chromosomes in the evaluation/objective function
        # (likelihood)
        # Decode the chromosomes into solutions first.
        decode()
        # Calculate the true fitness value for each solution from the objective
        # function.
        likel_list = [1.3, 0.2, 0.7, 5.6, 0.35, 20., 2.7, 0.64, 1.5, 6.6, 8.2, 1.3, 7.4, 3.0, 0.1, 1.2, 0.35, 0.21, 3.1]
        # Sort list in place. Minimum value goes first.
        likel_list.sort()


        #### Selection ###
        
        # Select np chromosomes for breeding from the CDF to generate the
        # intermediate population.
        select_chrom = selection(likel_list, breed_prob)