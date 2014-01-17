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



def gen_algor(N_gen, fdif):
    '''
    Main function.
    
    Genetic algorithm adapted to find the best fit model-obervation.
    N_ga: number of generations to process.
    fdif: Fitness differential. Establishes the 'selection pressure' for the
    algorithm.
    '''
    
    # Translate space of all possible solutions to binary format.
    
    # Pick N initial solutions (genotypes/chromosomes) at random from each
    # list of the parameters possible values.
    
    # Encode these initial parameters values into binary arrays and combine the
    # encoded parameters into N binary chromosomes.
    
    
    # Begin genetic algorithm process.
    for i in range(N_gen):
    
        ### Evaluation/fitness ###
        
        # Evaluate the N random chromosomes in the evaluation/objective function
        # (likelihood)
        # Decode the chromosomes into solutions first.
        decode()
        # Calculate the true fitness value for each solution from the objective
        # function.
        likel_list = [1.3, 0.2, 0.7, 5.6, 0.35, 20., 2.7, 0.64, 1.5, 6.6, 8.2, 1.3, 7.4, 3.0, 0.1, 1.2, 0.35, 0.21, 3.1]
        np = len(likel_list)    
        # Sort list in place. Minimum value goes first.
        likel_list.sort()
        
        # Rank-based breeding probability. Independent of the fitness values,
        # only depends on the toal number of chromosomes and the fitness
        # differential fdif.
        breed_prob = [1./np + fdif*(np+1.-2.*(i+1.))/(np*(np+1.)) for i in range(np)]
       
        
        #### Selection ###
        
        # Select np chromosomes for breeding from the CDF to generate the
        # intermediate population.
        select_chrom = selection(likel_list, breed_prob)
    
        
    
    gen_algor()
        
        #### Breeding ### 
        
        # Apply crossover operation on random chromosomes.
        
        # Apply mutation operation on random chromosomes.
        
        
        ### Decoding###
        
        # Translate binary form to normal naming to process this new generation.