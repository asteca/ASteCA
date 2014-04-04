#####################################################################################
##### OCAAT input parameters file ####
#
#
# This file contains the values of the parameters used by OCAAT.
# Any empty line or lines that begin with a # character will be ignored
# and as many as needed can be added to this file.
#
# The blocks can be moved around as the user sees fit.
#
# * DO NOT modify the ID strings in the first column of each line, these
# are used by the code to identify the parameters in that line.
# * DO NOT write less or more values than those requested in each line
# or the code will halt.
#
# Set parameters values below.
#####################################################################################


# Mode used to run the code: a (auto) / s (semi) / m (manual):
# - Automatic mode does not require user intervention of any kind.
# - Semi takes it input values from the 'clusters_input.dat' file.
# - Manual mode will ask for some values (center, radius, background) but
#   will also give you the chance to have them automatically calculated.
MO  s


# Move file once it is processed?
# The cluster's photometric file will be moved to the folder
# defined below if this flag is set to 'True'.
MF  false
# Path where the sub-folders are moved once processed if the flag
# above is set to 'True'.
CP_d  /path/


# Create output figure?
# True / False. If 'False', no output figure will be generated.
#    flag_make_plot  format  DPI
MP             True     png  150


# Column numbers identifying the needed data as stored in the clusters'
# photometric data files.
# id_star x_coor y_coor magnitude e_mag color e_col
PD      0      1      2         3     4     5     6


#################################################
# Photometric system parameters.

# CMD selection. Select a CMD to be used by the code.
# 1: (B-V) vs V   ; UBVI (Johnson)
# 2: (V-I) vs V   ; UBVI (Johnson)
# 3: (C-T1) vs T1 ; Washington
#    selection
CMD          1

# Select set of isochrones to use.
# iso_select: MAR (Marigo) / PAR (PARSEC)
#    iso_select
PS          MAR

# Maximum and minimum axis values for the CMD plots.
#   col_min col_max mag_min mag_max
MM       -1       4       7      30
#################################################


# Center.
# Bin widths for the cluster's center determination, in pixels. The
# minimum value will be used to establish the center if 'MO' is
# in auto mode.
# mode: auto / manual. If 'auto' is selected the code will calculate the
# width of the bins as 1%, 2% and 3% of the range spanned by
# either x,y axis, whichever spans the minimum. If 'manual' is set,
# the values will be taken from bin0 bin1 bin3, otherwise these values
# are of no importance.
#   mode  bin0 bin1 bin2 bin3
CC  auto    35   50   75  100


# Radius.
# Number of points that must be found below a minimum delta around the
# background value. This condition indicates that the star density
# has stabilized.
# mode: auto / manual. If 'auto' is selected the 'delta_step' and 'points'
# parameters will set as 5% and 10% of the total number of points in the
# radial density profile, respectively.
#    mode   delta_step(%) points
CR   auto            1       5


# Errors.
# Rejection of stars according to their errors.
# be: bright end, the magnitude value that will be added to the
# brightest star to determine the range where stars will be accepted
# if their errors are lower than 'be_e' (see below).
# be_e: maximum error value that stars in the bright end range can
# have without being rejected.
# e_max: maximum error value a star can have in either its magnitude or
# color without being rejected. 
#   be  be_e  e_max
ER  2.   0.05    0.3


# Field regions.
# regions: number of field regions around the cluster that the code will
# attempt to define.
#   regions
GR       5


# p-value test.
# Test to determine the probability of the cluster being a physical entity
# instead of a random overdensity of stars, based on comparing its CMD with
# field regions CMDs.
# flag: True / False, switch to apply or skip the test.
# runs: auto / manual, the number of times the cluster-field regions will be
# compared. If 'auto' is set, the code will use an integer value of
# 100/(areas-1). If 'manual' is set then 'num' will be used, otherwise this
# number is irrelevant.
#   flag     runs   num
PV  false     auto    1


# Bayesian decontamination algorithm.
# mode: auto / manual / read / skip.
# - If 'auto' is set the code will use the number of runs 'runs'
# as defined internally ie: 1000.
# - If 'manual' is set this value is taken from this file (below),
# otherwise this value is irrelevant.
# - If 'read' is set the code will attempt to read the probabilities from
# an existing file. This file should have the same format as the output
# '*_memb.dat' file created by the code and should be located in the folder
# defined by the ID 'CP0' above.
# - If 'skip' is selected all stars will be assigned an equal probability
# of 1. and the algorithm skipped altogether.
#   mode     runs
DA  skip     5


#################################################
# Best fit.
# run: True / False. Flag that determines if the best synthetic cluster
# fitting process is to be run.
# algorithm: select algorithm to use in best fitting process,
# 'brute' (for Brute Force) / 'genet' (for Genetic Algorithm)
# If 'brute' is selected, 'boot_num' is irrelevant.
# boot_num: number of times the bootstrap with replacement
# process will run. Minimum value is 2.
#     run  algorithm   boot_num
BF   false     genet         10

# Ranges and steps for the parameters to be fitted.
# Metallicity     z_min  z_max
PS_m             0.0005   0.03
# Age    age_min age_max
PS_a         6.6    10.1
# E(B-V)    e_bv_min e_bv_max e_bv_step
#PS_e              0.     0.4      0.01
#PS_e              0.35     0.35      0.05
PS_e              0.     1.      0.1
# Distance modulus    dis_mod_min dis_mod_max dis_mod_step
#PS_d                          18.       19.         0.05
#PS_d                          10.5       12.5         0.1
#MASSCLEAN
PS_d                          8.5       14.         0.2

# Synthetic cluster.
# IMF:  chabrier_2001 / kroupa_1993 / kroupa_2002
#             IMF   total_mass   frac_binar  mass_ratio
SC  chabrier_2001         5000          0.5         0.7

# Genetic Algorithm parameters.
# n_pop: number of chromosomes in the population.
# n_gen: number of generations to process.
# fdif: Fitness differential. Establishes the 'selection pressure' for the
# algorithm.
# p_cros: crossover probability.
# (Crossover operator)
# cr_sel: select 1-point (1P) or 2-point (2P) crossover.
# p_mut: mutation probability.
# (Elitism operator)
# n_el: number of best solutions to pass unchanged to the next generation
# n_ei: number of generations allowed to run with the same best solution
# before applying the Extinction/Immigration operator.
# (Exit switch)
# n_es:number of times the Extinction/Immigration operator is allowed to run
# returning the same best solution before exiting the generations loop.
#
#   n_pop  n_gen fdif p_cross cr_sel p_mut n_el n_ei n_es
GA    100  1000    1.    0.85     2P  0.01    1   75   10
#################################################
