# Modelling of the paper of Burger 2002 
# study of a range of temporal variation 

source(file="mod_burger.R", encoding ="UTF-8")

#### Parameters ####

# amplitude is the size of the range of values for the optimum genotype
# period is the length between 2 period
# magnitude important is the noise add to the sin curve describing the optimum 
# genotype
# strength_selec isstrenght of stabilization (noted s)
# viability of an individual with genotype G
# mut_rate is the common mutation rate (equally to pass from gamete j to k)
# num_loci the number of loci on the gamete
# t_max is the time when the simulation stop (in the paper stop when stable)
# rec_rate the recombination rate



amplitude <- 1
period <- 100
stoch_magnitude <- 0.2
strength_selec <- 0.5
mut_rate <- 10^(-6)
num_loci <- 2
t_max <- 200

genetic