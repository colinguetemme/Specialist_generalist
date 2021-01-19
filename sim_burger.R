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

amplitude <- 0.5
period <- 24
stoch_magnitude <- 0
strength_selec <- 5
mut_rate <- 10^(-6)

num_loci <- 5
t_max <- 50

# One simulation over the t_max period

initialisation <- genetic_sim(num_loci = num_loci)
gamette_values <- initialisation[[1]]
gamette_distr <- initialisation[[2]]
gamette_scheme <- initialisation [[3]]
rec_rate <- initialisation[[4]]

dtf_rec <- recombination(num_loci = num_loci, rec_rate = rec_rate)

opti_G <- get_optimum(t_max = t_max, amplitude = amplitude, period = period, 
            stoch_magnitude = stoch_magnitude)

distr_over_time <- matrix(0, ncol = length(gamette_distr), nrow = t_max + 1)
distr_over_time[1, ] <- gamette_distr
for (t in 1:t_max){
  gamette_distr <- new_distributions(t, num_loci, gamette_distr, gamette_values,
                                     opti_G, dtf_rec, mut_rate)
  distr_over_time[t+1, ] <- gamette_distr
}


