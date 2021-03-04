# Modelling of the paper of Burger 2002 
# study of a range of temporal variation 

library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))  

rm(list = ls())
source(file="mod_burger.R", encoding ="UTF-8")

#### Parameters ####

# amplitude is the size of the range of values for the optimum genotype
# period is the length between 2 period
# magnitude important is the noise add to the sin curve describing the optimum 
# genotype
# strength_selec is strenght of stabilization (noted s)
# viability of an individual with genotype G
# mut_rate is the common mutation rate (equally to pass from gamete j to k)
# num_loci the number of loci on the gamete
# t_max is the time when the simulation stop (in the paper stop when stable)
# rec_rate the recombination rate

# few basic quick model are available here

# 1 - THE BASIC MODEL with no stochasticity, no mutation, no recombination, and just
# two loci
# 2 - THE FULL MODEL where you add all what is missing in the previous and 4 loci

quick_params_selection <- 1

if (quick_params_selection == 1){
  amplitude <- 0.5
  period <- 100
  stoch_magnitude <- 0.1
  strength_selec <- 10
  mut_rate <- 0.0002
  num_loci <- 4
  t_max <- 4 * period
}

if (quick_params_selection == 2){
  amplitude <- 0.5
  period <- 50
  stoch_magnitude <- 0.2
  strength_selec <- 5
  mut_rate <- 5*10^(-5)
  num_loci <- 4
  t_max <- 2 * period
}

if (quick_params_selection == 3){
  amplitude <- 0.5
  period <- 100
  stoch_magnitude <- 0
  strength_selec <- 5
  mut_rate <- 5*10^(-5)
  num_loci <- 6
  t_max <- 2 * period
}


# One simulation over the t_max period

one_simul <- function(amplitude, period, stoch_magnitude, strength_selec,
                      mut_rate, num_loci, t_max){
  
  initialisation <- genetic_sim(num_loci = num_loci)
  gamete_values <- initialisation[[1]]
  gamete_distr <- initialisation[[2]]
  gamete_scheme <- initialisation [[3]]
  rec_rate <- initialisation[[4]]
  
  #Maybe put inside genetic_sim
  Vmax <- 0.5 * sum(gamete_values^2)

  
  dtf_rec <- recombination(num_loci = num_loci, rec_rate = rec_rate)

  opti_G <- get_optimum(t_max = t_max, amplitude = amplitude, period = period, 
                        stoch_magnitude = stoch_magnitude)
  
  distr_over_time <- matrix(0, ncol = length(gamete_distr), nrow = t_max+1)
  distr_over_time[1, ] <- gamete_distr
  mean_fitness_over_time <- rep(0, t_max)
  for (t in 1:t_max){
    
    new_results <- new_distributions(num_loci, distr_over_time[t, ], gamete_values,
                                       gamete_scheme, opti_G[t], dtf_rec, mut_rate)
    distr_over_time[t+1, ] <- new_results[[1]]
    mean_fitness_over_time[t] <- new_results[[2]]
  }
  
  results <- distr_over_time[(t_max-period):t_max, ]
  
  var_vec <- NULL
  for (i in 1:nrow(results)){
    var_vec <- c(var_vec, var(results[i,] * gamete_values))
  }

  return(list(distr_over_time, mean_fitness_over_time, opti_G, gamete_values))
}


simulation <- one_simul(amplitude = amplitude, period = period,
          stoch_magnitude = stoch_magnitude, strength_selec = strength_selec,
          mut_rate = mut_rate, num_loci =  num_loci, t_max =  t_max)



# 
plot(simulation[[3]], type = 'l')
lines(simulation[[2]], col = 'blue')

plot(simulation[[3]], type = 'l')
lines(rowSums(t(t(simulation[[1]])*simulation[[4]]*2)), col = 'red')

# Simulations for the 4000 generations over the t_max period


# Just put this to know when a simulation is over if needed
# system("xdg-open 'https://www.youtube.com/watch?v=0jgrCKhxE1s'")
