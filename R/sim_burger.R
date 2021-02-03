# Modelling of the paper of Burger 2002 
# study of a range of temporal variation 

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

quick_params_selection <- 2

if (quick_params_selection == 1){
  amplitude <- 0.5
  period <- 24
  stoch_magnitude <- 0
  strength_selec <- 1
  mut_rate <- 0
  num_loci <- 2
  t_max <- 50
}

if (quick_params_selection == 2){
  amplitude <- 0.5
  period <- 50
  stoch_magnitude <- 0.2
  strength_selec <- 20
  mut_rate <- 5*10^(-5)
  num_loci <- 4
  t_max <- 100
}

if (quick_params_selection == 3){
  amplitude <- 0.5
  period <- 100
  stoch_magnitude <- 0
  strength_selec <- 5
  mut_rate <- 5*10^(-5)
  num_loci <- 6
  t_max <- 200
}


# One simulation over the t_max period

one_simul <- function(amplitude, period, stoch_magnitude, strength_selec,
                      mut_rate, num_loci, t_max){
  
  initialisation <- genetic_sim(num_loci = num_loci)
  gamete_values <- initialisation[[1]]
  gamete_distr <- initialisation[[2]]
  gamete_scheme <- initialisation [[3]]
  rec_rate <- initialisation[[4]]
  
  dtf_rec <- recombination(num_loci = num_loci, rec_rate = rec_rate)

  opti_G <- get_optimum(t_max = t_max, amplitude = amplitude, period = period, 
                        stoch_magnitude = stoch_magnitude)
  
  distr_over_time <- matrix(0, ncol = length(gamete_distr), nrow = t_max+1)
  distr_over_time[1, ] <- gamete_distr
  mean_fitness_over_time <- rep(0, t_max)
  for (t in 1:t_max){
    opti_Gt <- opti_G[t]
    new_results <- new_distributions(num_loci, distr_over_time[t, ], gamete_values,
                                       gamete_scheme, opti_Gt, dtf_rec, mut_rate)
    distr_over_time[t+1, ] <- new_results[[1]]
    mean_fitness_over_time[t] <- new_results[[2]]
  }
  return(list(distr_over_time[-1,], mean_fitness_over_time, opti_G, gamete_values))
}


simulation <- one_simul(amplitude = amplitude, period = period,
          stoch_magnitude = stoch_magnitude, strength_selec = strength_selec,
          mut_rate = mut_rate, num_loci =  num_loci, t_max =  t_max)


plot(simulation[[3]], type = 'l')
lines(rowSums(t(t(simulation[[1]])*simulation[[4]]*2)), col = 'red')

plot(simulation[[3]], type = 'l')
lines(simulation[[2]], col = 'blue')


# Just put this to know when a simulation is over if needed
# system("xdg-open 'https://www.youtube.com/watch?v=0jgrCKhxE1s'")
