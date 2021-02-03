# Modelling of the paper of Burger 2002 
# study of a range of temporal variation 

# Functions
# To get a binary vector
# Number is a positive integer that will be transform in binary
# n is the number of bit in the vector

# > num_to_bin(23, 4)
# [1] 0 1 1 1

num_to_bin <- function(number, n) {
  bin_vector <- rev(as.numeric(intToBits(number)))
  if(missing(n)) {
    return(bin_vector)
  } else {
    bin_vector[-(1:(length(bin_vector) - n))]
  }
}


# initialisation modelling

# will return  a list with respectively the gamete phenotypic values, the
# initial distribution of gametes, the actual allelic composition of each 
# gamete and the recombination rate between each adjacent pair of loci

genetic_sim <- function(num_loci, haploid = F){
  
  nb_gamete <- 2^num_loci
  
  rec_rate <- runif(num_loci-1, 0, 0.5)
  ifelse(haploid, haploid <- 1, haploid <- 0.5)
  loci <- runif(num_loci, 0, 1)
  
  # scale the loci values to get a heterozygote of phenotypic value 0.5.
  loci <- haploid * loci / sum(loci)
  
  # initialize the gamete values and scheme (allelic composition)
  gamete_values <- rep(0, nb_gamete)
  gamete_scheme <- list(rep(0, nb_gamete))
  
  for (i in 1:(nb_gamete)){
    gamete_scheme[[i]] <- num_to_bin(i-1, n = num_loci)
    gamete_values[i] <- sum(loci[as.logical(num_to_bin(i-1, n = num_loci))])
  }
  gamete_distr <- runif(nb_gamete, 0, 1)
  gamete_distr <- gamete_distr / sum(gamete_distr)
  
  return(list(gamete_values, gamete_distr, gamete_scheme, rec_rate))
}




# Calculation of the optimal phenotype over time

get_optimum <- function(t_max, amplitude, period, stoch_magnitude){
  opti_G <- 1:t_max
  opti_G <- lapply(opti_G,
                   function(t) {rnorm(1, 0.5 + amplitude * sin(2*pi*t/period),
                                              stoch_magnitude * amplitude)})
  return(unlist(opti_G))
}


### recombination probabilities table ###

# will transform two gametes considering a list of recombination position
permut_gamete <- function(j_bin, k_bin, rec_vec){
  
  # Consider the first value of the vector of recombination position
  rec_spot <- rec_vec[1]

  j_bin_rec <- c(j_bin[1:rec_spot], k_bin[(rec_spot+1):length(j_bin)])
  k_bin_rec <- c(k_bin[1:rec_spot], j_bin[(rec_spot+1):length(k_bin)])
  
  if (length(rec_vec) > 1) {
    # Recall the function for the next value of the rec_vec
    permut_gamete(j_bin_rec, k_bin_rec, rec_vec[-1])
  } else {
    return(list(j_bin_rec, k_bin_rec))
  }
}

jk_i_recombination <- function(i_bin, j_bin, k_bin, rec_rate){
  R_jk_to_i <- 0
  n <- length(rec_rate)
  l <- rep(list(0:1), n)
  l <- expand.grid(l)
  l <- t((1:(num_loci-1))*t(l))
  for (t in 1:nrow(l)){
    
    rec_vec = l[t, ]

    if (sum(rec_vec)>0){
      rec_vec <- unlist(lapply(rec_vec, function(x) {x[x!=0]}))
      jk <- permut_gamete(j_bin, k_bin, rec_vec)

    } else {
      jk <- list(j_bin, k_bin)
      rec_vec <- -(1:(num_loci-1))
    }
    if (all(jk[[1]] == i_bin)){
      R_jk_to_i <- R_jk_to_i + 0.5 * prod(rec_rate[rec_vec]) * prod(1-rec_rate[-rec_vec])
    }
    if (all(jk[[2]] == i_bin)){
      R_jk_to_i <- R_jk_to_i + 0.5 * prod(rec_rate[rec_vec]) * prod(1-rec_rate[-rec_vec])
    }
  } 

  return(R_jk_to_i)
}




# Recombination of every cbn

recombination <- function(num_loci, rec_rate){
  
  prob_rec <- matrix(0, ncol = 2^num_loci, nrow = 2^(num_loci*2))
  index <- 0
  nb_gamete <- 2^num_loci
  names_cmb <- NULL
  
  for (j in 1:nb_gamete - 1){
    
    j_bin <- num_to_bin(j, num_loci)
    
    for (k in j:(nb_gamete - 1)){
      
      k_bin <- num_to_bin(k, num_loci)
      names_cmb = c(names_cmb, paste(j+1,',',k+1))
      index <- index + 1
      
      for (i in 1:nb_gamete-1){
        
        i_bin <- num_to_bin(i, num_loci)
        prob_rec[index, i+1] <- jk_i_recombination(i_bin, j_bin, k_bin,
                                                   rec_rate)
      }
    }
  }
  
  prob_rec <- prob_rec[1:length(names_cmb), ]
  dtf_rec <- data.frame(t(prob_rec))
  names(dtf_rec) <- names_cmb
  return(dtf_rec)
}

# Calculation of all pi_star

new_distributions <- function(num_loci, gamete_distr, gamete_values,
                              gamete_scheme, opti_Gt, dtf_rec, mut_rate){
  nb_gamete <- 2^num_loci
  p_star <- rep(0, nb_gamete)

  mean_fitness <- 0
  cool = NULL
  for (j in 1:nb_gamete){
    for (k in j:nb_gamete){

      mean_fitness <- mean_fitness + exp(-strength_selec *
                      (gamete_values[j] + gamete_values[k] - opti_Gt)^2) *
                      gamete_distr[j] * gamete_distr[k]
      cool <- c(cool, gamete_distr[j] * gamete_distr[k])
    }
  }

  for (i in 1:nb_gamete){
    for (j in 1:nb_gamete){
      for (k in j:nb_gamete){
        
        G <- gamete_values[j] + gamete_values[k]
        viability_G <- (exp(-strength_selec * (G - opti_Gt)^2) * gamete_distr[j] *
                          gamete_distr[k] * dtf_rec[i, paste(j,',',k)]) / mean_fitness
        p_star[i] <- p_star[i] + viability_G
      } 
    }
  }

  for (i in 1:nb_gamete) {
    
    tot_mut <- 0
    
    for (j in (1:nb_gamete)[-i]){
      mutation_ij <- mut_rate ^ (sum(gamete_scheme[[i]] != gamete_scheme[[j]]))
      tot_mut <- tot_mut + p_star[j]*mutation_ij - p_star[i]*mutation_ij
    }
    
    gamete_distr[i] <- p_star[i] + tot_mut
  }
  return(list(gamete_distr, mean_fitness))
}

new_distributions_haploid <- function(num_loci, gamete_distr, gamete_values,
                          gamete_scheme, opti_Gt, mut_rate){
  nb_gamete <- 2^num_loci
  p_star <- rep(0, nb_gamete)
  

  for (i in 1:nb_gamete){

      gamete_distr[i] <- exp(-strength_selec * (gamete_values[i] - opti_Gt)^2) * gamete_distr[i]

  }

  
  for (i in 1:nb_gamete) {
    
    tot_mut <- 0
    
    for (j in (1:nb_gamete)[-i]){
      mutation_ij <- mut_rate ^ (sum(gamete_scheme[[i]] != gamete_scheme[[j]]))
      tot_mut <- tot_mut + gamete_distr[j]*mutation_ij - gamete_distr[i]*mutation_ij
    }
    
    gamete_distr[i] <- gamete_distr[i] + tot_mut
  }
  gamete_distr <- gamete_distr / sum(gamete_distr)
  return(list(gamete_distr))
}


