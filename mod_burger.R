# Modelling of the paper of Burger 2002 
# study of a range of temporal variation 

# Functions
# To get a binary vector
num_to_bin <- function(number, n) {
  bin_vector <- rev(as.numeric(intToBits(number)))
  if(missing(n)) {
    return(bin_vector)
  } else {
    bin_vector[-(1:(length(bin_vector) - n))]
  }
}


# initialisation modelling

genetic_sim <- function(num_loci, haploid = F){
  rec_rate <- runif(num_loci - 1, 0, 0.5)
  ifelse(haploid, haploid <- 1, haploid <- 0.5)
  loci <- runif(num_loci, 0, 1)
  loci <- haploid * loci / sum(loci)
  
  nb_gamette <- 2^num_loci
  gamette_values <- rep(0, nb_gamette)
  gamette_scheme <- list(rep(0, nb_gamette))

  for (i in 1:(nb_gamette)){
    gamette_scheme[[i]] <- num_to_bin(i-1, n = num_loci)
    gamette_values[i] <- sum(loci[as.logical(num_to_bin(i-1, n = num_loci))])
  }
  gamette_distr <- runif(nb_gamette, 0, 1)
  gamette_distr <- gamette_distr / sum(gamette_distr)
  return(list(gamette_values, gamette_distr, gamette_scheme, rec_rate))
}




# Calculation of the optimal phenotype over time

get_optimum <- function(t_max, amplitude, period, stoch_magnitude){
  opti_G <- 1:t_max
  opti_G <- lapply(opti_G,
                   function(t) {rnorm(1, 0.5 + amplitude * sin(2*pi*t/period),
                                              stoch_magnitude * amplitude)})
  return(unlist(opti_G))
  }


# recombination

permut_gamette <- function(j_bin, k_bin, rec_vec){
  rec_spot <- rec_vec[1]

  j_bin_rec <- c(j_bin[1:rec_spot], k_bin[(rec_spot+1):length(j_bin)])
  k_bin_rec <- c(k_bin[1:rec_spot], j_bin[(rec_spot+1):length(k_bin)])
  if (length(rec_vec) > 1) {
    permut_gamette(j_bin_rec, k_bin_rec, rec_vec[-1])
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
      jk <- permut_gamette(j_bin, k_bin, rec_vec)

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
  nb_gamette <- 2^num_loci
  names_cmb <- NULL
  
  for (j in 1:nb_gamette - 1){
    
    j_bin <- num_to_bin(j, num_loci)
    
    for (k in j:(nb_gamette - 1)){
      
      k_bin <- num_to_bin(k, num_loci)
      names_cmb = c(names_cmb, paste(j+1,',',k+1))
      index <- index + 1
      
      for (i in 1:nb_gamette-1){
        
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

new_distributions <- function(num_loci, gamette_distr, gamette_values,
                              gamette_scheme, opti_Gt, dtf_rec, mut_rate){
  nb_gamette <- 2^num_loci
  p_star <- rep(0, nb_gamette)

  mean_fitness <- 0
  
  for (j in 1:nb_gamette){
    for (k in j:nb_gamette){
      
      mean_fitness <- mean_fitness + exp(-strength_selec *
                      (gamette_values[j] + gamette_values[k] - opti_Gt)^2) *
                      gamette_distr[j] * gamette_distr[k]
      
    }
  }

  for (i in 1:nb_gamette){
    for (j in 1:nb_gamette){
      for (k in j:nb_gamette){
        
        G <- gamette_values[j] + gamette_values[k]
        viability_G <- (exp(-strength_selec * (G - opti_Gt)^2) * gamette_distr[j] *
                          gamette_distr[k] * dtf_rec[i, paste(j,',',k)]) / mean_fitness
        p_star[i] <- p_star[i] + viability_G
      } 
    }
  }

  for (i in 1:nb_gamette) {
    
    tot_mut <- 0
    
    for (j in (1:nb_gamette)[-i]){
      mutation_ij <- mut_rate ^ (sum(gamette_scheme[[i]] != gamette_scheme[[j]]))
      tot_mut <- tot_mut + p_star[j]*mutation_ij - p_star[i]*mutation_ij
    }
    
    gamette_distr[i] <- p_star[i] + tot_mut
  }
  return(list(gamette_distr, mean_fitness))
}

new_distributions_haploid <- function(num_loci, gamette_distr, gamette_values,
                          gamette_scheme, opti_Gt, mut_rate){
  nb_gamette <- 2^num_loci
  p_star <- rep(0, nb_gamette)
  

  for (i in 1:nb_gamette){

      gamette_distr[i] <- exp(-strength_selec * (gamette_values[i] - opti_Gt)^2) * gamette_distr[i]

  }

  
  for (i in 1:nb_gamette) {
    
    tot_mut <- 0
    
    for (j in (1:nb_gamette)[-i]){
      mutation_ij <- mut_rate ^ (sum(gamette_scheme[[i]] != gamette_scheme[[j]]))
      tot_mut <- tot_mut + gamette_distr[j]*mutation_ij - gamette_distr[i]*mutation_ij
    }
    
    gamette_distr[i] <- gamette_distr[i] + tot_mut
  }
  gamette_distr <- gamette_distr / sum(gamette_distr)
  return(list(gamette_distr))
}


