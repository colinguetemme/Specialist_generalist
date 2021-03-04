library("rstudioapi")
library("ggplot2")
setwd(dirname(getActiveDocumentContext()$path))  


stats <- read.csv("stats/a50_L50_d0_s5_u0_n2.csv", sep = "\t", header = TRUE)
distr <-  read.csv("full_distr/a50_L50_d0_s5_u0_n2.csv", sep = "\t", header = TRUE)

num_loci <- 2
num_gamete <- 2^num_loci
L <- 50



distr <- distr[(1+(0*((num_gamete+1)*L))):(3*((num_gamete+1)*L)), ]
ggplot(distr, aes(x = time, y = prop*gamete_value, color = gamete_value, group = gamete_num)) + geom_line() + facet_grid(.~sim_num)
distr <-  read.csv("full_distr/a50_L50_d0_s5_u0_n2.csv", sep = "\t", header = TRUE)

distr$prop[distr$gamete_num >= 0] <- distr$prop[distr$gamete_num >= 0]*num_gamete
ggplot(distr, aes(x = time, y = prop, color = gamete_value, group = gamete_num)) + geom_line() + facet_grid(.~sim_num)
ggplot(stats, aes(x = stat, y = value)) + geom_boxplot()
