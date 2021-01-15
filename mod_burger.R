
# amplitude is the size of the range of values for the optimum genotype
amplitude <- 1
# period is the length between 2 period
period <- 100
# magnitude important is the noise add to the sin curve describing the optimum 
# genotype
magnitude <- 0.2

# t is the current time
t <- 1

# strength_selec isstrenght of stabilization (noted s)
strength_selec <- 0.5
# viability of an individual with genotype G
G <- 0.3
viability_G <- exp(-strength_selec * (G - opti_G)^2)

# opti_G is the optimum genotype (noted theta)
opti_G <- rep(0, 200)
for (t in 1:200) {
  opti_G[t] <- rnorm(1, 0.5 + amplitude * sin(2*pi*t/period), magnitude * amplitude)
}
plot(opti_G, type = "l")

# mut_rate is the common mutation rate (equally to pass from gamete j to k)
mut_rate <- 10^(-6)

# num_loci the number of loci on the gamete
num_loci <- 2
loci <- runif(num_loci, 0, 1)
loci <- 0.5 * loci / sum(loci)

#rec_rate the recombination rate
rec_rate <- runif(num_loci - 1, 0, 0.5)

prop <- runif(1 + 0.5 * num_loci * (num_loci + 1), 0, 1)
prop <- prop / sum(prop)

