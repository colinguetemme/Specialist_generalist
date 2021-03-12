# define one gene
# suppose only uniform distribution
# and completely random model
one_gene <- function(){
  low_lim <- runif(1)
  hi_lim <- runif(1)
  if (hi_lim < low_lim){
    lim <- low_lim
    low_lim <- hi_lim
    hi_lim <- lim
  }
  value <- runif(1)
  return(c(low_lim, hi_lim, value))
}




# plotting this allele
x <- seq(0, 1, 0.0001)
plot(x, allele_value * (low_lim < x) * (x < hi_lim), ylim = c(0, 1))

# Let's define multiple loci
num_loci <- 10000
first_row <- one_gene()
mult_genes <- data.frame("low_lim" = first_row[1], "hi_lim" = first_row[2],
                         "value" = first_row[3])
for (n in 1:(num_loci-1)){
  mult_genes <- rbind(mult_genes, one_gene())
}

a <- 0;
for (fit in x){
  a <- c(a, sum(mult_genes$value[(mult_genes$low_lim < fit) * (mult_genes$hi_lim > fit)]))
}
# Suppose a time changing env
plot(a)

