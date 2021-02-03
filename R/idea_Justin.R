genetic_set <- function(number_of_genes){
  mean_range <- runif(number_of_genes, -0.9, 0.9)
  a <- rnorm(number_of_genes, 0, 0.2)
  range <- data.frame(index = 1:number_of_genes, min_range = mean_range-(a/2),
                      max_range = mean_range+(a/2), gen_value = runif(number_of_genes, 0, 1))
}

num_gen <- 10
chou <- genetic_set(num_gen)
plot(NULL, xlim = c(-1, 1), ylim = c(0, 1))
for (i in 1:num_gen){
  lines(c(chou$min_range[i], chou$max_range[i]), c(chou$gen_value[i], chou$gen_value[i]))
}

