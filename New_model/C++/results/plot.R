#' @title Plot the results of the model of Justin
#' 
#' Will plot the variance of the niche

library("rstudioapi")
library("ggplot2")
library("plotly") # might get rid of it if focus on variance

# Get the current directory
setwd(dirname(getActiveDocumentContext()$path))  

# Resize the margin not mandatory
graphics.off()
par("mar")
par(mar=c(2,2,2,2))

# get the data
## In row, the timestep, in col the genotypic value
data_ind <- as.matrix(read.csv("data/ind.csv", header = TRUE, check.names=FALSE))
data_env <- read.csv("data/env.csv", header = F)

# //!// Do not know why one extra col with NA
data_ind <- data_ind[, -ncol(data_ind)]

# Just get the variance of the niche (fitness per parameter value)
# of the population
var_ind <- apply(data_ind, 1, var)

 
mean_histo <- NULL
for (i in 1:nrow(data_ind)){
  histo <- data_ind[i, ] * 1:ncol(data_ind)
  mean_histo <- c(mean_histo, sum(histo)/sum(data_ind[i, ]))
}



fitness_val <- mean_ind[]
plot(var_ind, type = 'l')

plot(data_env[seq(1, 100, 5), 2], type = "l", ylim = c(0,1))
lines(mean_histo/1000, col = "red")
 
# Get the 3D plot of the niche changes through time
fig <- plot_ly(z = ~data_ind)
fig <- fig %>% add_surface()
fig

