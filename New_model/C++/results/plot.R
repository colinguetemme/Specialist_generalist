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
data_ind <- as.matrix(read.csv("data/ind.csv", header = TRUE, check.names=FALSE))
data_env <- read.csv("data/env.csv", header = F)

# //!// Do not know why one extra col with NA
data_ind <- data_ind[, -ncol(data_ind)]

# Just get the variance of the niche (fitness per parameter value)
# of the population
var_ind <- apply(data_ind, 1, var)
plot(var_ind, type = 'l')

# Get the 3D plot of the niche changes through time
fig <- plot_ly(z = ~data_ind)
fig <- fig %>% add_surface()
fig

p <- plot_ly(z = data_ind, type = "surface")
