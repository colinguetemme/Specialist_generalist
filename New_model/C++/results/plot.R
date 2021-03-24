library("rstudioapi")
library("ggplot2")
library("plotly")

setwd(dirname(getActiveDocumentContext()$path))  

graphics.off()
par("mar")
par(mar=c(1,1,1,1))

data_ind <- as.matrix(read.csv("data/ind.csv", header = TRUE, check.names=FALSE))
data_env <- read.csv("data/env.csv", header = F)

# //!//
data_ind <- data_ind[, -ncol(data_ind)]
var_ind <- apply(data_ind, 1, var)

plot(var_ind)

fig <- plot_ly(z = ~data_ind)
fig <- fig %>% add_surface()
fig

p <- plot_ly(z = data_ind, type = "surface")
