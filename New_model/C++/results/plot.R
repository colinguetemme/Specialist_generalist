library("rstudioapi")
library("ggplot2")
setwd(dirname(getActiveDocumentContext()$path))  

graphics.off()
par("mar")
par(mar=c(1,1,1,1))

data_ind <- as.matrix(read.csv("data/ind.csv", header = TRUE, check.names=FALSE))
data_env <- read.csv("data/env.csv", header = F)

fig <- plot_ly(z = ~data_ind)
fig <- fig %>% add_surface()
fig


p <- plot_ly(z = data_ind, type = "surface")
