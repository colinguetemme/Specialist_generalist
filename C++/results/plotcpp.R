library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))  

data <- read.csv("n4s10.csv", sep = "\t", header = FALSE)

plot(data$V1, data$V2, type = 'l')
lines(data$V3 * 3, col = "red")
lines(data$V16 * 3, col = "blue")
lines(data$V9 * 3, col = "green")
