source('synthesis_data_generation.R')
# library(SemiMarkov)

D <- 4
n1 <- 30
n2 <- 30
M <- 5

df <- generate_dataset_H0(D, n1, n2, M)
mtrans <- matrix("W", nrow=D, ncol=D)
diag(mtrans) <- FALSE

