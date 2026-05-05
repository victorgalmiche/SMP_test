library(doParallel)
library(foreach)

D <- 10; n1 <- 30; n2 <- 30; M <- 5; nb_datasets <- 500

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

results <- foreach(i = 1:nb_datasets, .combine = rbind, .errorhandling = 'remove') %dopar% {
  source('src/synthesis_data_generation.R')
  source('src/two_samples_test.R')
  
  theta <- generate_theta(D)
  df    <- generate_dataset_H0(theta, D, n1, n2, M)
  
  p_asymp  <- likelihood_ratio_test(df, D, n1, n2)
  p_permutation <- permutation_test(df, D, n1, n2)
  
  c(p_asymp=p_asymp, p_permutation=p_permutation)
}

stopCluster(cl)

p_asymp  <- results[, "p_asymp"]
p_permutation <- results[, "p_permutation"]

# Plot
par(mar = c(4, 4, 2, 1))
plot(ecdf(unlist(p_asymp)),
     main = "n1=n2=30 and D=10",
     xlab = "p-value",
     ylab = "Cumulative probability",
     col  = "red",
     do.points=FALSE)
lines(ecdf(p_permutation), col = "green", do.points=FALSE)
abline(a = 0, b = 1, col = "black")
legend("bottomright",
       legend = c("Chi^2", "Permutation", "Uniforme"),
       col    = c("red", "green", "black"),
       lty    = c(1, 1, 1))

