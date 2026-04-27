library(doParallel)
library(foreach)

D <- 4; n1 <- 30; n2 <- 30; M <- 5; exp_size <- 100

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

p_values <- foreach(i = 1:exp_size, .combine = c) %dopar% {
  source('src/synthesis_data_generation.R')
  source('src/two_samples_test.R')
  theta <- generate_theta(D)
  df    <- generate_dataset_H0(theta, D, n1, n2, M)
  likelihood_ratio_test(df, D, n1, n2)
}

stopCluster(cl)

par(mar = c(4, 4, 2, 1))
plot(ecdf(p_values),
     main = "ECDF des p-valeurs sous H0",
     xlab = "p-value",
     ylab = "F(x)",
     col  = "steelblue")
abline(a = 0, b = 1, col = "red", lty = 2)  # référence uniforme

# Niveau empirique à 5%
mean(p_values < 0.05)  # doit être ≈ 0.05


