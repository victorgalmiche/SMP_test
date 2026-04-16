source('src/synthesis_data_generation.R')
source('src/two_samples_test.R')

D <- 4
n1 <- 30
n2 <- 30
M <- 5

library(parallel)
cl <- makeCluster(detectCores() - 1, rscript_args = "--vanilla")
clusterExport(cl, ls())  # export all functions/variables

p_values <- parSapply(cl, 1:exp_size, function(i) {
  repeat {
    theta <- generate_theta(D)
    df    <- generate_dataset_H0(theta, D, n1, n2, M)
    result <- tryCatch(parametric_bootstrap(df, D, n1, n2), error = function(e) NA)
    if (!is.na(result)) break
  }
  result
})
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


