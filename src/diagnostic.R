source('src/synthesis_data_generation.R')
source('src/two_samples_test.R')

# Generate one test dataset
set.seed(123)
D <- 4
n1 <- 30
n2 <- 30
M <- 5

theta <- generate_theta(D)
cat("Generated theta:\n")
cat("alpha:", theta$alpha, "\n")
cat("P matrix:\n")
print(theta$P)
cat("omega:\n")
print(theta$omega)

df <- generate_dataset_H0(theta, D, n1, n2, M)

cat("\nDataset summary:\n")
cat("Rows:", nrow(df), "\n")
cat("Time range:", range(df$time), "\n")
cat("Any invalid times?", any(!is.finite(df$time) | df$time <= 0), "\n")

# Fit estimates
cat("\n=== GLOBAL FIT ===\n")
global_est <- mle_fit(df, D)
cat("Global log-likelihood:", global_est$log_likelihood, "\n")
cat("Is finite?", is.finite(global_est$log_likelihood), "\n")
cat("omega:\n")
print(global_est$estimator$omega)
cat("P matrix:\n")
print(global_est$estimator$P)

# Split and fit
df1 <- subset(df, id <= n1)
df2 <- subset(df, id > n1)

cat("\n=== SAMPLE 1 FIT ===\n")
est1 <- mle_fit(df1, D)
cat("Sample 1 log-likelihood:", est1$log_likelihood, "\n")
cat("Is finite?", is.finite(est1$log_likelihood), "\n")

cat("\n=== SAMPLE 2 FIT ===\n")
est2 <- mle_fit(df2, D)
cat("Sample 2 log-likelihood:", est2$log_likelihood, "\n")
cat("Is finite?", is.finite(est2$log_likelihood), "\n")

# Compute lambda
lambda <- 2*(est1$log_likelihood + est2$log_likelihood - global_est$log_likelihood)
cat("\n=== LAMBDA ===\n")
cat("est1$ll + est2$ll:", est1$log_likelihood + est2$log_likelihood, "\n")
cat("global_est$ll:", global_est$log_likelihood, "\n")
cat("lambda = 2*(est1$ll + est2$ll - global$ll):", lambda, "\n")
cat("Is finite?", is.finite(lambda), "\n")
cat("Is negative?", lambda < 0, "\n")

# Compute p-value
dof <- D^2 + D - 1
p_val <- 1 - pchisq(lambda, df = dof)
cat("p-value:", p_val, "\n")
