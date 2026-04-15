source('src/synthesis_data_generation.R')
source('src/mle_estimation.R')

D <- 7
n1 <- 30
n2 <- 30
M <- 5

likelihood_ratio_test <- function(df, D, n1, n2){
  df1 <- subset(df, id<=n1)
  df2 <- subset(df, id>n1)
  
  global_est <- mle_fit(df, D)
  est1 <- mle_fit(df1, D)
  est2 <- mle_fit(df2, D)
  
  lambda <- 2*(est1$log_likelihood + 
                 est2$log_likelihood - 
                 global_est$log_likelihood)
  
  # if (!is.finite(lambda) || lambda < 0) return(NA)
  
  dof <- D^2+D-1
  p_asymp <- 1-pchisq(lambda, df=dof)
  p_asymp
}

exp_size <- 100

# Relancer en capturant les messages d'erreur
errors <- character(exp_size)
p_values <- numeric(exp_size)

for (i in 1:exp_size) {
  repeat {
    dataset <- generate_dataset_H0(D, n1, n2, M)
    df <- dataset$data
    
    result <- tryCatch(
      likelihood_ratio_test(df, D, n1, n2),
      error = function(e) NA
    )
    
    if (!is.na(result)) break  # recommencer si échec
  }
  p_values[i] <- result
}

# Voir les messages d'erreur distincts
table(errors[errors != ""])

par(mar = c(4, 4, 2, 1))
plot(ecdf(p_values),
     main = "ECDF des p-valeurs sous H0",
     xlab = "p-value",
     ylab = "F(x)",
     col  = "steelblue")
abline(a = 0, b = 1, col = "red", lty = 2)  # référence uniforme

# Niveau empirique à 5%
mean(p_values < 0.05)  # doit être ≈ 0.05


parametric_bootstrap <- function(df, D, n1, n2, R=100) {
  df1 <- subset(df, id<=n1)
  df2 <- subset(df, id>n1)
  
  global_est <- mle_fit(df, D)
  est1 <- mle_fit(df1, D)
  est2 <- mle_fit(df2, D)
  
  Tl <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  
  for (r in 1:R){
    
  }
}

