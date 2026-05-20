source('src/synthesis_data_generation.R')
source('src/mle_estimation.R')

# Classical LR test (also called chi2)
likelihood_ratio_test <- function(df1, df2, D, law_sojourn='gamma'){
  # Combining the two groups
  df <- rbind(df1, df2)
  
  # Computing the different estimation
  global_est <- mle_fit(df, D, law_sojourn)
  est1 <- mle_fit(df1, D, law_sojourn)
  est2 <- mle_fit(df2, D, law_sojourn)
  
  # And the lambda statistic
  lambda <- 2*(est1$log_likelihood + 
                 est2$log_likelihood - 
                 global_est$log_likelihood)
  
  # if (!is.finite(lambda) || lambda < 0) return(NA)
  
  # The number of degrees of freedom depends on D
  dof <- switch(law_sojourn, 
                gamma=D^2+D-1,
                weibull=D^2+D-1,
                exponential=D^2-1)
  p_asymp <- 1-pchisq(lambda, df=dof)
  p_asymp
}

# Parametric bootstrap test
parametric_bootstrap <- function(df1, df2, D, law_sojourn='gamma', R=100) {
  df <- rbind(df1, df2)
  
  global_est <- mle_fit(df, D, law_sojourn)
  est1 <- mle_fit(df1, D, law_sojourn)
  est2 <- mle_fit(df2, D, law_sojourn)
  
  Tl <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  theta_hat <- global_est$estimator
  
  T_star <- numeric(R)
  for (r in 1:R){
    df_bootstrap <- generate_dataset_H0(theta_hat, law_sojourn, n1, n2, M)
    df1_bootstrap <- subset(df_bootstrap, id<=n1)
    df2_bootstrap <- subset(df_bootstrap, id>n1)
    
    global_est <- mle_fit(df_bootstrap, D, law_sojour)
    est1 <- mle_fit(df1_bootstrap, D, law_sojourn)
    est2 <- mle_fit(df2_bootstrap, D, law_sojourn)
    
    T_star[r] <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  }
  p_boot <- mean(T_star<=Tl)
  p_boot
}

# Permutation test
permutation_test <- function(df1, df2, D, law_sojourn='gamma', R=100) {
  df <- rbind(df1, df2)
  n1 <- length(unique(df1$id))
  n2 <- length(unique(df2$id))
  
  global_est <- mle_fit(df, D, law_sojourn)
  est1 <- mle_fit(df1, D, law_sojourn)
  est2 <- mle_fit(df2, D, law_sojourn)
  
  Tl <- global_est$log_likelihood - est1$log_likelihood - est2$log_likelihood
  T_star <- numeric(R)
  
  for (r in 1:R) {
    sample1_id   <- sample(n1 + n2, n1)
    df1_permuted <- subset(df, id %in% sample1_id)
    df2_permuted <- subset(df, !(id %in% sample1_id))
    
    est1 <- mle_fit(df1_permuted, D, law_sojourn)
    est2 <- mle_fit(df2_permuted, D, law_sojourn)
    
    T_star[r] <- global_est$log_likelihood - est1$log_likelihood - est2$log_likelihood
  }
  
  mean(T_star <= Tl, na.rm=TRUE)
}
