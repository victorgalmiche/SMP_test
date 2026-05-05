source('src/synthesis_data_generation.R')
source('src/mle_estimation.R')

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


parametric_bootstrap <- function(df, D, n1, n2, R=100) {
  df1 <- subset(df, id<=n1)
  df2 <- subset(df, id>n1)
  
  global_est <- mle_fit(df, D)
  est1 <- mle_fit(df1, D)
  est2 <- mle_fit(df2, D)
  
  Tl <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  theta_hat <- global_est$estimator
  
  T_star <- numeric(R)
  for (r in 1:R){
    df_bootstrap <- generate_dataset_H0(theta_hat, D, n1, n2, M)
    df1_bootstrap <- subset(df_bootstrap, id<=n1)
    df2_bootstrap <- subset(df_bootstrap, id>n1)
    
    global_est <- mle_fit(df_bootstrap, D)
    est1 <- mle_fit(df1_bootstrap, D)
    est2 <- mle_fit(df2_bootstrap, D)
    
    T_star[r] <- global_est$log_likelihood -est1$log_likelihood - est2$log_likelihood
  }
  p_boot <- mean(T_star<=Tl)
  p_boot
}


permutation_test <- function(df, D, n1, n2, R=100) {
  df1 <- subset(df, id <= n1)
  df2 <- subset(df, id > n1)
  
  global_est <- mle_fit(df, D)
  est1       <- mle_fit(df1, D)
  est2       <- mle_fit(df2, D)
  
  Tl <- global_est$log_likelihood - est1$log_likelihood - est2$log_likelihood
  T_star <- numeric(R)
  
  if (is.na(Tl)) stop("Tl is NA: mle_fit failed on observed data or full data")
  
  
  for (r in 1:R) {
    sample1_id   <- sample(n1 + n2, n1)
    df1_permuted <- subset(df, id %in% sample1_id)
    df2_permuted <- subset(df, !(id %in% sample1_id))
    
    est1 <- mle_fit(df1_permuted, D)
    est2 <- mle_fit(df2_permuted, D)
    
    T_star[r] <- global_est$log_likelihood - est1$log_likelihood - est2$log_likelihood
  }
  
  mean(T_star <= Tl, na.rm=TRUE)
}
