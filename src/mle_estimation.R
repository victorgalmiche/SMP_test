# We dispose of a dataframe respecting the semiMarkov library canvas
### LOG-LIKELIHOOD (DECOMPOSEE) ###
log_likelihood_P <- function(df, P){
  sum(log(P[cbind(df$state.h, df$state.j)]))
}

log_likelihood_omega <- function(df, omega){
  sum(dweibull(df$time, 
               shape=omega[df$state.h, 'eta'], 
               scale=omega[df$state.h, 'beta'],
               log=TRUE))
}

### MLE pour P matrice de transition ###
mle_P <- function(df, D) {
  # Count observed transitions i -> j
  counts <- matrix(0, nrow = D, ncol = D)
  idx <- cbind(df$state.h, df$state.j)
  counts <- tabulate(
    (df$state.h - 1) * D + df$state.j, 
    nbins = D * D
  ) |> matrix(nrow = D, ncol = D)
  
  # Force diagonal to 0 (semi-Markov: no self-transitions)
  diag(counts) <- 0
  
  # Normalize each row by off-diagonal sum (closed-form MLE under constraint)
  row_sums <- rowSums(counts)
  row_sums[row_sums == 0] <- 1  # avoid NaN for unvisited states
  counts / row_sums
}

#### ENCODAGE DES PARAMETRES ###
decode_omega <- function(params, D) {
  matrix_params <- matrix(params, nrow = D)
  matrix_params <- pmax(pmin(matrix_params, 10), -10)  # clamp dans [-10, 10]
  omega <- exp(matrix_params)
  colnames(omega) <- c('eta', 'beta')
  return(omega)
}



neg_ll_params_omega <- function(params_omega, df, D) {
  omega <- decode_omega(params_omega, D)
  if (any(omega <= 0) || any(!is.finite(omega))) return(1e10)
  ll <- log_likelihood_omega(df, omega)
  if (!is.finite(ll)) return(1e10)
  -ll
}

### SMARTER INIT FOR OMEGA (avoids uniroot, uses moment approximation) ###
init_params_omega <- function(df, D) {
  params <- numeric(2 * D)
  for (s in 1:D) {
    times <- df$time[df$state.h == s]
    if (length(times) < 2) next
    m  <- mean(times)
    cv <- sd(times) / m
    # Closed-form approximation: eta ≈ (cv)^{-1.086} (Menon 1963 approximation)
    eta_init  <- max(0.1, cv^(-1.086))
    beta_init <- max(1e-6, m / gamma(1 + 1/eta_init))
    params[s]     <- log(eta_init)
    params[s + D] <- log(beta_init)
  }
  params
}



### MAXIMUM LIKELIHOOD ESTIMATION ###
mle_fit <- function(df, D){
  init_states  <- tapply(df$state.h, df$id, head, 1) # Extraction du premier état de chaque chaîne
  alpha_hat <- as.vector(table(init_states)) # Nombre pour chaque état
  alpha_hat <- alpha_hat/sum(alpha_hat)
  
  P_hat <- mle_P(df, D)
  ll_P <- log_likelihood_P(df, P_hat)
  
  res_omega <- optim(
    par = init_params_omega(df, D), 
    fn = neg_ll_params_omega, 
    df = df, 
    D = D, 
    method = "Nelder-Mead", 
    control = list(maxit=5000)
  )
  
  omega_hat <- decode_omega(res_omega$par, D)
  
  theta_hat <- list(alpha=alpha_hat, P=P_hat, omega=omega_hat)
  ll <- sum(log(alpha_hat[init_states])) + ll_P - res_omega$value
  
  return(list(estimator=theta_hat, log_likelihood=ll))
}