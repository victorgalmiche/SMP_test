# We dispose of a dataframe respecting the semiMarkov library canvas
### LOG-LIKELIHOOD (DECOMPOSEE) ###
log_likelihood_P <- function(df, P){
  P_safe <- pmax(P, 1e-10) # To avoid log(0)
  sum(log(P_safe[cbind(df$state.h, df$state.j)]))
}


log_likelihood_omega <- function(df, omega){
  sum(dgamma(df$time,
             shape=omega[df$state.h, 'a'],
             rate=omega[df$state.h, 'lambda'],
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
  ) |> matrix(nrow = D, ncol = D, byrow=TRUE)
  
  # Force diagonal to 0 (semi-Markov: no self-transitions)
  diag(counts) <- 0
  
  # Normalize each row by off-diagonal sum (closed-form MLE under constraint)
  row_sums <- rowSums(counts)
  row_sums[row_sums == 0] <- Inf  # avoid NaN for unvisited states
  counts / row_sums
}


### MLE pour omega ###
mle_omega <- function(df, D){
  omega <- matrix(1, nrow=D, ncol=2, dimnames=list(1:D, c("a", "lambda")))
  for (i in 1:D){
    times_i <- df$time[df$state.h==i]
    
    if (length(times_i) < 2) {
      warning(paste('State', i, 'has less than 2 valid observations'))
      next
    }
    
    # Objective function: negative log likelihood
    nll <- function(pars){
      if (pars[1] <= 0 || pars[2] <= 0) return(1e10)
      ll <- sum(dgamma(times_i, shape=pars[1], rate=pars[2], log=TRUE))
      if (!is.finite(ll)) return(1e10)
      -ll
    }
    
    # Starting values using method of moments
    mean_x <- mean(times_i)
    var_x <- var(times_i)
    
    # Handle invalid variance
    if (var_x <= 0 || !is.finite(var_x)) {
      var_x <- (mean_x / 4)^2
    }
    
    shape_start <- mean_x^2 / var_x
    rate_start <- mean_x / var_x
    
    # Bound starting values to reasonable range
    shape_start <- pmax(0.01, pmin(shape_start, 500))
    rate_start <- pmax(0.001, pmin(rate_start, 500))
    
    # Test starting values
    if (!is.finite(nll(c(shape_start, rate_start)))) {
      warning(paste('State', i, 'invalid starting values'))
      next
    }
    
    # Optimize with L-BFGS-B (has parameter bounds)
    result <- tryCatch(
      optim(c(shape_start, rate_start), nll, 
            method='L-BFGS-B',
            lower=c(0.001, 0.001),
            upper=c(1000, 1000),
            control=list(maxit=5000)),
      error = function(e) {
        warning(paste('State', i, 'optimization error'))
        list(par=c(shape_start, rate_start), convergence=1)
      }
    )
    
    omega[i, 'a'] <- result$par[1]
    omega[i, 'lambda'] <- result$par[2]
  }
  omega
}

mle_omega_nm <- function(df, D){
  omega <- matrix(1, nrow=D, ncol=2, dimnames=list(1:D, c("a", "lambda")))
  for (i in 1:D){
    times_i <- df$time[df$state.h==i]
    
    if (length(times_i) < 2) {
      warning(paste('State', i, 'has not enough observations'))
      next
    }
    
    # Objective function: negative log likelihood
    nll <- function(pars){
      if (pars[1] <= 0 || pars[2] <= 0) return (1e10) # Penalize invalid parameters
      - sum(dgamma(times_i, shape=pars[1], rate=pars[2], log=TRUE))
    }
    
    # Starting w/ values based on the method of moments
    mean_x <- mean(times_i)
    var_x <- var(times_i)
    
    # Ensure var_x is positive and not too small
    if (var_x <= 0 || is.na(var_x)) {
      var_x <- mean_x^2 / 10  # Use fallback if variance is invalid
    }
    
    # Calculate method of moments starting values with safety checks
    shape_start <- mean_x^2 / var_x
    rate_start <- mean_x / var_x
    
    # Ensure starting values are positive and reasonable
    shape_start <- max(0.1, min(shape_start, 100))
    rate_start <- max(0.01, min(rate_start, 100))
    start <- c(shape=shape_start, rate=rate_start)
    
    # Optimization w/ Nelder-Mead
    result <- tryCatch(
      optim(start, nll, 
            method='Nelder-Mead', 
            control=list(maxit=1000, reltol=1e-8)),
      error = function(e) list(par=start, convergence=1)
    )
    
    omega[i, 'a'] <- result$par[1]
    omega[i, 'lambda'] <- result$par[2]
  }
  omega
}


### MAXIMUM LIKELIHOOD ESTIMATION ###
mle_fit <- function(df, D){
  init_states  <- tapply(df$state.h, df$id, head, 1) # Extraction du premier état de chaque chaîne
  alpha_hat <- as.vector(table(init_states)) # Nombre pour chaque état
  alpha_hat <- alpha_hat/sum(alpha_hat)
  
  P_hat <- mle_P(df, D)
  ll_P <- log_likelihood_P(df, P_hat)
  
  omega_hat <- mle_omega(df, D)
  ll_omega <- log_likelihood_omega(df, omega_hat)
  
  theta_hat <- list(alpha=alpha_hat, P=P_hat, omega=omega_hat)
  ll <- sum(log(alpha_hat[init_states])) + ll_P + ll_omega
  
  return(list(estimator=theta_hat, log_likelihood=ll))
}