

softmax <- function(x) {
  e <- exp(x - max(x))  # stabilité numérique
  e / sum(e)
}


decode_P <- function(params, D){
  # params un vecteur de longueur Dx(D-1) (logits)
  matrix_params <- matrix(params, nrow=D)
  coef <- t(apply(matrix_params, 1, softmax))
  P <- matrix(0, nrow=D, ncol=D)
  for (i in 1:D){
    P[i, (1:D)[-i]] <- coef[i, ]
  }
  return(P)
}


neg_ll_params_P <- function(params_P, df, D) {
  P <- decode_P(params_P, D)
  ll <- log_likelihood_P(df, P)
  if (!is.finite(ll)) return(1e10)
  -ll
}

# MLE FOR WEIBULL

log_likelihood_omega <- function(df, omega){
  sum(dweibull(df$time,
               shape=omega[df$state.h, 'eta'],
               scale=omega[df$state.h, 'beta'],
               log=TRUE))
}

decode_omega <- function(params, D) {
  matrix_params <- matrix(params, nrow = D)
  matrix_params <- pmax(pmin(matrix_params, 10), -10)  # clamp dans [-10, 10]
  omega <- exp(matrix_params)
  colnames(omega) <- c('eta', 'beta')
  return(omega)
}

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


neg_ll_params_omega <- function(params_omega, df, D) {
  omega <- decode_omega(params_omega, D)
  if (any(omega <= 0) || any(!is.finite(omega))) return(1e10)
  ll <- log_likelihood_omega(df, omega)
  if (!is.finite(ll)) return(1e10)
  -ll
}


#### ENCODAGE DES PARAMETRES ###
decode_omega <- function(params, D) {
  matrix_params <- matrix(params, nrow = D)
  omega <- exp(matrix_params)
  colnames(omega) <- c('a', 'lambda')
  return(omega)
}



library(fitdistrplus)

mle_omega <- function(df, D) {
  # Initialize matrix of parameters
  omega <- matrix(1, nrow = D, ncol = 2,
                  dimnames = list(1:D, c("a", "lambda")))
  
  for (i in 1:D) {
    times_i <- df$time[df$state.h == i]
    
    if (length(times_i) < 2) {
      warning(paste("État", i, ": pas assez d'observations"))
      next
    }
    
    fit <- fitdist(times_i, distr = "gamma", method = "mle",
                   lower = c(shape = 1e-6, rate = 1e-6))
    
    omega[i, "a"]      <- fit$estimate["shape"]
    omega[i, "lambda"] <- fit$estimate["rate"]
  }
  omega
}

mle_gamma_closed <- function(x) {
  s <- log(mean(x)) - mean(log(x))   # toujours > 0
  # Approximation de Choi & Wette (1969) pour shape
  shape <- (3 - s + sqrt((s-3)^2 + 24*s)) / (12*s)
  # Affinage par Newton (1-2 itérations suffisent)
  for (k in 1:5) {
    shape <- shape - (log(shape) - digamma(shape) - s) /
      (1/shape - trigamma(shape))
  }
  rate <- shape / mean(x)
  c(a = shape, lambda = rate)
}

mle_omega_closed <- function(df, D) {
  t(sapply(1:D, function(i) {
    times_i <- df$time[df$state.h == i]
    if (length(times_i) < 2) return(c(a=1, lambda=1))
    mle_gamma_closed(times_i)
  }))
}

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

