# We dispose of a dataframe respecting the semiMarkov library canvas
### LOG-LIKELIHOOD (DECOMPOSEE) ###
log_likelihood_P <- function(df, P){
  sum(log(P[cbind(df$state.h, df$state.j)]))
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
  ) |> matrix(nrow = D, ncol = D)
  
  # Force diagonal to 0 (semi-Markov: no self-transitions)
  diag(counts) <- 0
  
  # Normalize each row by off-diagonal sum (closed-form MLE under constraint)
  row_sums <- rowSums(counts)
  row_sums[row_sums == 0] <- 1  # avoid NaN for unvisited states
  counts / row_sums
}


### MLE pour omega ###
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


### MAXIMUM LIKELIHOOD ESTIMATION ###
mle_fit <- function(df, D){
  init_states  <- tapply(df$state.h, df$id, head, 1) # Extraction du premier état de chaque chaîne
  alpha_hat <- as.vector(table(init_states)) # Nombre pour chaque état
  alpha_hat <- alpha_hat/sum(alpha_hat)
  
  P_hat <- mle_P(df, D)
  ll_P <- log_likelihood_P(df, P_hat)
  
  # omega_hat <- mle_omega(df, D)
  omega_hat <- mle_omega_closed(df, D)
  ll_omega <- log_likelihood_omega(df, omega_hat)
  
  theta_hat <- list(alpha=alpha_hat, P=P_hat, omega=omega_hat)
  ll <- sum(log(alpha_hat[init_states])) + ll_P + ll_omega
  
  return(list(estimator=theta_hat, log_likelihood=ll))
}