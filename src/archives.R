

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

