source('synthesis_data_generation.R')
# library(SemiMarkov)

D <- 4
n1 <- 30
n2 <- 30
M <- 5

df <- generate_dataset_H0(D, n1, n2, M)
mtrans <- matrix("W", nrow=D, ncol=D)
diag(mtrans) <- FALSE

log_likelihood <- function(theta, df, D) {
  ll <- 0
  
  # 1. Contribution des états initiaux
  first_states <- df$state.h[df$id == unique(df$id)[1]]  # premier état de chaque chaîne
  init_states  <- tapply(df$state.h, df$id, head, 1)
  ll <- ll + sum(log(theta$alpha[init_states]))
  
  # 2. Contribution des transitions
  ll <- ll + sum(log(theta$P[cbind(df$state.h, df$state.j)]))
  
  # 3. Contribution des temps de séjour (Weibull)
  ll <- ll + sum(dweibull(
    df$time,
    shape = theta$omega[df$state.h, "eta"],
    scale = theta$omega[df$state.h, "beta"],
    log   = TRUE
  ))
  
  return(ll)
}

softmax <- function(x) {
  e <- exp(x - max(x))  # stabilité numérique
  e / sum(e)
}

# Vecteur de paramètres → theta (décodage)
decode_theta <- function(params, D) {
  
  alpha <- softmax(params[1:D])
  
  P_raw <- matrix(params[(D+1):(D+D^2)], nrow=D)
  diag(P_raw) <- -Inf
  P <- t(apply(P_raw, 1, softmax))
  
  log_omega <- matrix(params[(D+D^2+1):length(params)], ncol=2)
  omega <- exp(log_omega)
  colnames(omega) <- c("eta", "beta")
  
  list(alpha=alpha, P=P, omega=omega)
}

# theta → vecteur de paramètres (encodage, pour initialisation)
encode_theta <- function(theta, D) {
  c(
    log(theta$alpha),          # softmax inverse ≈ log
    as.vector(log(theta$P + 1e-10)),  # idem
    as.vector(log(theta$omega))
  )
}

# Negative log-likelihood (optim minimise)
neg_ll <- function(params, df, D) {
  theta <- decode_theta(params, D)
  -log_likelihood(theta, df, D)
}

# Optimisation par Nelder-Mead
fit_mle <- function(df, D, theta_init=NULL) {
  
  if (is.null(theta_init)) {
    params_init <- rep(0, D + D^2 + D*2)
  } else {
    params_init <- encode_theta(theta_init, D)
  }
  
  res <- optim(
    par     = params_init,
    fn      = neg_ll,
    df      = df,
    D       = D,
    method  = "Nelder-Mead",
    control = list(maxit=5000, reltol=1e-10)
  )
  
  list(
    theta    = decode_theta(res$par, D),  # paramètres estimés
    loglik   = -res$value,                # log-vraisemblance maximisée
    convergence = res$convergence,        # 0 = succès
    counts   = res$counts                 # nb d'évaluations
  )
}

# Utilisation
res <- fit_mle(df, D)
res$convergence  # doit être 0
res$loglik
res$theta
