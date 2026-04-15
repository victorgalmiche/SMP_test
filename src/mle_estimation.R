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

#### ENCODAGE DES PARAMETRES ###
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

decode_omega <- function(params, D) {
  matrix_params <- matrix(params, nrow = D)
  matrix_params <- pmax(pmin(matrix_params, 10), -10)  # clamp dans [-10, 10]
  omega <- exp(matrix_params)
  colnames(omega) <- c('eta', 'beta')
  return(omega)
}



# theta → vecteur de paramètres (encodage, pour initialisation)
encode_theta <- function(theta, D) {
  c(
    log(theta$alpha),          # softmax inverse ≈ log
    as.vector(log(theta$P + 1e-10)),  # idem
    as.vector(log(theta$omega))
  )
}

neg_ll_params_P <- function(params_P, df, D) {
  P <- decode_P(params_P, D)
  ll <- log_likelihood_P(df, P)
  if (!is.finite(ll)) return(1e10)
  -ll
}

neg_ll_params_omega <- function(params_omega, df, D) {
  omega <- decode_omega(params_omega, D)
  
  # Weibull exige eta > 0 et beta > 0 (garanti par exp, mais vérifier les extrêmes)
  if (any(omega <= 0) || any(!is.finite(omega))) return(1e10)
  
  ll <- log_likelihood_omega(df, omega)
  
  if (!is.finite(ll)) return(1e10)
  -ll
}

# Ajout d'une fonction pour initialiser les paramètres
init_params_omega <- function(df, D) {
  params <- numeric(2 * D)
  for (s in 1:D) {
    times <- df$time[df$state.h == s]
    if (length(times) < 2) next  # garder 0 si pas assez de données
    
    m  <- mean(times)
    cv <- sd(times) / m  # coefficient de variation
    
    # Approximation Weibull par moments
    # cv ≈ sqrt(gamma(1+2/eta)/gamma(1+1/eta)^2 - 1)
    # On résout numériquement
    eta_init <- tryCatch({
      uniroot(function(eta) {
        sqrt(gamma(1 + 2/eta) / gamma(1 + 1/eta)^2 - 1) - cv
      }, interval = c(0.1, 20))$root
    }, error = function(e) 1.0)  # fallback
    
    beta_init <- m / gamma(1 + 1/eta_init)
    
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
  
  res_P <- optim(
    par     = rep(0, D*(D-1)),
    fn      = neg_ll_params_P,
    df      = df,
    D       = D,
    method  = "Nelder-Mead",
    control = list(maxit=5000)
  )
  
  P_hat <- decode_P(res_P$par, D)
  
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
  ll <- sum(log(alpha_hat[init_states])) - res_P$value - res_omega$value
  
  return(list(estimator=theta_hat, log_likelihood=ll))
}