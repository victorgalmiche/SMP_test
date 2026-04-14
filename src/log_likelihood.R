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

decode_omega <- function(params, D){
  # params vecteur de longeur Dx2 mais avec log(omega)
  matrix_params <- matrix(params, nrow=D)
  omega <- exp(params)
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
  - log_likelihood_P(df, P)
}

neg_ll_params_omega <- function(params_omega, df, D) {
  omega <- decode_omega(params_omega, D)
  - log_likelihood_omega(df, omega)
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
    par = rep(0, 2*D), 
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