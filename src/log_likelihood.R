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

decode_P <- function(params){
  # params de taille Dx(D-1) (logits)
  D <- nrow(params)
  coef <- t(apply(params, 1, softmax))
  P <- matrix(0, nrow=D, ncol=D)
  for (i in 1:D){
    P[i, (1:D)[-i]] <- coef[i, ]
  }
  return(P)
}

decode_omega <- function(params){
  # params de taille Dx2 mais avec log(omega)
  omega <- exp(params)
  colnames(omega) <- c('eta', 'beta')
  return(omega)
}

# Vecteur de paramètres -> theta (décodage)
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