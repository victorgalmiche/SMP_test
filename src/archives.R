

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