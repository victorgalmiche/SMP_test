
# Fonction permettant de générer les paramètres theta définissant un SMP
generate_theta <- function(D){
  # proba initiales
  alpha <- runif(D)
  alpha <- alpha/sum(alpha)
  
  # matrice de transition 
  P <- matrix(runif(D^2), nrow=D)
  diag(P) <- 0
  P <- P/rowSums(P)
  
  # temps de séjour - Weibull distribution
  # eta <- rexp(D, rate=0.8) # shape 
  # beta <- rexp(D, rate=0.15) # scale 
  # omega <- cbind(eta, beta)
  
  # temps de séjour - Gamma distribution 
  a <- rexp(D, rate=0.5) # shape
  lambda <- rexp(D, rate=3) # rate = 1/scale
  omega <- cbind(a, lambda)
  
  return (list(alpha=alpha, P=P, omega=omega))
}

# Fonction permettant de générer un dataset
generate_dataset_H0 <- function(theta, D, n1, n2, M){
  n <- n1 + n2
  
  # Génération chaînes
  states <- matrix(0, nrow=n, ncol=M+1) # On ajoute une colonne de manière artificielle
  states[,1] <- sample.int(D, size=n, replace=TRUE, prob=theta$alpha)
  for (j in 2:(M+1)){
    current_states <- states[, j-1]
    probs <- theta$P[current_states,]
    states[, j] <- apply(probs, 1, function(p) sample.int(D, size=1, prob=p))
  } 
  
  # Temps de séjour
  # sojournTimes <- apply(states, c(1,2), function(s) rweibull(1, shape=theta$omega[s, 'eta'], scale=theta$omega[s, 'beta']))
  sojournTimes <- apply(states, c(1,2), function(s) {
    max(rgamma(1, shape=theta$omega[s, 'a'], rate=theta$omega[s, 'lambda']),
        .Machine$double.xmin)
  })

  
  # DataFrame respectant le format de SemiMarkov
  id <- rep(1:n, each=M)
  state.h <- as.vector(t(states[, 1:M]))
  state.j <- as.vector(t(states[, 2:(M+1)]))
  time <- as.vector(t(sojournTimes[, 1:M]))
  df <- data.frame(id=id, state.h=state.h, state.j=state.j, time=time)
  
  return(df)
} 

