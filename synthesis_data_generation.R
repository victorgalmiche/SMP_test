# Number of states
D <- 4 

# Reproducibilité
set.seed(42)

# Vecteur alpha des proba initiales
alpha <- runif(D)
alpha <- alpha/sum(alpha)

# Matrice de transition 
P <- matrix(runif(D^2), nrow=D)
diag(P) <- 0
P <- P/rowSums(P)


# Génération de plusieurs chaînes de Markov 
M <- 5 # Taille  des chaînes
N <- 500 # Nombre de chaînes
states <- matrix(0, nrow=N, ncol=M)
states[,1] <- sample.int(D, size=N, replace=TRUE, prob=alpha)
for (j in 2:M){
  current_states <- states[, j-1]
  probs <- P[current_states,]
  states[, j] <- apply(probs, 1, function(p) sample.int(D, size=1, prob=p))
} 



# Paramètres des gamma - temps de séjour
a <- rexp(D) # Shape
lambda <- rexp(D) # Rate
omega <- cbind(a, lambda)

# Temps de séjour
sojournTimes <- apply(states, c(1,2), function(s) rgamma(1, shape=a[s], rate=lambda[s]))


