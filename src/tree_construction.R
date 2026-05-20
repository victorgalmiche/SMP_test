### REGRESSION TREE CONSTRUCTION
# Split criterion: p-value of the two-sample test

# Data structure:
# - trajectories: data.frame, each row is a trajectory
# - covariates: data.frame, each row are the covariate for an individual


# Conversion function from trajectories to dataframe format used by two_sample_test.R
traj_to_df <- function(trajectories) {
  res <- list()
  
  for (i in seq_len(nrow(trajectories))) {
    states <- as.integer(trajectories[i, ]) # Extracting the trajectory of i
    changes <- c(TRUE, diff(states) != 0) # Breakpoint detection
    episodes <- cumsum(changes) # Episode numbering
    durations <- table(episodes) # Duration of each episode
    episode_states <- states[changes] # Corresponding state
    
    # Construction of the resulting dataframe
    res[[i]] <- data.frame(
      id = i,
      state = episode_states,
      time = as.integer(durations)
    )
  }
  
  do.call(rbind, res)
}

# Find the best split among the different covariates
# We suppose that each covariates are binary categorical (2 modalities)
find_best_split <- function(df, covariates, min_leaf, pvalue_algo) {
  best <- list(pval = 1, var = NULL, left_level = NULL)
  
  for (var in names(covariates)) {
    x <- covariates[[var]]
    levs <- levels(x)   
    left_ids <- which(x==levs[1])
    
    df_left <- subset(df, id %in% left_ids)
    df_right <- subset(df, !(id %in% left_ids))
    
    if (length(unique(df_left$id)) < min_leaf || 
        length(unique(df_right$id)) < min_leaf) next 
    
    pval <- pvalue_algo(df_left, df_right)
    
    if (pval < best$pval)
      best <- list(pval = pval, var = var, left_level = levs[1], right_level = levs[2])
  }
  best
}


# Building the tree recursively 
build_tree <- function(df, covariates, pvalue_algo, min_obs = 20, min_leaf = 5,
                       alpha = 0.05, max_depth = 5, depth = 0) {
  pop_size <- length(unique(df$id))
  
  # Stopping criterion
  if (depth >= max_depth || pop_size < min_obs)
    return(list(type = "leaf", n = pop_size))
  
  best <- find_best_split(df, covariates, min_leaf, pvalue_algo)
  
  # No significant split
  if (is.null(best$var) || best$pval >= alpha)
    return(list(type = "leaf", n = pop_size))
  
  # Split and recursively construct the subtree
  left_ids <- which(covariates[[best$var]] == best$left_level)
  right_ids <- which(covariates[[best$var]] == best$right_level)
  df_left <- subset(df, id %in% left_ids)
  df_right <- subset(df, id %in% right_ids)
  
  list(
    type = "node",
    var = best$var,
    left_level = best$left_level,
    right_level = best$right_level,
    pval = best$pval,
    n = pop_size,
    left = build_tree(df_left, covariates, pvalue_algo, min_obs, min_leaf, alpha, max_depth, depth + 1),
    right = build_tree(df_right, covariates, pvalue_algo, min_obs, min_leaf, alpha, max_depth, depth + 1)
  )
}

# Print function 
print_tree <- function(node, indent = 0) {
  pad <- strrep("  ", indent)
  if (node$type == "leaf") {
    cat(pad, "└─ Leaf:  (n =", node$n, ")\n")
  } else {
    cat(pad, "├─ [", node$var, "] p =", formatC(node$pval, digits = 3, format = "e"), "\n")
    cat(pad, "  ├─ ==", node$left_level,  "\n"); print_tree(node$left,  indent + 2)
    cat(pad, "  └─ ==", node$right_level, "\n"); print_tree(node$right, indent + 2)
  }
}

# Use example w/ mvad data
source('src/two_samples_test.R')
library(TraMineR)
data(mvad)
trajectories <- mvad[, 17:86]
covariates <- mvad[, 3:14]
traj_df <- traj_to_df(trajectories)
D <- 6

# Function for p_value computation
alg <- function(df1, df2){
  likelihood_ratio_test(df1, df2, D, 'exponential')
}

# Tree construction
tree <- build_tree(traj_df, covariates, alg)

# Affichage
cat("=== ARBRE DE TRAJECTOIRES ===\n\n")
print_tree(tree)

