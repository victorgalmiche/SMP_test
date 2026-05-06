# =============================================================================
# ARBRE DE RÉGRESSION POUR TRAJECTOIRES DE VIE
# Critère de split : p-value d'un two-sample test (fonction personnalisable)
# =============================================================================
# STRUCTURE DES DONNÉES ATTENDUE :
#   - trajectories : liste ou data.frame, une ligne = une trajectoire
#   - covariates   : data.frame, une ligne = un individu, colonnes = covariables
#   - Les deux objets doivent avoir le même nombre de lignes (même ordre)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. CALCUL p-value
# -----------------------------------------------------------------------------

source('src/two_samples_test.R')

traj_to_semimarkov <- function(trajectories, state_map = NULL) {
  
  if (is.null(state_map)) {
    all_states <- sort(unique(as.character(unlist(trajectories))))
    state_map  <- setNames(seq_along(all_states), all_states)
  }
  
  results <- lapply(seq_len(nrow(trajectories)), function(i) {
    states <- state_map[as.character(unlist(trajectories[i, ]))]
    
    rle_out <- rle(states)
    lengths <- rle_out$lengths
    values  <- rle_out$values
    
    n_transitions <- length(values) - 1
    if (n_transitions == 0) return(NULL)
    
    data.frame(
      id      = i,
      state.h = rle_out$values[-length(rle_out$values)],
      state.j = rle_out$values[-1],
      time    = rle_out$lengths[-length(rle_out$lengths)],
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, results)
  attr(df, "state_map") <- state_map
  df
}

# -----------------------------------------------------------------------------
# 2. FONCTIONS UTILITAIRES
# -----------------------------------------------------------------------------

# Trouve le meilleur split pour une covariable numérique
best_split_numeric <- function(covar_values, trajectories, test_fn, min_node_size) {
  n <- length(covar_values)
  sorted_vals <- sort(unique(covar_values))
  
  # Les seuils candidats sont les moyennes entre valeurs consécutives
  thresholds <- (sorted_vals[-length(sorted_vals)] + sorted_vals[-1]) / 2
  
  best_pval  <- 1
  best_thresh <- NA
  
  for (thresh in thresholds) {
    idx_left  <- which(covar_values <= thresh)
    idx_right <- which(covar_values >  thresh)
    
    if (length(idx_left) < min_node_size || length(idx_right) < min_node_size) next
    
    pval <- tryCatch(
      {p <- test_fn(trajectories[idx_left, , drop = FALSE], 
                   trajectories[idx_right, , drop = FALSE])
      if (is.na(p) || is.nan(p)) 1 else p},
      error = function(e) 1
    )
    
    if (pval < best_pval) {
      best_pval   <- pval
      best_thresh <- thresh
    }
  }
  list(pval = best_pval, threshold = best_thresh)
}

# Trouve le meilleur split pour une covariable catégorielle
best_split_categorical <- function(covar_values, trajectories, test_fn, min_node_size) {
  levels_covar <- unique(covar_values)
  if (length(levels_covar) < 2) return(list(pval = 1, groups = NULL))
  
  # Génère toutes les bipartitions non triviales des niveaux
  generate_bipartitions <- function(lvls) {
    n <- length(lvls)
    parts <- list()
    for (k in 1:(2^(n-1) - 1)) {
      bits <- as.integer(intToBits(k))[1:n]
      left  <- lvls[bits == 1]
      right <- lvls[bits == 0]
      parts[[length(parts) + 1]] <- list(left = left, right = right)
    }
    parts
  }
  
  partitions <- generate_bipartitions(levels_covar)
  best_pval  <- 1
  best_groups <- NULL
  
  for (part in partitions) {
    idx_left  <- which(covar_values %in% part$left)
    idx_right <- which(covar_values %in% part$right)
    
    if (length(idx_left) < min_node_size || length(idx_right) < min_node_size) next
    
    pval <- tryCatch(
      {p <- test_fn(trajectories[idx_left, , drop = FALSE], 
                    trajectories[idx_right, , drop = FALSE])
      if (is.na(p) || is.nan(p)) 1 else p},
      error = function(e) 1
    )
    
    if (pval < best_pval) {
      best_pval   <- pval
      best_groups <- part
    }
  }
  list(pval = best_pval, groups = best_groups)
}


# -----------------------------------------------------------------------------
# 3. CONSTRUCTION RÉCURSIVE DE L'ARBRE
# -----------------------------------------------------------------------------

build_tree <- function(
    trajectories,   # matrice/data.frame des trajectoires (n x T)
    covariates,     # data.frame des covariables (n x p)
    indices,        # indices courants dans le nœud
    test_fn,        # fonction two_sample_test
    alpha        = 0.05,   # seuil de p-value pour accepter un split
    max_depth    = 5,      # profondeur maximale
    min_node_size = 10,    # taille minimale d'un nœud
    depth        = 0       # profondeur courante (usage interne)
) {
  
  node <- list(
    indices  = indices,
    n        = length(indices),
    depth    = depth,
    is_leaf  = TRUE,
    split    = NULL,
    left     = NULL,
    right    = NULL
  )
  
  # Conditions d'arrêt
  if (depth >= max_depth || length(indices) < 2 * min_node_size) {
    return(node)
  }
  
  traj_node  <- trajectories[indices, , drop = FALSE]
  covar_node <- covariates[indices, , drop = FALSE]
  
  best_pval  <- 1
  best_var   <- NULL
  best_split_info <- NULL
  
  # Cherche le meilleur split parmi toutes les covariables
  for (var in names(covar_node)) {
    vals <- covar_node[[var]]
    
    if (is.numeric(vals) || is.integer(vals)) {
      res <- best_split_numeric(vals, traj_node, test_fn, min_node_size)
      split_type <- "numeric"
    } else {
      res <- best_split_categorical(vals, traj_node, test_fn, min_node_size)
      split_type <- "categorical"
    }
    
    if (res$pval < best_pval) {
      best_pval <- res$pval
      best_var  <- var
      best_split_info <- c(res, list(type = split_type))
    }
  }
  
  # Applique le split si la p-value est significative
  if (!is.null(best_var) && best_pval <= alpha) {
    vals <- covar_node[[best_var]]
    
    if (best_split_info$type == "numeric") {
      thresh   <- best_split_info$threshold
      idx_left  <- indices[vals <= thresh]
      idx_right <- indices[vals >  thresh]
      split_desc <- list(
        variable  = best_var,
        type      = "numeric",
        threshold = thresh,
        pvalue    = best_pval
      )
    } else {
      left_grp  <- best_split_info$groups$left
      right_grp <- best_split_info$groups$right
      idx_left  <- indices[vals %in% left_grp]
      idx_right <- indices[vals %in% right_grp]
      split_desc <- list(
        variable    = best_var,
        type        = "categorical",
        left_levels = left_grp,
        pvalue      = best_pval
      )
    }
    
    if (length(idx_left) >= min_node_size && length(idx_right) >= min_node_size) {
      node$is_leaf <- FALSE
      node$split   <- split_desc
      node$left    <- build_tree(trajectories, covariates, idx_left,
                                 test_fn, alpha, max_depth, min_node_size, depth + 1)
      node$right   <- build_tree(trajectories, covariates, idx_right,
                                 test_fn, alpha, max_depth, min_node_size, depth + 1)
    }
  }
  
  node
}


# -----------------------------------------------------------------------------
# 4. AFFICHAGE DE L'ARBRE
# -----------------------------------------------------------------------------

print_tree <- function(node, prefix = "", is_last = TRUE) {
  connector  <- if (is_last) "\u2514\u2500\u2500 " else "\u251c\u2500\u2500 "
  child_pref <- if (is_last) "    " else "\u2502   "
  
  if (node$is_leaf) {
    cat(prefix, connector,
        "[FEUILLE] n = ", node$n,
        " | profondeur = ", node$depth, "\n", sep = "")
  } else {
    sp <- node$split
    if (sp$type == "numeric") {
      split_str <- paste0(sp$variable, " <= ", round(sp$threshold, 3))
    } else {
      split_str <- paste0(sp$variable, " in {", paste(sp$left_levels, collapse = ", "), "}")
    }
    cat(prefix, connector,
        "[SPLIT] ", split_str,
        "  (p = ", format(sp$pvalue, digits = 3, scientific = TRUE), ")",
        "  n = ", node$n, "\n", sep = "")
    
    print_tree(node$left,  paste0(prefix, child_pref), is_last = FALSE)
    print_tree(node$right, paste0(prefix, child_pref), is_last = TRUE)
  }
}


# -----------------------------------------------------------------------------
# 5. PRÉDICTION (assignation d'un individu à une feuille)
# -----------------------------------------------------------------------------

predict_node <- function(node, new_covariates) {
  # new_covariates : data.frame d'une seule ligne
  if (node$is_leaf) return(node)
  
  sp  <- node$split
  val <- new_covariates[[sp$variable]]
  
  go_left <- if (sp$type == "numeric") {
    val <= sp$threshold
  } else {
    val %in% sp$left_levels
  }
  
  if (go_left) predict_node(node$left,  new_covariates)
  else         predict_node(node$right, new_covariates)
}


# -----------------------------------------------------------------------------
# 6. EXTRACTION DES FEUILLES (pour analyse post-hoc)
# -----------------------------------------------------------------------------

get_leaves <- function(node, leaves = list()) {
  if (node$is_leaf) {
    leaves[[length(leaves) + 1]] <- node
  } else {
    leaves <- get_leaves(node$left,  leaves)
    leaves <- get_leaves(node$right, leaves)
  }
  leaves
}


# =============================================================================
# EXEMPLE D'UTILISATION
# =============================================================================

data(mvad)
trajectories <- mvad[, 17:86]
covariates <- mvad[, 1:14]
traj_df <- traj_to_semimarkov(trajectories)
D <- 6
# Calculer le state_map une fois sur toutes les données
all_states <- sort(unique(as.character(unlist(trajectories))))
state_map  <- setNames(seq_along(all_states), all_states)

# Le passer au wrapper
make_test_fn <- function(test_fn, state_map) {
  function(traj_group1, traj_group2) {
    sm1 <- traj_to_semimarkov(traj_group1, state_map)
    sm2 <- traj_to_semimarkov(traj_group2, state_map)
    test_fn(sm1, sm2, D)
  }
}

# Construction de l'arbre
tree <- build_tree(
  trajectories  = trajectories,
  covariates    = covariates,
  indices       = 1:712,
  test_fn       = make_test_fn(likelihood_ratio_test, state_map),
  alpha         = 0.05,
  max_depth     = 4,
  min_node_size = 15
)

# Affichage
cat("=== ARBRE DE TRAJECTOIRES ===\n\n")
print_tree(tree)

# Feuilles
leaves <- get_leaves(tree)
cat("\n\nNombre de feuilles :", length(leaves), "\n")
cat("Tailles des feuilles :", sapply(leaves, `[[`, "n"), "\n")

# Prédiction pour un nouvel individu
new_ind <- data.frame(age = 35, gender = "F", diplome = "bac+5", statut = "emploi")
leaf <- predict_node(tree, new_ind)
cat("\nNouvel individu assigné à la feuille de profondeur", leaf$depth,
    "avec", leaf$n, "individus\n")
