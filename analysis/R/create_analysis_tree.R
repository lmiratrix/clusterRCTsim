

library(rpart)
library( rpart.plot )
require(RColorBrewer)

# Historical sibling of create_analysis_tree_minmax() in
# create_analysis_tree_minmax.R (same directory). Both are independent
# rpart-based tree-fitting functions with different pruning APIs (this one:
# num_leaves + a legacy column name_map + a signed color scale; the other:
# min_leaves/max_leaves). They are NOT interchangeable and NOT drop-in
# replacements for each other - this one is used by
# analysis/Simulation_WriteUp.Rmd and empirical/empirical_results_report.Rmd;
# the other is used by fingerprint/fingerprint_simulation_report.Rmd.
# Deliberately kept as two files/functions rather than merged so neither
# report's tree output changes silently.

#' Make a regression tree with the simulation results
#'
create_analysis_tree <- function(data, outcome_var,
                                 predictor_vars = c("is_FI", "is_MLM", "is_cw", "is_pw", "blk_type", 
                                                    "J.bar", "cor_J_tau", "sd.clu", "ICC.clu", 
                                                    "cv.clu", "num.cov", "p.bar"),
                                 num_leaves = 10, palette_name = "Greens", 
                                 plot_tree = TRUE, initial_cp = 0.001,
                                 tree_title = NULL,
                                 force_scale = NULL,
                                 signed = FALSE ) {
  
  # Check if Times New Roman is available, fallback to serif
  available_fonts <- names(pdfFonts())
  if ("Times" %in% available_fonts) {
    title_family <- "Times"
  } else if ("serif" %in% available_fonts) {
    title_family <- "serif"
  } else {
    title_family <- ""  # Use default
  }
  
  # Convert 11pt to cex (assuming 12pt default)
  title_cex <- 11/12
  
  name_map <- c(
    block_type = "blk_type",
    size.impact.cor.cluster = "cor_J_tau",
    num.cov.cluster = "num.cov",
    cv.cluster = "cv.clu",
    sd.cluster = "sd.clu",
    ICC.cluster = "ICC.clu"
  )
  predictor_vars <- if_else(predictor_vars %in% names(name_map),
                            name_map[predictor_vars],
                            predictor_vars)
  
  name_map <- name_map[names(name_map) %in% names(data)]
  name_map <- set_names(names(name_map), name_map)
  data <- rename(data, !!!name_map)
  
  
  
  # Create formula with the specified outcome variable and predictors
  predictors <- paste(predictor_vars, collapse = " + ")
  form <- as.formula(paste(outcome_var, "~", predictors))
  
  # Initial tree fit
  tree <- rpart(form, data = data, cp = initial_cp)
  
  # CP table analysis
  cp_table <- tree$cptable %>% as_tibble()
  
  # Find optimal CP using 1-SE rule
  min_xerror <- which.min(cp_table$xerror)
  xerror_se <- cp_table$xstd[[min_xerror]]
  xerror_target <- cp_table$xerror[[min_xerror]] + xerror_se
  cp_best <- min(cp_table$CP[cp_table$xerror <= xerror_target])
  
  # Alternative CP based on number of leaves
  row <- pmin(nrow(cp_table),
              1 + sum(!(cp_table[, "nsplit"] >= (num_leaves - 1))))
  cp_leaves <- cp_table[row, "CP"]
  
  # Prune tree using the more conservative CP
  final_cp <- max(cp_leaves, cp_best)
  tree_p <- prune(tree, cp = final_cp)
  
  # Color mapping for nodes
  
  
  # Create colors
  if ( signed ) {
    node_means <- tree_p$frame$yval
    if ( !is.null( force_scale ) ) {
      rng = force_scale
    } else {
    rng = max( abs( node_means ) )
    }
    norm_means <- (node_means / rng ) * 50
    cols <- colorRampPalette(c("green", "white", "red"))(131)
    node_colors <- cols[round(norm_means) + 66]
  } else {
    node_means <- tree_p$frame$yval
    if ( sign( min(node_means) ) != sign( max(node_means) ) ) {
      warning( "All node means have the same sign; consider using signed=TRUE for better color mapping." )
    }
    
    rng = (max(node_means) - min(node_means))
    if ( !is.null( force_scale ) ) {
      rng = force_scale
    }
    norm_means <- ((node_means - min(node_means)) / rng) * 100
    norm_means <- round(norm_means) + 1
    
    cols <- colorRampPalette(brewer.pal(9, palette_name))(131)
    node_colors <- cols[round(norm_means) + 1]
  }
  
  # Leaf function for clean display
  leaf_fun <- function(x, labs, digits, varlen) {
    ifelse(x$frame$var == "<leaf>",
           sprintf("%.2f", x$frame$yval),
           labs)
  }
  
  # Plot if requested
  if (plot_tree) {
    # Use custom title if provided, otherwise default
    plot_title <- if (is.null(tree_title)) paste("Decision Tree for", outcome_var) else tree_title
    
    # Set graphics parameters for title formatting
    old_par <- par(cex.main = title_cex, font.main = 1, family = title_family)
    on.exit(par(old_par))
    
    prp(tree_p,
        type = 5,
        box.col = node_colors,
        border.col = "black",
        node.fun = leaf_fun,
        uniform = TRUE,
        yspace = 0.2,
        varlen = -10,
        under = TRUE,
        faclen = 0,
        main = plot_title)
  }
  
  # Return results
  return(list(
    tree = tree_p,
    cp_table = cp_table,
    final_cp = final_cp,
    cp_best = cp_best,
    formula = form,
    node_colors = node_colors,
    data = data,
    predictor_vars = predictor_vars
  ))
}

# Example usage:
# Use default title
# result1 <- create_analysis_tree(sruns.graph.bias, "abs(bias.s.estimand)")

# Use custom title
# result2 <- create_analysis_tree(sruns.graph.bias, "bias.s.estimand", 
#                                tree_title = "Bias Analysis Tree")

# Use custom title with other parameters
# result3 <- create_analysis_tree(sruns.graph.bias, "some_outcome",
#                                predictor_vars = c("is_FI", "is_MLM", "J.bar"),
#                                tree_title = "Custom Analysis: My Tree Title",
#                                palette_name = "Blues")




library(tidyverse)
library(partykit)
library(RColorBrewer)

#' Make a regression tree with the simulation results 
#'
create_analysis_tree_v2 <- function(data, outcome_var, 
                                 predictor_vars = c("is_FI", "is_MLM", "is_cw", "is_pw", "blk_type", 
                                                    "J.bar", "cor_J_tau", "sd.clu", "ICC.clu", 
                                                    "cv.clu", "num.cov", "p.bar"),
                                 num_leaves = 10,
                                 palette_name = "Greens", 
                                 plot_tree = TRUE,
                                 mincriterion = 0.95,
                                 minsplit = 20,
                                 minbucket = 7,
                                 tree_title = NULL,
                                 force_scale = NULL,
                                 signed = FALSE,
                                 show_pval = FALSE,
                                 show_ids = FALSE) {
  
  available_fonts <- names(pdfFonts())
  if ("Times" %in% available_fonts) {
    title_family <- "Times"
  } else if ("serif" %in% available_fonts) {
    title_family <- "serif"
  } else {
    title_family <- ""
  }
  
  title_cex <- 11 / 12
  
  name_map <- c(
    block_type = "blk_type",
    size.impact.cor.cluster = "cor_J_tau",
    num.cov.cluster = "num.cov",
    cv.cluster = "cv.clu",
    sd.cluster = "sd.clu",
    ICC.cluster = "ICC.clu"
  )
  
  predictor_vars <- if_else(
    predictor_vars %in% names(name_map),
    name_map[predictor_vars],
    predictor_vars
  )
  
  name_map <- name_map[names(name_map) %in% names(data)]
  name_map <- set_names(names(name_map), name_map)
  data <- rename(data, !!!name_map)

  # Convert to plain data.frame and coerce logical columns to factor.
  # ctree with multiway = TRUE crashes on tibbles and logical predictors
  # (partykit bug: object 'X' not found in .ctree_test_1d).
  df <- as.data.frame(data)
  for (v in names(df)) {
    if (is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
  }

  predictors <- paste(predictor_vars, collapse = " + ")
  form <- as.formula(paste(outcome_var, "~", predictors))

  ctrl <- ctree_control(
    mincriterion = mincriterion,
    minsplit = minsplit,
    minbucket = minbucket,
    maxdepth = Inf,
    multiway = TRUE
  )

  tree <- ctree(form, data = df, control = ctrl)

  # Greedy weakest-link pruning: collapse the least-significant inner split
  # until n_leaves <= num_leaves.  Uses nodeprune() from partykit.
  while (length(nodeids(tree, terminal = TRUE)) > num_leaves) {
    inner_ids <- setdiff(nodeids(tree), nodeids(tree, terminal = TRUE))
    if (length(inner_ids) == 0) break

    # Each inner node stores the criterion (1 - p-value) used to make the split.
    # The weakest split has the smallest criterion value (highest p-value).
    criteria <- sapply(inner_ids, function(id) {
      info <- nodeapply(tree, ids = id, FUN = info_node)[[1]]
      crit <- info$criterion
      # criterion may be a named vector; take the max (the statistic that drove the split)
      if (is.null(crit)) return(Inf)
      max(crit, na.rm = TRUE)
    })

    weakest_id <- inner_ids[which.min(criteria)]
    tree <- nodeprune(tree, ids = weakest_id)
  }

  # Compute mean outcome per terminal node only (predict() returns terminal ids)
  terminal_ids <- nodeids(tree, terminal = TRUE)
  node_assignments <- predict(tree, type = "node")
  terminal_means <- sapply(terminal_ids, function(id) {
    mean(df[[outcome_var]][node_assignments == id], na.rm = TRUE)
  })

  if (signed) {
    rng <- if (!is.null(force_scale)) force_scale else max(abs(terminal_means), na.rm = TRUE)
    norm_means <- pmax(pmin((terminal_means / rng) * 50, 65), -65)
    cols <- colorRampPalette(c("green", "white", "red"))(131)
    terminal_colors <- cols[round(norm_means) + 66]
  } else {
    rng <- if (!is.null(force_scale)) force_scale else diff(range(terminal_means, na.rm = TRUE))
    if (rng == 0) {
      terminal_colors <- rep(brewer.pal(9, palette_name)[5], length(terminal_means))
    } else {
      norm_means <- ((terminal_means - min(terminal_means, na.rm = TRUE)) / rng) * 100
      norm_means <- pmax(pmin(round(norm_means) + 1, 131), 1)
      cols <- colorRampPalette(brewer.pal(9, palette_name))(131)
      terminal_colors <- cols[norm_means]
    }
  }

  n_leaves <- length(terminal_ids)

  if (plot_tree) {
    plot_title <- if (is.null(tree_title)) {
      paste("Decision Tree for", outcome_var)
    } else {
      tree_title
    }

    old_par <- par(cex.main = title_cex, font.main = 1, family = title_family)
    on.exit(par(old_par))

    # terminal_panel receives a partynode; look up color by node id
    t_ids   <- terminal_ids
    t_means <- terminal_means
    t_cols  <- terminal_colors
    n_asn   <- node_assignments
    colored_panel <- function(obj, ...) {
      nid     <- id_node(obj)
      idx     <- which(t_ids == nid)
      mean_val <- t_means[idx]
      n_obs    <- sum(n_asn == nid, na.rm = TRUE)
      grid::grid.rect(gp = grid::gpar(fill = t_cols[idx], col = "black"))
      grid::grid.text(sprintf("%.3f\nn=%d", mean_val, n_obs),
                      gp = grid::gpar(fontsize = 9))
    }

    plot(tree, main = plot_title, terminal_panel = colored_panel,
         ip_args = list(pval = show_pval, id = show_ids))
  }

  list(
    tree = tree,
    formula = form,
    node_colors = terminal_colors,
    terminal_ids = terminal_ids,
    terminal_means = terminal_means,
    data = df,
    predictor_vars = predictor_vars,
    n_leaves = n_leaves,
    control = ctrl
  )
}

