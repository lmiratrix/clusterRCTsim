
library( rpart.plot )
library( rpart )

# Historical sibling of create_analysis_tree() in create_analysis_tree.R (same
# directory). Both are independent rpart-based tree-fitting functions with
# different pruning APIs (this one: min_leaves/max_leaves; the other:
# num_leaves + a legacy column name_map + a signed color scale). They are NOT
# interchangeable and NOT drop-in replacements for each other - this one is
# used by fingerprint/fingerprint_simulation_report.Rmd; the other is used by
# analysis/Simulation_WriteUp.Rmd and empirical/empirical_results_report.Rmd.
# Deliberately kept as two files/functions rather than merged so neither
# report's tree output changes silently.

#' Make a regression tree with the simulation results
#'
create_analysis_tree_minmax <- function( data,
                                  outcome_var, 
                                  predictor_vars,
                                  min_leaves = 1,
                                  max_leaves = Inf, 
                                  palette_name = "Greens", 
                                  plot_tree = TRUE, 
                                  initial_cp = 0.001,
                                  tree_title = NULL,
                                  include_counts = FALSE
                                  ) {
  
  # Create formula with the specified outcome variable and predictors
  predictors <- paste(predictor_vars, collapse = " + ")
  form <- as.formula(paste(outcome_var, "~", predictors))
  
  # Initial tree fit
  tree <- rpart(form, data = data, cp = initial_cp)
  
  # CP table analysis
  cp_table <- tree$cptable %>% 
    as_tibble()
  
  # Find optimal CP using 1-SE rule
  min_xerror <- which.min(cp_table$xerror)
  xerror_se <- cp_table$xstd[[min_xerror]]
  xerror_target <- cp_table$xerror[[min_xerror]] + xerror_se
  cp_best <- max(cp_table$CP[cp_table$xerror <= xerror_target])
  
  # Alternative CP based on passed max number of leaves
  row <- pmin(nrow(cp_table),
              1 + sum(!(cp_table[, "nsplit"] >= (max_leaves - 1))))
  cp_leaves <- cp_table[row, "CP"]
  
  # Alternative CP based on passed min number of leaves
  row <- pmin( nrow(cp_table),
               sum( cp_table[, "nsplit"] <= min_leaves))
  cp_min_leaves <- cp_table[row, "CP"]
  
  
  # Prune tree using the more conservative CP
  final_cp <- min( cp_min_leaves, max(cp_leaves, cp_best) )
  tree_p <- prune(tree, cp = final_cp)
  
  # Color mapping for nodes
  node_means <- tree_p$frame$yval
  norm_means <- (node_means - min(node_means)) / 
    (max(node_means) - min(node_means)) * 100
  norm_means <- round(norm_means) + 1
  
  # Create colors
  cols <- colorRampPalette(RColorBrewer::brewer.pal(9, palette_name))(131)
  node_colors <- cols[round(norm_means) + 1]
  
  # Leaf function for clean display
  leaf_fun <- function(x, labs, digits, varlen) {
    fr <- x$frame
    is_leaf <- fr$var == "<leaf>"
    if ( include_counts ) {
      labs[is_leaf] <- sprintf("%.2f (%d)", 
                             fr$yval[is_leaf],
                             round(fr$n[is_leaf]) )
    } else {
      labs[is_leaf] <- sprintf("%.2f", fr$yval[is_leaf])
    }
    labs
  } 
  
  # Plot if requested
  if (plot_tree) {
    # Use custom title if provided, otherwise default
    plot_title <- if (is.null(tree_title)) paste("Decision Tree for", outcome_var) else tree_title
    
    prp(tree_p,
        box.col = node_colors,
        border.col = "black",
        node.fun = leaf_fun,
        varlen = -10,
        main = plot_title)
  }
  
  # Return results
  return( invisible( list(
    tree = tree_p,
    cp_table = cp_table,
    final_cp = final_cp,
    cp_best = cp_best,
    formula = form,
    node_colors = node_colors,
    data = data,
    predictor_vars = predictor_vars
  )) )
}
