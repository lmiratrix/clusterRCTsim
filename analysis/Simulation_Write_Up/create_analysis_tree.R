#' Make a regression tree with the simulation results 
#' 
create_analysis_tree <- function(data, outcome_var, 
                                 predictor_vars = c("is_FI", "is_MLM", "is_cw", "is_pw", "blk_type", 
                                                    "J.bar", "cor_J_tau", "sd.sch", "ICC.sch", 
                                                    "cv.sch", "num.cov", "p.bar"),
                                 num_leaves = 10, palette_name = "Greens", 
                                 plot_tree = TRUE, initial_cp = 0.001,
                                 tree_title = NULL) {
  
  
  name_map <- c(
    block_type = "blk_type",
    size.impact.cor.school = "cor_J_tau",
    num.cov.school = "num.cov",
    cv.school = "cv.sch",
    sd.school = "sd.sch",
    ICC.school = "ICC.sch"
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
  node_means <- tree_p$frame$yval
  norm_means <- (node_means - min(node_means)) / 
    (max(node_means) - min(node_means)) * 100
  norm_means <- round(norm_means) + 1
  
  # Create colors
  cols <- colorRampPalette(brewer.pal(9, palette_name))(131)
  node_colors <- cols[round(norm_means) + 1]
  
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
    
    prp(tree_p,
        box.col = node_colors,
        border.col = "black",
        node.fun = leaf_fun,
        varlen = -10,
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