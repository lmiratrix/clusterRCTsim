
#' Make a boxplot with custom whiskers and colors
#' 
create_custom_boxplot <- function(data, y_var, 
                                  title = NULL, y_label = NULL,
                                  vline_at = NULL, log_scale = FALSE) {
  # Default title and y-axis label if not provided
  if(is.null(title)) {
    title <- paste("Distribution of", y_var, "by Method")
  }
  
  if(is.null(y_label)) {
    y_label <- y_var
  }
  
  # Create the plot
  p <- ggplot(data, aes(x = method_ordered, y = .data[[y_var]], color = weight_biased))
  
  # Add vertical reference line first
  if(!is.null(vline_at)) {
    p <- p + geom_hline(yintercept = vline_at, linetype = "dashed", color = "black", alpha = 0.7)
  }
  
  # First, add the faint min/max lines in the background
  p <- p + geom_linerange(
    stat = "summary",
    fun.data = function(x) {
      data.frame(ymin = min(x), ymax = max(x))
    },
    alpha = 0.2,
    linewidth = 0.5
  )
  
  # Add the actual 1% and 99% whiskers
  p <- p + stat_summary(
    fun.data = function(x) {
      percentiles <- quantile(x, probs = c(0.01, 0.99))
      return(data.frame(ymin = percentiles[1], ymax = percentiles[2]))
    },
    geom = "errorbar",
    width = 0.5,
    position = position_dodge(width = 0.75)
  )
  
  # Add the custom boxplot for 10% to 90%
  p <- p + stat_summary(
    fun.data = function(x) {
      r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
      names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
      r
    },
    geom = "boxplot",
    alpha = 0.7,
    width = 0.5,
    position = position_dodge(width = 0.75)
  )
  
  # Complete the plot formatting
  p <- p + labs(
    title = title,
    x = "",
    y = y_label
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    ) +
    coord_flip()
  
  # Apply log scale if requested
  if(log_scale) {
    p <- p + scale_y_log10()
  }
  
  # Add the color scale
  p <- p + scale_color_manual(
    values = c(
      "Cluster_0" = "darkblue", 
      "Cluster_1" = "lightblue",
      "Person_0" = "darkred", 
      "Person_1" = "lightcoral"
    ),
    name = "Weight and Bias",
    labels = c(
      "Cluster_0" = "Cluster (Unbiased)",
      "Cluster_1" = "Cluster (Biased)",
      "Person_0" = "Person (Unbiased)",
      "Person_1" = "Person (Biased)"
    )
  )
  
  return(p)
}

#this is the custom boxplot but with colors based on different categories 
#' Make a boxplot with custom whiskers and colors
#' 
create_custom_boxplot2 <- function(data, y_var, 
                                  title = NULL, y_label = NULL,
                                  vline_at = NULL, log_scale = FALSE) {
  # Default title and y-axis label if not provided
  if(is.null(title)) {
    title <- paste("Distribution of", y_var, "by Method")
  }
  
  if(is.null(y_label)) {
    y_label <- y_var
  }
  
  # Create new variables for coloring
  data <- data %>%
    mutate(
      # Weight type: from the weight variable
      weight_type = weight,
      
      # Hierarchical method categorization
      method_type = case_when(
        is_MLM ~ "MLM",
        !is_MLM & is_FI ~ "FI", 
        TRUE ~ "FE"
      ),
      
      # Combined coloring variable
      color_group = paste(weight_type, method_type, sep = "_")
    )
  
  # Create the plot
  p <- ggplot(data, aes(x = method_ordered, y = .data[[y_var]], color = color_group))
  
  # Add vertical reference line first
  if(!is.null(vline_at)) {
    p <- p + geom_hline(yintercept = vline_at, linetype = "dashed", color = "black", alpha = 0.7)
  }
  
  # First, add the faint min/max lines in the background
  p <- p + geom_linerange(
    stat = "summary",
    fun.data = function(x) {
      data.frame(ymin = min(x), ymax = max(x))
    },
    alpha = 0.2,
    linewidth = 0.5
  )
  
  # Add the actual 2.5% and 97.5% whiskers
  p <- p + stat_summary(
    fun.data = function(x) {
      percentiles <- quantile(x, probs = c(0.025, 0.975))
      return(data.frame(ymin = percentiles[1], ymax = percentiles[2]))
    },
    geom = "errorbar",
    width = 0.5,
    position = position_dodge(width = 0.75)
  )
  
  # Add the custom boxplot for 10% to 90%
  p <- p + stat_summary(
    fun.data = function(x) {
      r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
      names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
      r
    },
    geom = "boxplot",
    alpha = 0.7,
    width = 0.5,
    position = position_dodge(width = 0.75)
  )
  
  # Complete the plot formatting
  p <- p + labs(
    title = title,
    x = "",
    y = y_label
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    ) +
    coord_flip()
  
  # Apply log scale if requested
  if(log_scale) {
    p <- p + scale_y_log10()
  }
  
  # Add the new color scale
  p <- p + scale_color_manual(
    values = c(
      "Cluster_MLM" = "darkblue", 
      "Cluster_FI" = "blue",
      "Cluster_FE" = "lightblue",
      "Person_MLM" = "darkred", 
      "Person_FI" = "red",
      "Person_FE" = "lightcoral"
    ),
    name = "Weight Type & Method",
    labels = c(
      "Cluster_MLM" = "Cluster - MLM",
      "Cluster_FI" = "Cluster - FI",
      "Cluster_FE" = "Cluster - FE", 
      "Person_MLM" = "Person - MLM",
      "Person_FI" = "Person - FI",
      "Person_FE" = "Person - FE"
    )
  )
  
  return(p)
}


#this is the custom boxplot but for the data overall and not by metho
create_custom_boxplot3 <- function(data, y_var, 
                                   title = NULL, y_label = NULL,
                                   vline_at = NULL, log_scale = FALSE) {
  # Default title and y-axis label if not provided
  if(is.null(title)) {
    title <- paste("Distribution of", y_var)
  }
  
  if(is.null(y_label)) {
    y_label <- y_var
  }
  
  # Create the plot - single boxplot with constant x value
  p <- ggplot(data, aes(x = factor(1), y = .data[[y_var]]))
  
  # Add vertical reference line first
  if(!is.null(vline_at)) {
    p <- p + geom_hline(yintercept = vline_at, linetype = "dashed", color = "black", alpha = 0.7)
  }
  
  # First, add the faint min/max lines in the background
  p <- p + geom_linerange(
    stat = "summary",
    fun.data = function(x) {
      data.frame(ymin = min(x), ymax = max(x))
    },
    alpha = 0.2,
    linewidth = 0.5,
    color = "black"
  )
  
  # Add the actual 2.5% and 97.5% whiskers
  p <- p + stat_summary(
    fun.data = function(x) {
      percentiles <- quantile(x, probs = c(0.025, 0.975))
      return(data.frame(ymin = percentiles[1], ymax = percentiles[2]))
    },
    geom = "errorbar",
    width = 0.5,
    color = "black"
  )
  
  # Add the custom boxplot for 10% to 90%
  p <- p + stat_summary(
    fun.data = function(x) {
      r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
      names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
      r
    },
    geom = "boxplot",
    alpha = 0.7,
    width = 0.5,
    fill = "gray80",
    color = "black"
  )
  
  # Complete the plot formatting
  p <- p + labs(
    title = title,
    x = "",
    y = y_label
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Hide x-axis labels since it's just one boxplot
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    ) +
    coord_flip()
  
  # Apply log scale if requested
  if(log_scale) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}
#this is a new function that creates a top n most of something table 

create_top_n_table <- function(data, 
                               column_to_rank_by, 
                               methods_to_include = NULL,
                               n = 5, 
                               abs_value = TRUE, 
                               descending = TRUE,
                               min_threshold = NULL,
                               columns_to_display = NULL) {
  
  # Make a copy of the data
  df <- data
  
  # Convert column_to_rank_by to symbol for dplyr operations
  col_sym <- sym(column_to_rank_by)
  
  # Filter methods if specified
  if (!is.null(methods_to_include)) {
    df <- df %>% filter(method %in% methods_to_include)
  }
  
  # Setup ranking
  if (abs_value) {
    # Use absolute value for ranking
    df <- df %>%
      group_by(method) %>%
      arrange(desc(abs(!!col_sym))) %>%
      slice_head(n = n)
  } else {
    # Use actual value for ranking
    if (descending) {
      df <- df %>%
        group_by(method) %>%
        arrange(desc(!!col_sym)) %>%
        slice_head(n = n)
    } else {
      df <- df %>%
        group_by(method) %>%
        arrange(!!col_sym) %>%
        slice_head(n = n)
    }
  }
  
  # Sort methods alphabetically
  df <- df %>% arrange(method)
  
  # Apply threshold if specified
  if (!is.null(min_threshold)) {
    if (abs_value) {
      df <- df %>% filter(abs(!!col_sym) > min_threshold)
    } else {
      if (descending) {
        df <- df %>% filter(!!col_sym > min_threshold)
      } else {
        df <- df %>% filter(!!col_sym < min_threshold)
      }
    }
  }
  
  # Select columns to display
  if (!is.null(columns_to_display)) {
    df <- df %>% select(all_of(columns_to_display))
  }
  
  # Convert to data frame
  df <- as.data.frame(df)
  
  # Wrap column names
  current_names <- names(df)
  for (i in 1:length(current_names)) {
    if (nchar(current_names[i]) > 7) {  # If name is longer than 7 characters
      # Find a good breaking point (space, underscore, etc.)
      mid_point <- ceiling(nchar(current_names[i]) / 2)
      
      # Look for a space near the middle
      space_pos <- as.numeric(regexpr("[ \\._]", substr(current_names[i], mid_point-2, mid_point+2)))
      
      if (space_pos > 0) {
        # Break at the space
        break_point <- mid_point - 3 + space_pos
        names(df)[i] <- paste0(
          substr(current_names[i], 1, break_point-1),
          "\n",
          substr(current_names[i], break_point+1, nchar(current_names[i]))
        )
      } else {
        # No space found, just break in the middle
        names(df)[i] <- paste0(
          substr(current_names[i], 1, mid_point),
          "\n",
          substr(current_names[i], mid_point+1, nchar(current_names[i]))
        )
      }
    }
  }
  
  
  # Return the data frame
  return(df)
}
