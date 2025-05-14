
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


