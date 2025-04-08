
# Utility functions for making the plots in the paper.
#
# by Colin


#' Make plot of differences and save to file
#' 
makeDotplot <- function(df, var, 
                        xlab       = NULL, 
                        ylab       = NULL, 
                        title      = NULL, 
                        subtitle   = NULL, 
                        binw       = NULL,
                        filename   = NULL,
                        coordratio = NULL,
                        axisfs     = NULL,
                        xscale     = NULL,
                        plotmargin = NULL,
                        tcbreaks   = NULL,
                        tcratio    = F,
                        dim        = NULL,
                        rmlegend   = F,
                        byest      = F) {
  
  p <- ggplot(df, aes(x = !!sym(var))) +
    geom_dotplot( aes( fill = low_df ),
                  binwidth   = binw,
                  dotsize    = 0.8,
                  stackratio = 1.2,
                  stackgroups = TRUE,
                  stroke = 1,
                  method     = 'histodot') +
    theme_minimal() +
    labs(x        = xlab,
         y        = ylab,
         title    = title,
         subtitle = subtitle,
   #      color = "Structure",
         fill = "df" ) +
    theme(axis.text.y        = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA)) +
    scale_color_manual(values = c("white", "black")) 
  
  if (byest == T) {
    p <- p +
      facet_grid(~ Estimand) +
      theme(strip.background = element_rect(fill = "lightgray"),
            panel.spacing = unit(2, "lines"))
  }
  
  if (rmlegend == T) {
    p <- p +
      theme(legend.position = 'none')
  }
  
  
  if (!is.null(xscale)) {
    p <- p +
      scale_x_continuous(breaks = xscale) +
      xlim(c(min(xscale), max(xscale)))
  }
  
  if (!is.null(plotmargin)) {
    p <- p +
      theme(plot.margin = margin(plotmargin))
  }
  
  if (!is.null(tcbreaks)) {
    if (tcratio) {
      p <- p +
        # coord_cartesian(xlim = c(min(tcbreaks), max(tcbreaks))) +
        scale_x_continuous(breaks = tcbreaks, 
                           labels = c(paste0(tcbreaks[1:(length(tcbreaks) - 1)], 'X'), paste0(tcbreaks[length(tcbreaks)], 'X+')))
      
    } else {
      p <- p +
        # coord_cartesian(xlim = c(min(tcbreaks), max(tcbreaks))) +
        scale_x_continuous(breaks = tcbreaks, 
                           labels = c(tcbreaks[1:(length(tcbreaks) - 1)], paste0(tcbreaks[length(tcbreaks)], '+')))
    }
  }
  
  if (!is.null(axisfs)) {
    p <- p +
      theme(axis.text  = element_text(size = axisfs * 0.8),
            axis.title = element_text(size = axisfs),
            strip.text = element_text(size = axisfs))
  }
  
  if (!is.null(coordratio)) {
    if (!is.null(tcbreaks)) {
      p <- p +
        coord_fixed(xlim = c(min(tcbreaks), max(tcbreaks)),
                    ratio = coordratio)
    } else {
      p <- p +
        coord_fixed(ratio = coordratio)
    }
  }
  
  if (!is.null(filename)) {
    if (!is.null(dim)) {
      ggsave( paste0(PLOT_PATH, filename, '.png'), plot = p,
              width = dim[1],
              height = dim[2])
    } else {
      ggsave( paste0(PLOT_PATH, filename, '.png'), plot = p)
    }
  }
  
  return(p)
}
