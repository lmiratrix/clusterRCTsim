
# Utility functions for making the plots in the paper.
#
# by Colin


# Utility for moving legend inside plot ----

legend_inside <- function(x = 0.98, y = 0.98) {
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(x, y),
    legend.justification = c(1, 1),
    legend.background = ggplot2::element_rect(
      fill = "white", colour = "grey40", linewidth = 0.4
    ),
    legend.margin = ggplot2::margin(5, 7, 5, 7)
  )
}


# My custom dot plot maker so I can have "x" marks ----

# A version of ggplot2's dot plot that places an X inside selected dots.
# Source this file once before calling geom_dotplot_marked().

marked_dot_grob <- function(dots, marked, cross_fraction = 0.55) {
  grid::gTree(
    dots = dots, marked = marked, cross_fraction = cross_fraction,
    cl = "marked_dot_grob"
  )
}

makeContent.marked_dot_grob <- function(x) {
  # Ask ggplot2's original dot-stack grob for the finished circles. Its
  # coordinates include binning, stacking, dotsize, stackratio, and facets.
  circles <- grid::makeContext(x$dots)
  if (!any(x$marked)) {
    return(grid::setChildren(x, grid::gList(circles)))
  }
  
  cx <- grid::convertX(circles$x, "mm", valueOnly = TRUE)[x$marked]
  cy <- grid::convertY(circles$y, "mm", valueOnly = TRUE)[x$marked]
  half_cross <- grid::convertUnit(circles$r, "mm", valueOnly = TRUE) *
    x$cross_fraction
  cross_style <- grid::gpar(col = "black", lwd = 0.6)
  
  grid::setChildren(x, grid::gList(
    circles,
    grid::segmentsGrob(
      cx - half_cross, cy - half_cross,
      cx + half_cross, cy + half_cross,
      default.units = "mm", gp = cross_style
    ),
    grid::segmentsGrob(
      cx - half_cross, cy + half_cross,
      cx + half_cross, cy - half_cross,
      default.units = "mm", gp = cross_style
    )
  ))
}

# grid needs the method registered to use it while rendering a PDF.
registerS3method(
  "makeContent", "marked_dot_grob", makeContent.marked_dot_grob,
  envir = asNamespace("grid")
)

GeomDotplotMarked <- ggplot2::ggproto(
  "GeomDotplotMarked", ggplot2::GeomDotplot,
  draw_group = function(data, panel_params, coord, ...) {
    dots <- ggplot2::GeomDotplot$draw_group(data, panel_params, coord, ...)
    marked <- as.character(data$mark) == "yes"
    marked_dot_grob(dots, !is.na(marked) & marked)
  }
)

geom_dotplot_marked <- function(...) {
  layer <- ggplot2::geom_dotplot(...)
  layer$geom <- GeomDotplotMarked
  layer
}


# The wrapper for making the plots ----



#' Make plot of differences and save to file
#' 
makeDotplot <- function(df, var, 
                        x_mark_var = NULL,
                        xlab       = NULL, 
                        ylab       = NULL, 
                        title      = NULL, 
                        subtitle   = NULL, 
                        binw       = NULL,
                        filename   = NULL,
                        coordratio = NULL,
                        axisfs     = 12,
                        xscale     = NULL,
                        plotmargin = NULL,
                        tcbreaks   = NULL,
                        tcratio    = F,
                        dim        = NULL,
                        rmlegend   = F,
                        byest      = F) {
  
  if ( !is.null( x_mark_var ) ) {
    df$is_under = df[[x_mark_var]]
  } else {
    df$is_under = "no"
  }
  
 
  
  p <- ggplot(df, aes(x = !!sym(var), mark = is_under)) +
    geom_dotplot_marked( aes( fill = low_df ),
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
         fill = "DF" ) +
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
    my_ggsave( p, filename, dim )
  }
  
  return(p)
}



my_ggsave <- function( p, filename, dim ) {
  out_dir <- dirname( paste0(PLOT_PATH, filename) )
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  if (is.null(dim)) {
    dim = c( NA, NA )
  }
  ggsave( paste0(PLOT_PATH, filename, '.pdf'), 
          plot = p,
          device = cairo_pdf,
          width = dim[1],
          height = dim[2])
  
}



