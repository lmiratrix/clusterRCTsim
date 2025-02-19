
# Simulation helper functions




#' This function takes a bunch of simulation results and groups
#' estimators based on whether they give the same point estimate (to a
#' high degree of precision).
#'
#' @param rps A data frame of simulation results. Columns are runID,
#'   method, and ATE_hat, at least.  Each row is a simulation result
#'   for a given method.  The runID groups the methods on the same
#'   dataset.
#' @param   Method will call two methods the same if they are within
#'   tolerance of each other.  This could cause chaining where members
#'   in a group are linked by another member but are not close.
#'   Beware this!
#' @return A list of groups of methods that give the same point
#'   estimate.
group_estimators <- function( rps, tolerance = 10^-7, return_matrix = FALSE ) {

    wid <- rps %>% dplyr::select( runID, method, ATE_hat ) %>%
        pivot_wider( names_from = method, values_from = ATE_hat ) %>%
        dplyr::select( -runID )


    #results <- combn(names(wid), 2, function(x) {
    #        mean_diff <- mean( abs( wid[[x[1]]] - wid[[x[2]]] ) )
    #       tibble(pair = paste(x, collapse = " - "), mean_difference = round( mean_diff, digits=10 ) )
    #  }, simplify = FALSE) %>% bind_rows()


    # Create a matrix where each element is the mean difference between two columns
    mean_diff_matrix <- outer(1:ncol(wid), 1:ncol(wid), Vectorize(function(x, y) {
        mean( abs( wid[[x]] - wid[[y]] ) )
    }))

    # Naming rows and columns for clarity
    colnames(mean_diff_matrix) <- names(wid)
    rownames(mean_diff_matrix) <- names(wid)

    mean_diff_matrix


    # Create a binary matrix where 1 indicates differences below the tolerance
    binary_matrix <- 0 + (mean_diff_matrix < tolerance)
    binary_matrix

    # Convert this binary matrix into a list of groups
    groups <- list()
    already_included <- c()

    ests = sort( rownames(binary_matrix) )
    for (i in ests) {
        if (!( i %in% already_included)) {
            # Find indices where the estimator is within tolerance of others
            in_group <- which(binary_matrix[i, ] == 1)
            group_names <- colnames(binary_matrix)[in_group]
            groups[[length(groups) + 1]] <- group_names
            already_included <- c(already_included, group_names)
        }
    }

    if ( return_matrix ) {
        return( mean_diff_matrix )
    } else {
        return( groups )
    }
}



#### Demo and testing code ####

if ( FALSE ) {

    set.seed( 1039 )

    library( tidyverse )
    library( clusterRCT )

    one_run <- function() {

        sim.data <- blkvar::generate_multilevel_data( n.bar = 50, J = 30, p = 0.3,
                                                      size.impact.correlate = 1,
                                                      variable.n = TRUE )
        sim.data <- mutate( sim.data,
                            Z = randomizr::cluster_ra(clusters = sid),
                            Yobs = ifelse( Z == 1, Y1, Y0 ) )

        c1 <- clusterRCT::compare_methods( Yobs ~ Z | sid,
                                           data=sim.data,
                                           include_method_characteristics = FALSE )

        c1
    }

    rps = map_df( 1:10, ~ one_run(), .id = "runID", .progress=TRUE )

    groups <- rps %>% group_estimators()
    groups


    # List the groups
    map_chr( groups, paste, collapse = ", " ) %>%
        as.data.frame() %>%
        mutate( grp = 1:n() ) %>%
        relocate( grp ) %>%
        knitr::kable()
}



