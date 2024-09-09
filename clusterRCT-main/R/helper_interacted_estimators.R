

# Calculate the Welch-Satterthwaite degrees of freedom
#
# Calculate for a weighted sum of variances, with individual degrees
# of freedom for each (estimated) variance.
#
# @param SE2 A vector of variances.
# @param weight A vector of weights.
#
# @return The Welch-Satterthwaite degrees of freedom.
Welch_Satterthwaite_df <- function( weight, SE2, df ) {

    return( sum( weight * SE2 ) ^ 2 / sum( ( weight * SE2 )^2 / df ) )
    # Formula from Wiki: https://en.wikipedia.org/wiki/Welch%E2%80%93Satterthwaite_equation

}


# Calculate ATE and SE given weights and individual block estimates
#
# The function takes in the estimated ATEs and their standard errors,
# and the sizes of the clusters. The function then calculates the ATE
# and SE using different weighting schemes.

# @return A tibble with the ATE, SE, and p-value for each weighting
#   scheme.
calc_agg_estimate <- function( wts, ATE_hats, SE2_hats, df ) {

    # normalize weights
    wts = wts / sum(wts)

    ATE_hat <- weighted.mean(ATE_hats, wts)
    SE_hat = NA
    if ( is.matrix(SE2_hats) ) {
        wts = as.matrix(wts, ncol=1)
        SE_hat = as.numeric( sqrt( t(wts) %*% SE2_hats %*% wts ) )
    } else {
        SE_hat = sqrt( sum(wts ^ 2 * SE2_hats) )
    }
    tibble( ATE_hat = ATE_hat,
            SE_hat = SE_hat,
            p_value = 2*pnorm( abs(ATE_hat) / SE_hat, lower.tail=FALSE ),
            df = Welch_Satterthwaite_df( weight=wts, SE2=SE_hat^2, df=df ) )
}


# @param sizes A data frame with columns blockID, J, n, and
#   (optionally) weight
calculate_avg_impacts <- function( ATE_hats,
                                   SE2_hats,
                                   df,
                                   sizes = NULL,
                                   clusterID = NULL,
                                   blockID = NULL ) {

    if ( is.null( sizes ) ) {
        stopifnot( !is.null( clusterID ) && !is.null( blockID ) )
    } else {
        stopifnot( is.null( clusterID ) && is.null( blockID ) )
    }

    if ( is.null( sizes ) ) {
        sizes <- count_block_sizes( clusterID, blockID )
    }

    stopifnot( all( sizes$blockID == names(ATE_hats) ) )

    K = nrow( sizes )

    # equal weighting
    ATE_eq <- calc_agg_estimate( rep(1,K), ATE_hats, SE2_hats, df=df )

    # weight by number of clusters
    ATE_cluster <- calc_agg_estimate(sizes$J, ATE_hats, SE2_hats, df=df )

    # weight by number of students
    ATE_indiv <- calc_agg_estimate(sizes$n, ATE_hats, SE2_hats, df=df )

    rsp <- bind_rows( Block = ATE_eq,
                      Cluster = ATE_cluster,
                      Person = ATE_indiv, .id="weight" )

    # if ( "weight" %in% names(sizes) ) {
    #     # weight by arbitrary weight, if given
    #     ATE_weight = calc_agg_estimate(sizes$weight, ATE_hats, SE2_hats)
    #     ATE_weight$weight = "Weight"
    #     rsp = bind_rows( rsp, ATE_weight )
    # }

    return( rsp )
}


#' Average block-level estimates to get overall ATE
#'
#' Also calculate associated standard error.
#'
#' @inheritParams interacted_linear_model_estimators
#' @param fitModel The model that has the interacted estimates in it.
#' @param method Prefix of the method.  Will add the weighting to
#'   stem.
#' @param aggregated TRUE means data is summarized at the Block level.
#'   False means it is not, and needs to be aggregated.
#' @export
generate_all_interacted_estimates <- function( fitModel, data,
                                               use_full_vcov = FALSE,
                                               SE_table = NULL,
                                               method = "LR_FI_CRVE",
                                               aggregated = FALSE,
                                               include_block_estimates = FALSE ) {
    J = length( unique( data$blockID ) )

    cc = NA
    if ( !( class(fitModel) %in% c( "lm", "lm_robust" ) ) ) {
        # Multilevel model!
        cc = fixef( fitModel )
    } else {
        cc = coef( fitModel )
    }

    ids <- grep( "Z:", names(cc) )
    stopifnot(length(ids) == J)

    VC <- vcov(fitModel)

    if ( class( fitModel ) %in% "lm_robust" ) {
        df = fitModel$df[ ids ]
    } else if ( class( fitModel ) == "lm" ) {
        # No good block-level degree of freedom due to homoskedasticity assumption
        df = rep( NA, length( ids ) )
    } else {
        # MLM
        ccc = summary( fitModel )$coefficients
        df = ccc[ids,"df"]
    }

    ATE_hats <- cc[ids]
    names(ATE_hats) = stringr::str_remove(names(ATE_hats), "Z:blockID")


    ests = NA
    SE2_hats = NA
    if ( !is.null( SE_table ) ) {
        stopifnot( all( names(ATE_hats) == SE_table$blockID ) )
        SE2_hats = SE_table$SE_hat^2
        ests <- calculate_avg_impacts( ATE_hats,
                                       SE2_hats,
                                       df = df,
                                       sizes = SE_table )
    } else {
        sizes = data
        if ( !aggregated ) {
            sizes <- count_block_sizes( data$clusterID, data$blockID )
        }

        if ( use_full_vcov ) {
            # This is the cautious way we don't need since we have 0s in
            # the off diagonal
            SE2_hats = diag( as.matrix( VC ) )[ids] # for table below
            ests <- calculate_avg_impacts( ATE_hats, VC[ids,ids],
                                           df = df,
                                           sizes = sizes )
        } else {
            SE2_hats <- diag( as.matrix( VC ) )[ids]
            ests <- calculate_avg_impacts( ATE_hats, SE2_hats,
                                           df = df,
                                           sizes = sizes )
        }
    }

    ests <- ests %>%
        mutate( method = paste0( method, "_", weight ) ) %>%
        relocate( method )

    # Add block-level estimates to the output as an attribute
    if ( include_block_estimates ) {
        tb = tibble( blockID = names(ATE_hats),
                     ATE_hat = ATE_hats,
                     SE_hat = sqrt( SE2_hats ) )
        attr( ests, "blocks" ) <- tb
    }

    # Now a check: if rank deficient, then give NAs for all ATEs
    # to signal something is wrong
    if ( any( is.na( cc ) ) ) {
        warning( "Some coefficients of linear model undefined, likely due to colinearity in covariates and block dummies", call. = FALSE )
        ests$ATE_hat = rep( NA, nrow(ests) )
        ests$SE_hat = rep( NA, nrow(ests) )
    }

    ests
}
