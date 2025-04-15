

# Canoncial weight estimators plug in weights to verify the weight
# formula in the paper.
#
# They should exactly correspond to design based (except for the MLM
# weight ones).


#' Estimate ATE via canonical weights on the cluster means
#'
#' These estimators use the cluster means and the table of weights
#' shown in the main paper to estimate the ATE across all combinations
#' of weights at the block and the cluster level.
#'
#' There are no standard errors for this method.  It is purely to
#' validate how the other estimators can be represented as given with
#' the weight formula in the paper.
#'
#' @inheritParams linear_model_estimators
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
canonical_weight_estimators <- function( formula,
                                         data = NULL,
                                         data_agg = NULL ) {

    stopifnot( !is.null( data ) )
    data = make_canonical_data( formula=formula, data=data )

    if ( is.null( data_agg ) ) {
        # Collapse to clusters
        data_agg = aggregate_data( data )
    }

    has_block = "blockID" %in% names( data )
    if ( !has_block ) {
        data_agg$blockID = "Uni"
    }

    calc_est <- function( data_agg, cw, bw ) {
        if ( length( bw ) == 1 ) {
            bw = rep( bw, nrow(data_agg) )
        }
        if ( length( cw ) == 1 ) {
            cw = rep( cw, nrow(data_agg) )
        }
        if ( length( bw ) != nrow(data_agg) ) {
            stop( "Block weight vector must be same length as data" )
        }
        if ( length( cw ) != nrow(data_agg) ) {
            stop( "Cluster weight vector must be same length as data" )
        }

        data_agg$cw = cw
        data_agg$bw = bw
        s_dat <- data_agg %>%
            group_by( blockID, Z, bw ) %>%
            summarise( Ybar = weighted.mean( Ybar, w = cw ), .groups = "drop" )
        stopifnot( nrow( s_dat ) == length( unique( data_agg$blockID ) ) * 2 )
        s_dat <- s_dat %>%
            group_by( Z ) %>%
            summarise( ATE = weighted.mean( Ybar, w = bw ) )
        if ( nrow( s_dat ) != 2 ) {
            return( NA )
        }
        s_dat$ATE[[2]] - s_dat$ATE[[1]]
    }

    data_agg <- data_agg %>%
        group_by( blockID ) %>%
        mutate( N = sum( n ),
                p = weighted.mean( Z, w=n ),
                pc = mean( Z ),
                J = n() )

    person_person = calc_est( data_agg,
                              cw = data_agg$n,
                              bw = data_agg$N )

    person_FE = calc_est( data_agg,
                          cw = data_agg$n,
                          bw = data_agg$N * data_agg$p * (1-data_agg$p) )

    person_block = calc_est( data_agg,
                             cw = data_agg$n,
                             bw = 1 )

    cluster_cluster = calc_est( data_agg,
                                cw = 1,
                                bw = data_agg$J )
    cluster_FE = calc_est( data_agg,
                           cw = 1,
                           bw = data_agg$J * data_agg$pc * (1-data_agg$pc) )

    cluster_block = calc_est( data_agg,
                              cw = 1,
                              bw = 1 )

    form = make_regression_formula( FE = has_block, cluster_RE = TRUE )
    M0 = lmer( form, data=data )
    tau2 = VarCorr(M0)$clusterID[1,1]
    sigma2 = sigma(M0)^2
    formRE = make_regression_formula( FE = FALSE, cluster_RE = TRUE )
    formRE = update( formRE, . ~ . + (1+Z|blockID) )
    M1 = lmer( formRE, data=data )
    eta2 = VarCorr(M1)$blockID[ 2, 2 ]

    cw_MLM = 1 / (tau2 + sigma2 / data_agg$J)
    bw_RIRC = 1 / (eta2 + tau2 / (data_agg$J * data_agg$pc * (1-data_agg$pc) ))

    MLM_FI_cluster = calc_est( data_agg,
                               cw = cw_MLM,
                               bw = data_agg$J )
    MLM_FE = calc_est( data_agg,
                       cw = cw_MLM,
                       bw = data_agg$N * data_agg$p * (1-data_agg$p) )
    MLM_RIRC = calc_est( data_agg,
                         cw = cw_MLM,
                         bw = bw_RIRC )
    MLM_FI_block = calc_est( data_agg,
                             cw = cw_MLM,
                             bw = 1 )

    rs <- tibble( method = c( "Person-Person", "Person-FE", "Person-Block",
                              "MLM-FI-Cluster", "MLM-FE", "MLM-RIRC", "MLM-Block",
                              "Cluster-Cluster", "Cluster-FE", "Cluster-Block" ),
                  weight = method,
                  ATE_hat = c( person_person, person_FE, person_block,
                               MLM_FI_cluster, MLM_FE, MLM_RIRC, MLM_FI_block,
                               cluster_cluster, cluster_FE, cluster_block ) )
    rs$method = paste0( "WT_", rs$method )

    return( rs )
}


