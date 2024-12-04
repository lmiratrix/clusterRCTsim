
#' Estimate ATEs for cluster RCT using data aggregation.
#'
#' This method aggregates data at the cluster level, and analyzes
#' using linear regression. Can optionally pass pre-aggregated data if
#' desired.
#'
#' This automatically does both weighting approaches.
#'
#' @inheritParams linear_model_estimators
#'
#' @export
aggregation_estimators <- function( formula,
                                    data = NULL,
                                    control_formula = NULL,
                                    control_interacted = FALSE,
                                    aggregated = FALSE ) {


    if ( is.null( control_formula ) ) {
        control_interacted = FALSE
    }

    # Utility to convert model results to easy to use tibble
    get_agg_ests <- function( M1, weight, name = "<Unknown>" ) {
        est1 <- M1$coefficients[["Z"]]
        se1 <- M1$std.error[["Z"]]
        pv1 <- M1$p.value[["Z"]]
        df1 <- M1$df[["Z"]]

        # Compile our results
        tibble(
            method = c( name ),
            weight = c( weight ),
            ATE_hat = c( est1 ),
            SE_hat = c( se1 ),
            p_value = c( pv1 ),
            df = c( df1 )
        )
    }

    if ( !is.null( formula ) && !aggregated ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula=control_formula )
    }

    datagg = data
    if ( !aggregated ) {
        datagg = aggregate_data( data = data, control_formula = control_formula)
        control_formula = attr( datagg, "control_formula" )
    }

    needFE = "blockID" %in% names( data )


    form = make_regression_formula( Yobs = "Ybar",
                                    FE = needFE,
                                    control_formula = control_formula,
                                    control_interacted = control_interacted )

    if ( control_interacted ) {
        datagg = center_controls( datagg, control_formula )
    }
    M3 <- lm_robust( form, data=datagg, se_type = "HC2" )
    Agg_FE_cluster = get_agg_ests( M3, "Cluster", ifelse( needFE, "LRa-FE-het", "LRa-het" ) )


    if ( control_interacted ) {
        datagg = center_controls( datagg, control_formula, weights = datagg$n )
    }
    M4 <- lm_robust_quiet( form, data=datagg, weights = n, se_type = "HC2" )
    Agg_FE_person = get_agg_ests( M4, "Person", ifelse( needFE, "LRapw-FE-het", "LRapw-het" ) )


    if ( !needFE ) {
        return( bind_rows( Agg_FE_cluster, Agg_FE_person ) )
    }



    # Interacted estimator that estimate within block and then average
    formI = make_regression_formula( Yobs = "Ybar",
                                     FE = needFE, interacted = TRUE,
                                     control_formula = control_formula,
                                     control_interacted = control_interacted )

    M5 <- lm_robust_quiet( formI, data=datagg, se_type = "HC2" )
    df = nrow( datagg ) - length( coef( M5 ) )
    aggd <- datagg %>% group_by( blockID ) %>%
        summarise( n = sum( n ),
                   J = n() )
    Agg_FI = generate_all_interacted_estimates( M5, aggd,
                                                aggregated = TRUE,
                                                use_full_vcov=TRUE,
                                                method = "LRa",
                                                weight = "cw",
                                                se_method = "het")
    Agg_FI$df = df

    M6 <- lm_robust_quiet( formI, data=datagg, se_type = "HC2", weights = n )

    Agg_wFI = generate_all_interacted_estimates( M6, aggd,
                                                 aggregated = TRUE,
                                                 use_full_vcov=TRUE,
                                                 method = "LRapw",
                                                 weight = "pw",
                                                 se_method = "het")
    Agg_wFI$df = df #M6$df[[1]]
    # Note: The df for heteroskedastic seems to be an overall single
    # number, so the satterwhite approximation of the above is
    # inappropriate.
    # Note2: Setting to just obs and covariates, ignoring weights

    # Drop SEs if there are singleton treated or control blocks.
    if ( has_singleton_clusters_agg( datagg ) || df <= 0 ) {
        Agg_FI$SE_hat = NA
        Agg_FI$p_value = NA
        Agg_FI$df = df
        Agg_wFI$SE_hat = NA
        Agg_wFI$p_value = NA
        Agg_wFI$df = df
    }

    # Compile our results
    res <- bind_rows( Agg_FE_cluster,
                      Agg_FE_person,
                      Agg_FI,
                      Agg_wFI )


    res
}




#### Testing/Demo code ####

if ( FALSE ) {

    data( fakeCRT )

    fakeCRT

    formula =  Yobs ~ T.x | S.id | D.id
    control_formula = ~ X.jk + C.ijk
    data = clusterRCT:::make_canonical_data(formula=formula, data=fakeCRT,
                                            control_formula = control_formula)

    aggregated = FALSE

    aggregation_estimators(formula, control_formula = control_formula,
                           data = fakeCRT )

    aggregation_estimators( Yobs ~ T.x | S.id | D.id,
                            data = fakeCRT )

    aggregation_estimators( Yobs ~ T.x | S.id,
                            data = fakeCRT )

}



