


#### Linear model estimators ####


# This helps catch warnings thrown when weights get weird with degrees
# of freedom issues.
lm_robust_quiet <- function( ... ) {
    M6 <- NULL

    # TODO: For small df datasets lm_robust can fall on its face.
    # This catches warnings when that occurs.
    withCallingHandlers(
        {
            M6 <- estimatr::lm_robust( ... )
        }, warning = function(w) {
            message("lm_robust threw warning due to df issues: ", conditionMessage(w))
            invokeRestart("muffleWarning")  # Suppresses the warning
        }
    )
    M6
}


#' Estimate ATEs for a cluster RCT using linear model
#'
#' @param formula Formula for outcome and treatment and nesting.  If
#'   NULL, data is assumed to be in canonical form (see vignette for
#'   further discussion).
#'
#' @importFrom estimatr lm_robust
#'
#' @export
linear_model_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL,
                                     weight = c( "Person", "Cluster" ) ) {

    weight = match.arg(weight)
    est_method = ifelse( weight == "Person", "LRi", "LRicw" )

    require( estimatr )

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }


    # Make our weights variable
    if (weight == "Person") {
        data$.weight <- 1
    } else {
        data <- data %>%
            group_by( clusterID ) %>%
            mutate( .weight = 1 / n() ) %>%
            ungroup()
    }


    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    form = make_regression_formula( FE = needFE,
                                    control_formula = control_formula )

    M2 <- lm_robust_quiet( form, data=data, clusters=clusterID, weights = .weight )
    est2 <- M2$coefficients[["Z"]]
    se2  <- M2$std.error[["Z"]]
    pv2 <- M2$p.value[["Z"]]
    df2 <- M2$df[[ "Z" ]]

    # Compile our results
    nm = paste0( est_method,
                 ifelse( needFE, "-FE-crve", "-crve" ) )

    tibble(
        method = c( nm ),
        ATE_hat = c( est2 ),
        weight = weight,
        SE_hat = c( se2 ),
        p_value = c( pv2 ),
        df = c( df2 )
    )

}



#' Estimate ATEs for a cluster RCT using linear model
#'
#' Estimate district-average effects using an interacted FE model and
#' then aggregate.
#'
#' @inheritParams linear_model_estimators
#' @param use_full_vcov TRUE/FALSE. When calculating standard errors,
#'   should the full variance-covariance matrix of the block-level
#'   estimates be used, or just the diagonal of standard errors?
#' @importFrom estimatr lm_robust
#'
#' @export
interacted_linear_model_estimators <- function( formula,
                                                data = NULL,
                                                control_formula = NULL,
                                                weight = c( "Person", "Cluster" ),
                                                use_full_vcov = FALSE ) {

    require( estimatr )
    weight = match.arg(weight)
    est_method = ifelse( weight == "Person", "LRi", "LRicw" )

    if ( !is.null( formula ) ) {
        data = clusterRCT:::make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }

    # Make our weights variable
    if (weight == "Person") {
        data$.weight <- 1
    } else {
        data <- data %>%
            group_by( clusterID ) %>%
            mutate( .weight = 1 / n() ) %>%
            ungroup()
    }


    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    if ( !needFE ) {
        return(     tibble(
            method = c(),
            ATE_hat = c(),
            SE_hat = c(),
            p_value = c()
        ) )
    }


    form = make_regression_formula( FE = TRUE, interacted = TRUE,
                                    control_formula = control_formula )

    M0.int <- lm_robust_quiet( form, data=data, clusters=clusterID, weights = .weight )

    ests <- generate_all_interacted_estimates( M0.int, data,
                                               method = est_method,
                                               weight = weight,
                                               se_method = "crve",
                                               use_full_vcov = use_full_vcov )

    ests
}




# Testing code ----


if ( FALSE ) {

    data( fakeCRT )
    fakeCRT

    linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                             control_formula = ~ X.jk + C.ijk,
                             data = fakeCRT )

    linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                             control_formula = ~ X.jk + C.ijk,
                             weight = "Cluster",
                             data = fakeCRT )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                        control_formula = ~ X.jk + C.ijk,
                                        data = fakeCRT )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                        control_formula = ~ X.jk + C.ijk,
                                        data = fakeCRT,
                                        weight = "Cluster" )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                        data = fakeCRT )

    # No interacted methods when no blocks
    interacted_linear_model_estimators( Yobs ~ T.x | S.id,
                                        data = fakeCRT )

    compare_methods( Yobs ~ T.x | S.id | D.id,
                     control_formula = ~ X.jk + C.ijk,
                     data = fakeCRT, include_method_characteristics = FALSE ) %>%
        knitr::kable( digits = 2 )
}



