


#### Linear model estimators ####





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
                                     control_formula = NULL ) {

    require( estimatr )

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }

    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    form = make_regression_formula( FE = needFE,
                                    control_formula = control_formula )

    M2 <- estimatr::lm_robust( form, data=data, clusters=clusterID )
    est2 <- M2$coefficients[["Z"]]
    se2  <- M2$std.error[["Z"]]
    pv2 <- M2$p.value[["Z"]]
    df2 <- M2$df[[ "Z" ]]

    # Compile our results
    nm = ifelse( needFE, "LR_FE_CRVE", "LR_CRVE" )
    tibble(
        method = c( nm ),
        ATE_hat = c( est2 ),
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
                                     use_full_vcov = FALSE ) {

    require( estimatr )

    if ( !is.null( formula ) ) {
        data = clusterRCT:::make_canonical_data( formula=formula, data=data, control_formula=control_formula )
    }

    # Linear regression (possibly with fixed effects if block exists)
    needFE = "blockID" %in% names(data)
    if ( !needFE ) {
        return(     tibble(
            method = c(),
            ATE_hat = c(  ),
            SE_hat = c(  ),
            p_value = c(  )
        ) )
    }
    form = make_regression_formula( FE = TRUE, interacted = TRUE,
                                    control_formula = control_formula )

    M0.int <- estimatr::lm_robust( form, data=data, clusters=clusterID )

    ests <- generate_all_interacted_estimates( M0.int, data,
                                               method = "LR_FI_CRVE",
                                               use_full_vcov = use_full_vcov )

    ests
}





if ( FALSE ) {

    data( fakeCRT )
    fakeCRT

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                        control_formula = ~ X.jk + C.ijk,
                             data = fakeCRT )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                        data = fakeCRT )

    interacted_linear_model_estimators( Yobs ~ T.x | S.id,
                                        data = fakeCRT )

}



