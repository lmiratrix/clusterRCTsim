



#' Cluster RCT with multi-level modeling
#'
#' This uses the lmerTest package to get the p-values, and fits a
#' model with assumed homoskedasticity, etc.  I.e., this is the
#' vanilla MLM that one would typically fit.
#'
#' @inheritParams linear_model_estimators
#'
#' @export
MLM_estimators <- function( formula,
                            data = NULL,
                            control_formula = NULL,
                            suppress_warnings = TRUE ) {

    require( lme4 )
    require( lmerTest )

    # If suppress warnings then wrap lmer so all messages and warnings
    # get buried.
    my_lmer <- lmerTest::lmer
    if ( suppress_warnings ) {
        my_lmer <- function( ... ) {
            rs_f <- purrr::quietly( lmerTest::lmer )
            rs_f( ... )$result
        }
    }

    get_mlm_ests <- function( M1, name = "MLM" ) {
        est1 <- lme4::fixef( M1 )[["Z"]]
        se1 <- arm::se.fixef( M1 )[["Z"]]
        ss = summary(M1)
        pv1 <- ss$coefficients["Z",5]
        df1 <- ss$coefficients["Z",3]
        # Compile our results
        tibble(
            method = c( name ),
            ATE_hat = c( est1 ),
            SE_hat = c( se1 ),
            p_value = c( pv1 ),
            df = c( df1 )
        )
    }

    if ( !is.null( formula ) ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula=control_formula )
    }

    needFE = "blockID" %in% names(data)

    if ( !needFE ) {
        # Do the only MLM that is possible, and bail.
        form = make_regression_formula( FE = FALSE,
                                        control_formula = control_formula,
                                        cluster_RE = TRUE )
        M1 <- my_lmer( form, data=data )

        return( get_mlm_ests( M1, "MLM") )
    }


    formFE = make_regression_formula( FE = TRUE,
                                      control_formula = control_formula,
                                      cluster_RE = TRUE )
    M1 <- my_lmer( formFE, data=data )
    MLM_FE = get_mlm_ests(M1, "MLM_FE")

    formRE = make_regression_formula( FE = FALSE,
                                      control_formula = control_formula,
                                      cluster_RE = TRUE )
    formRE = update( formRE, . ~ . + (1|blockID) )
    M1 <- my_lmer( formRE, data=data )
    MLM_RE = get_mlm_ests(M1, "MLM_RE")


    formFI = make_regression_formula( FE = TRUE, interacted = TRUE,
                                      control_formula = control_formula,
                                      cluster_RE = TRUE )
    M1 <- my_lmer( formFI, data=data )
    MLM_FI = clusterRCT:::generate_all_interacted_estimates( M1, data,
                                                             use_full_vcov = TRUE,
                                                             method = "MLM_FI" )
    # Hack conservative DF calculation:
    # df = #clusters - 2 * #blocks - #covariates
    MLM_FI$df = length(unique(data$clusterID)) - length( fixef(M1) )

    formRIRC = make_regression_formula( FE = FALSE,
                                        control_formula = control_formula,
                                        cluster_RE = TRUE )
    formRIRC = update( formRIRC, . ~ . + (1+Z|blockID) )
    M_RIRC <- my_lmer( formRIRC, data=data )
    MLM_RIRC = get_mlm_ests(M_RIRC, "MLM_RIRC")


    formFIRC = make_regression_formula( FE = TRUE,
                                        control_formula = control_formula,
                                        cluster_RE = TRUE )
    formFIRC = update( formFIRC, . ~ 0 + . + (0+Z|blockID) )
    M_FIRC <- my_lmer( formFIRC, data=data )
    MLM_FIRC = get_mlm_ests(M_FIRC, "MLM_FIRC")


    res = bind_rows( MLM_FE,
                     MLM_RE,
                     MLM_FI,
                     MLM_RIRC,
                     MLM_FIRC )


    res


}




if ( FALSE ) {

    data( fakeCRT )

    fakeCRT

    formula =  Yobs ~ T.x | S.id | D.id
    control_formula = ~ X.jk + C.ijk
    data = clusterRCT:::make_canonical_data(formula=formula, data=fakeCRT,
                                            control_formula = control_formula)

    MLM_estimators(formula, control_formula = control_formula,
                   data = fakeCRT )

    MLM_estimators( Yobs ~ T.x | S.id | D.id,
                    data = fakeCRT )

    MLM_estimators( Yobs ~ T.x | S.id,
                    data = fakeCRT )

}



