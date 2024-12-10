


if ( FALSE ) {
    # make the characteristics table


    # datapasta::tribble_paste(a)

    if ( FALSE ) {
        cc %>% dplyr::select( method ) %>%
            mutate( weight = "person",
                    population = "super",
                    biased = 0 ) %>% datapasta::tribble_paste()
    }


    mc <- tibble::tribble(
        ~method, ~reg_weight,  ~weight,
        "LRi",   "", "Person",
        "LRi", "cw", "Cluster",
        "LRa","", "Cluster",
        "LRa","pw", "Person" )
    mc
    mc <- expand_grid( mc,
                       FE = c("", "FE", "FIpw", "FIcw", "FIbw" ),
                       SE = c("crve", "het", "db" ) ) %>%
        mutate( block_weight = ifelse( str_detect( FE, "FI" ),
                                       str_replace( FE, "FI", "" ),
                                       "" ) ) %>%
        filter( SE != "crve" | method != "LRa",
                SE != "het" | method != "LRi" )

    mc
    mc$biased = mc$FE == "FE"

    mc <- mutate( mc,
                  blocked = ifelse( FE == "", 0, 1 ),
                  method = paste0( method, reg_weight, "-", FE, "-", SE ),
                  method = stringr::str_replace( method, "--", "-"),
                  weight = make_weight_names( weight, block_weight ),
                  disfavored = weight %in% c( "Cluster-Person", "Person-Cluster") ) %>%
        dplyr::select( -reg_weight, -block_weight, -SE, -FE )

    mc$disfavored = as.numeric( mc$disfavored )
    mc$biased = as.numeric( mc$biased )
    print( mc, n=100 )

    others <- tibble::tribble(
        ~method,           ~weight, ~biased, ~blocked, ~disfavored,

        "MLM-FE",         "Cluster",       1,        1,           0,
        "MLM-FIRC",   "Cluster-Block",       1,        1,           1,
        "MLM-FIbw",   "Cluster-Block",       1,        1,           1,
        "MLM-FIcw",         "Cluster",       1,        1,           0,
        "MLM-FIpw",         "Cluster",       1,        1,           1,
        "MLM-RE",         "Cluster",       1,        1,           0,
        "MLM-RIRC",   "Cluster-Block",       1,        1,           1,

        "MLM",         "Cluster",       1,        0,           0,
        "DB_HT",          "Person",       0,        1,           1,
        "DB_Raj",          "Person",       0,        1,           1
    )
    others$population = NULL

    mc = bind_rows( mc, others ) %>%
        arrange( blocked, method )

    print( mc, n=100 )

    datapasta::tribble_paste( mc )


}

#' Get table of characteristics of all the methods implemented in this
#' package.
#'
#' @return Dataframe with columns for method, weight, population,
#'   biased, and disfavored (whether we do not like the estimator due
#'   to odd weighting or instability).
#' @importFrom dplyr bind_rows filter left_join mutate select
#' @examples
#' method_characteristics()
#'
#' @export
#'
method_characteristics <- function() {

    mc <-
        tibble::tribble(
            ~method,           ~weight, ~biased, ~blocked, ~disfavored,
            "LRa-db",         "Cluster",       0,        0,           0,
            "LRa-het",         "Cluster",       0,        0,           0,
            "LRapw-db",          "Person",       0,        0,           0,
            "LRapw-het",          "Person",       0,        0,           0,
            "LRi-crve",          "Person",       0,        0,           0,
            "LRi-db",          "Person",       0,        0,           0,
            "LRicw-crve",         "Cluster",       0,        0,           0,
            "LRicw-db",         "Cluster",       0,        0,           0,
            "MLM",         "Cluster",       1,        0,           0,
            "DB_HT",          "Person",       0,        1,           1,
            "DB_Raj",          "Person",       0,        1,           1,
            "LRa-FE-db",         "Cluster",       1,        1,           0,
            "LRa-FE-het",         "Cluster",       1,        1,           0,
            "LRa-FIbw-db",   "Cluster-Block",       0,        1,           1,
            "LRa-FIbw-het",   "Cluster-Block",       0,        1,           1,
            "LRa-FIcw-db", "Cluster-Cluster",       0,        1,           0,
            "LRa-FIcw-het", "Cluster-Cluster",       0,        1,           0,
            "LRa-FIpw-db",  "Cluster-Person",       0,        1,           1,
            "LRa-FIpw-het",  "Cluster-Person",       0,        1,           1,
            "LRapw-FE-db",          "Person",       1,        1,           0,
            "LRapw-FE-het",          "Person",       1,        1,           0,
            "LRapw-FIbw-db",    "Person-Block",       0,        1,           1,
            "LRapw-FIbw-het",    "Person-Block",       0,        1,           1,
            "LRapw-FIcw-db",  "Person-Cluster",       0,        1,           1,
            "LRapw-FIcw-het",  "Person-Cluster",       0,        1,           1,
            "LRapw-FIpw-db",   "Person-Person",       0,        1,           0,
            "LRapw-FIpw-het",   "Person-Person",       0,        1,           0,
            "LRi-FE-crve",          "Person",       1,        1,           0,
            "LRi-FE-db",          "Person",       1,        1,           0,
            "LRi-FIbw-crve",    "Person-Block",       0,        1,           1,
            "LRi-FIbw-db",    "Person-Block",       0,        1,           1,
            "LRi-FIcw-crve",  "Person-Cluster",       0,        1,           1,
            "LRi-FIcw-db",  "Person-Cluster",       0,        1,           1,
            "LRi-FIpw-crve",   "Person-Person",       0,        1,           0,
            "LRi-FIpw-db",   "Person-Person",       0,        1,           0,
            "LRicw-FE-crve",         "Cluster",       1,        1,           0,
            "LRicw-FE-db",         "Cluster",       1,        1,           0,
            "LRicw-FIbw-crve",   "Cluster-Block",       0,        1,           1,
            "LRicw-FIbw-db",   "Cluster-Block",       0,        1,           1,
            "LRicw-FIcw-crve", "Cluster-Cluster",       0,        1,           0,
            "LRicw-FIcw-db", "Cluster-Cluster",       0,        1,           0,
            "LRicw-FIpw-crve",  "Cluster-Person",       0,        1,           1,
            "LRicw-FIpw-db",  "Cluster-Person",       0,        1,           1,
            "MLM-FE",         "Cluster",       1,        1,           0,
            "MLM-FIRC",   "Cluster-Block",       1,        1,           1,
            "MLM-FIbw",   "Cluster-Block",       1,        1,           1,
            "MLM-FIcw",         "Cluster",       1,        1,           0,
            "MLM-FIpw",         "Cluster",       1,        1,           1,
            "MLM-RE",         "Cluster",       1,        1,           0,
            "MLM-RIRC",   "Cluster-Block",       1,        1,           1
        )


    # Should we let estimators define that in the guts of the code?
    # Many of the estimators do.
    mc$weight = NULL

    return( mc )
}



#' Compare different estimators for a (possibly blocked) cluster
#' randomized trial.
#'
#' This function calculates the point estimates and SE estimates for a
#' variety of the estimators implemented by this package.
#'
#' @param Yobs vector observed outcomes (or column name in data)
#' @param Z vector of assignment indicators (1==treated) (or column
#'   name in data)
#' @param B vector of block ids (or column name in data)
#' @param blockID block ids (variable name as string if data frame
#'   passed) (if randomization blocks are nested in block).
#' @param data frame holding Y, Z, B and (possibly a column with name
#'   specified by blockID).
#' @param include_MLM Include MLM estimators
#' @param include_DB Include Design-Based estimators (taken from
#'   RCTYes documentation and prior literature).
#' @param include_LM Include Linear Model-based estimators (including
#'   Huber-White SEs, etc.)
#' @param include_agg Include estimators applied to aggregated data,
#'   aggregated at the cluster level.
#' @param include_dumb Include "dumb" estimators (i.e., those
#'   interacted estimators that weight by both person and cluster or
#'   vice versa).
#' @param control_formula What variables to control for, in the form
#'   of "~ X1 + X2".
#' @param patch_data If TRUE impute all missing covariates with mean
#'   imputation (adding dummy variables as needed) via the
#'   `patch_data()` method.  Will drop all rows with missing outcome,
#'   treatment, or clustering info.  If FALSE do not do this.
#' @param include_method_characteristics Include details of the
#'   methods (target estimands and sampling framework assumed) in the
#'   return value.
#'
#' @return Dataframe of point estimates and standard errors for each
#'   method considered. If \code{include_method_characteristics=TRUE}
#'   also add some features of the methods as additional columns.
#' @export
compare_methods <- function(formula,
                            data = NULL,
                            control_formula = NULL,
                            weight = c( "individual", "cluster" ),
                            include_MLM = TRUE,
                            include_DB = TRUE,
                            include_LM = TRUE,
                            include_agg = TRUE,
                            warn_missing = TRUE,
                            include_dumb = FALSE,
                            include_disfavored = FALSE,
                            patch_data = TRUE,
                            handle_singleton_blocks = c( "drop", "pool", "fail" ),
                            include_method_characteristics = TRUE ) {

    handle_singleton_blocks = match.arg(handle_singleton_blocks)

    if ( !is.null( formula ) ) {

        if ( !inherits( formula, "formula" ) ) {
            stop( "formula must be a formula object" )
        }

        lhs_vars <- formula.tools::lhs.vars(formula)
        if (length(lhs_vars) > 1) {
            res <- purrr::map(
                lhs_vars,
                function(v) {
                    # see https://stackoverflow.com/questions/37824625/how-to-modify-the-left-side-of-a-formula
                    outcome <- as.name(v)
                    formula1 <- as.formula(bquote(.(outcome) ~ .(formula.tools::rhs(formula))))
                    res <- compare_methods(
                        formula1,
                        data = data,
                        control_formula = control_formula,
                        weight = weight,
                        include_MLM = include_MLM,
                        include_DB = include_DB,
                        include_LM = include_LM,
                        include_agg = include_agg,
                        include_method_characteristics = include_method_characteristics,
                        patch_data = patch_data,
                        handle_singleton_blocks = handle_singleton_blocks,
                        warn_missing = warn_missing ) %>%
                        dplyr::mutate(outcome = v, .before=method)
                }) %>%
                purrr::list_rbind()
            return(res)
        }


        data = make_canonical_data( formula=formula, data=data,
                                    control_formula = control_formula,
                                    warn_missing = warn_missing,
                                    patch_data = patch_data )
        if ( patch_data ) {
            control_formula = attr( data, "control_formula" )
        }

    } else if ( patch_data ) {
        data = patch_data_set( NULL, data, control_formula = control_formula,
                               warn_missing = warn_missing )
        control_formula = attr( data, "control_formula" )
    }


    # handle blocks with singletons
    if ( ("blockID" %in% names(data)) && handle_singleton_blocks != "fail" ) {
        data = patch_singleton_blocks( Yobs ~ Z | clusterID | blockID, data=data,
                                       drop_data = handle_singleton_blocks == "drop",
                                       warn_missing = warn_missing )
    }


    n <- nrow( data )
    check_data_integrity( data )


    # aggregate data once for efficiency.
    aggdat = aggregate_data( data = data, control_formula = control_formula)
    control_formula_agg = attr( aggdat, "control_formula" )


    summary_table <- data.frame()

    if ( "blockID" %in% colnames(data) && length( unique( data$blockID ) ) == 1 ) {
        data$blockID = NULL
    }
    has_blocks <- "blockID" %in% colnames(data)


    if (include_LM) {
        lms <- linear_model_estimators(formula = NULL, data = data,
                                       control_formula = control_formula )
        lmsW <- linear_model_estimators(formula = NULL, data = data,
                                        weight = "Cluster",
                                        control_formula = control_formula )
        summary_table = dplyr::bind_rows(summary_table, lms, lmsW)

        if ( has_blocks ) {
            lms_int = interacted_linear_model_estimators( formula = NULL, data = data,
                                                          control_formula = control_formula )
            lms_intW = interacted_linear_model_estimators( formula = NULL, data = data,
                                                           weight = "Cluster",
                                                           control_formula = control_formula )

            summary_table = dplyr::bind_rows(summary_table,
                                             lms_int, lms_intW)
        }
    }


    if (include_MLM) {
        mlms <- MLM_estimators(formula = NULL, data = data,
                               control_formula = control_formula,
                               include_disfavored = include_disfavored )
        summary_table <- dplyr::bind_rows( summary_table, mlms )
    }


    if (include_agg) {
        lms <- aggregation_estimators( formula = NULL, data = aggdat,
                                       control_formula = control_formula_agg,
                                       aggregated = TRUE )
        summary_table = dplyr::bind_rows(summary_table, lms)
    }

    if (include_DB) {
        db_res_i <- design_based_estimators(formula = NULL, data = aggdat,
                                            control_formula = control_formula_agg,
                                            weight = "Person",
                                            aggregated = TRUE)
        db_res_c <- design_based_estimators(formula = NULL, data = aggdat,
                                            control_formula = control_formula_agg,
                                            weight = "Cluster",
                                            aggregated = TRUE)

        db_res_ii <- design_based_estimators_individual(formula = NULL, data = data,
                                                        control_formula = control_formula,
                                                        weight = "Person")
        db_res_ci <- design_based_estimators_individual(formula = NULL, data = data,
                                                        control_formula = control_formula,
                                                        weight = "Cluster")
        summary_table = dplyr::bind_rows( summary_table,
                                          db_res_i, db_res_c,
                                          db_res_ii, db_res_ci )

        if ( include_disfavored ) {
            db_middleton <- middleton_aronow_estimator(formula = NULL, data = aggdat,
                                                       control_formula = control_formula_agg,
                                                       aggregated = TRUE)

            summary_table = dplyr::bind_rows( summary_table, db_middleton )
        }
    }

    if ( !include_dumb ) {
        summary_table <- summary_table %>%
            filter( !grepl( "Cluster-Person|Person-Cluster", weight ) )
    }

    mc <- method_characteristics()
    mc$blocked = NULL
    mc$weight = NULL
    summary_table <- left_join( summary_table, mc, by = "method" )

    if ( !include_disfavored ) {
        summary_table <- filter( summary_table, disfavored == 0 )
    }
    summary_table$weight = stringr::str_replace( summary_table$weight, "Cluster-Cluster", "Cluster" )
    summary_table$weight = stringr::str_replace( summary_table$weight, "Person-Person", "Person" )

    summary_table <- summary_table %>%
        relocate( method, weight )

    # If not desired, remove info on the methods (e.g., what estimand they are targeting)
    if (!include_method_characteristics) {
        summary_table <- summary_table %>%
            dplyr::select( -weight, -biased, -disfavored )
    }

    if ( nrow( summary_table ) > 0 ) {
        summary_table = tibble::remove_rownames( summary_table )
    }

    summary_table <- summary_table %>%
        arrange( method )

    return(summary_table)
}




# Testing code ---

if ( FALSE ) {

    library( tidyverse )
    library( clusterRCT )
    data( fakeCRT )
    fakeCRT

    design_based_estimators_individual( Yobs ~ T.x | S.id | D.id,
                                        control_formula = ~ X.jk + C.ijk,
                                        data = fakeCRT, weight="Cluster" )


    compare_methods( Yobs ~ T.x | S.id | D.id,
                     control_formula = ~ X.jk + C.ijk,
                     include_dumb = FALSE,
                     data = fakeCRT ) %>%
        knitr::kable( digits = 2 )

}
