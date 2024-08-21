


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

    # datapasta::tribble_paste(a)

    if ( FALSE ) {
        cc %>% dplyr::select( method ) %>%
            mutate( weight = "person",
                    population = "super",
                    biased = 0 ) %>% datapasta::tribble_paste()
    }

    mc <- tibble::tribble(
        ~method,   ~weight, ~population, ~biased,  ~disfavored,
        "LR_FE_CRVE", "person",          "super",       0,  0,
        "LR_FI_CRVE_Block", "person/block",         "super",       1, 0,
        "LR_FI_CRVE_Cluster", "person",         "super",       0,1,
        "LR_FI_CRVE_Person", "person",          "super",       0,0,
        "MLM_FE", "cluster",          "super",       1,0,
        "MLM_RE", "cluster",          "super",       1,0,
        "MLM_FI_Block", "cluster/block",         "super",       1,0,
        "MLM_FI_Cluster", "cluster",         "super",       1,0,
        "MLM_FI_Person", "cluster",          "super",       1,1,
        "MLM_RIRC", "cluster/block",         "super",       1,0,
        "MLM_FIRC", "cluster/block",         "super",       1,0,
        "Agg_FE_Cluster", "cluster",         "super",       0,0,
        "Agg_FE_Person", "person",          "super",       0,0,
        "Agg_FI_Cluster_Block", "cluster/block",        "super",       1,0,
        "Agg_FI_Cluster_Cluster", "cluster",         "super",       0,0,
        "Agg_FI_Person_Block", "person/block",         "super",       1,0,
        "Agg_FI_Person_Person", "person",          "super",       0,0,
        "DB_FI_Person_Block", "person/block",         "super",       0,0,
        "DB_FI_Person_Person", "person",          "super",       0,0,
        "DB_FE_Person", "person",         "super",       1,0,
        "DB_FI_Cluster_Block", "cluster/block",         "super",       1,0,
        "DB_FI_Cluster_Cluster", "cluster",         "super",       0,0,
        "DB_FE_Cluster", "cluster",         "super",       0,0,
        "DB_HT", "person",         "super",       0,1,
        "DB_Raj", "person",         "super",       0,1,

        # Methods that are disfavored due to odd weighting
        "DB_FI_Cluster_Person", "cluster",         "super",       1,1,
        "DB_FI_Person_Cluster", "person",         "super",       1,1,
        "Agg_FI_Cluster_Person", "cluster",         "super",       1,1,
        "Agg_FI_Person_Cluster", "person",         "super",       1,1,

        # Non-blocked estimators
        "LR_CRVE", "person",      "super",       0,0,
        "MLM", "cluster",     "super",       1,0,
        "Agg_Cluster", "cluster",     "super",       0,0,
        "Agg_Person", "person",      "super",       0,0,
        "DB_Person", "person",      "super",       0,0,
        "DB_Cluster", "cluster",     "super",       0,0,
    )


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
                            patch_data = TRUE,
                            handle_singleton_blocks = c( "drop", "pool", "fail" ),
                            include_method_characteristics = TRUE ) {

    handle_singleton_blocks = match.arg(handle_singleton_blocks)

    if ( !is.null( formula ) ) {

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


    # handle singleton blocks
    if ( ("blockID" %in% names(data)) && handle_singleton_blocks != "fail" ) {
        data = patch_singleton_blocks( Yobs ~ Z | clusterID | blockID, data=data,
                                       drop_data = handle_singleton_blocks == "drop",
                                       warn_missing = warn_missing )
    }


    n <- nrow( data )
    check_data_integrity( data )


    # aggregate data once for efficiency.
    aggdat = aggregate_data(data, control_formula)
    control_formula_agg = attr( aggdat, "control_formula" )


    summary_table <- data.frame()

    has_blocks <- "blockID" %in% colnames(data)

    if (include_LM) {
        lms <- linear_model_estimators(formula = NULL, data = data,
                                       control_formula = control_formula )
        summary_table = dplyr::bind_rows(summary_table, lms)

        if ( has_blocks ) {
            lms_int = interacted_linear_model_estimators( formula = NULL, data = data,
                                                          control_formula = control_formula )
            summary_table = dplyr::bind_rows(summary_table, lms_int)
        }
    }


    if (include_MLM) {
        mlms <- MLM_estimators(formula = NULL, data = data,
                               control_formula = control_formula )
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

        db_middleton <- middleton_aronow_estimator(formula = NULL, data = aggdat,
                                                   control_formula = control_formula_agg,
                                                   aggregated = TRUE)
        summary_table = dplyr::bind_rows( summary_table, db_res_i, db_res_c, db_middleton )
    }

    if ( !include_dumb ) {
        summary_table <- summary_table %>%
            filter( !grepl( "Person_Cluster", method ),
                    !grepl( "Cluster_Person", method ),
                    !grepl( "LR_FI_CRVE_Cluster|MLM_FI_Person", method ) )
    }

    # Add info on the methods (e.g., what estimand they are targeting)
    summary_table$weight = NULL
    if (include_method_characteristics) {
        mc <- method_characteristics()

        summary_table <- left_join( summary_table, mc, by = "method" )
    }

    if ( nrow( summary_table ) > 0 ) {
        summary_table = tibble::remove_rownames( summary_table )
    }

    return(summary_table)
}
