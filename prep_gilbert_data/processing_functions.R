
# Processing functions

library( tidyverse )


#' i is the ID of the dataset to read in (if numeric).  If i is a
#' string, it is the path to the file
read_in <- function(i) {
    if ( is.numeric(i) ) {
        dat <- read_csv(glue::glue(here::here("prep_gilbert_data/data/raw/crt_data_{i}.csv")),
                        show_col_types = FALSE)
    } else {
        dat <- read_csv(i, show_col_types = FALSE)
    }
    #dat$source = i
    dat
}



#' Do some quick standardization on given dataset and produce list of
#' all covariates. Also identify block ID (or NULL if not blocked)
#'
#' NOTE: This will look for global variable "dataset_list" to pull a
#' dataset from, if the dataset is not directly passed.
#'
#' @param i If not a dataset, then pull this from the passed dataset
#'   list of datasets.
grab_data <- function( i, i_name, dataset_list = NULL, size_threshold = NULL ) {
    if ( is.data.frame(i) ) {
        cur_dat = i
        i = i_name
    } else {
        if ( is.numeric(i) ) {
            i = as.character(i)
        }
        stopifnot( !is.null( dataset_list ) )

        cur_dat = dataset_list[[i]]
        i_name = as.character(i)
    }

    cli::cli_alert_info( "Working on {i} ({nrow(cur_dat)} rows)" )

    # NOTE: Some datasets didn't have cluster ID either?  Simple randomization?
    if ( !( "cluster_id" %in% colnames(cur_dat) ) ) {
        cli::cli_alert( glue::glue( "No cluster_id in dataset {i} ({nrow(cur_dat)} rows)") )
        return ( NULL )
    }

    if ( !is.null( size_threshold ) && nrow( cur_dat ) > size_threshold ) {
        cli::cli_alert( glue::glue( "Skipping {i} since it is SO BIG ({nrow(cur_dat)} rows)--do this for real in future" ) )
        return ( NULL )
    }

    if ( all( is.na( cur_dat$y_std ) ) ) {
        # standardize according to control group mean and SD
        control_desc <- cur_dat |>
            filter(treat == 0) |>
            summarise(sd = sd(y, na.rm=TRUE),
                      mean = mean(y,na.rm=TRUE))

        cur_dat <- cur_dat |>
            mutate(y_std = (y - control_desc[[1,2]])/(control_desc[[1,1]]))
    }

    # get covariates list
    cov_list <- cur_dat |>
        select(std_baseline, starts_with("cov_")) |>
        colnames()

    cur_dat <- filter( cur_dat, !is.na( cluster_id ) )

    blk = "block_id"
    if ( !(blk %in% colnames(cur_dat) )) {
        blk = NULL
    } else {
        cur_dat$cluster_id = paste0( cur_dat$block_id, "_", cur_dat$cluster_id )
    }


    return( list( name = i_name, data = cur_dat, cov_list = cov_list, blk = blk ) )
}



# load libraries ----


# install the clusterRCT package from GitHub
# devtools::install_github("https://github.com/lmiratrix/clusterRCT")

library(tidyverse)
library(haven)
library(tools)
library(writexl)
library(clusterRCT)


# Helper functions for compiling results ----



get_best_covariate <- function( covars_list, outcome, data ) {

    # Calculate correlations of covariate to outcome
    correlations <- sapply(covars_list, function(covariate) {
        if ( is.factor(data[[covariate]]) || is.character(data[[covariate]]) ) {
            return(0)
        }
        if ( sd( data[[covariate]], na.rm=TRUE ) == 0 ) {
            cli::cli_alert( "Zero variance in covariate {covariate} for {proj}/{outcome}" )
            return(0)
        }
        cor(data[[outcome]], data[[covariate]], use = "complete.obs")
    })

    # Find the covariate with the highest correlation
    max_cor = max(abs(correlations), na.rm=TRUE)
    best_covariate = NA
    if ( is.na( max_cor ) ) {
        cli::cli_alert( "No max correlation and thus no best covariate found" )
    } else {
        best_covariate <- names(which.max(abs(correlations)))
    }
    return( list( best_covariate = best_covariate, max_cor = max_cor ) )
}


# Functions to process data (compare or describe) ----


# the functions below run compare_methods() and describe_clusterRCT()
# for a given outcome and set of covariates. writing the functions in
# this generalized way allows the code below to be streamlined and
# avoid repetition

#' compare_methods() helper function for data processing
#'
#' Call compare_methods after doing some selection of variables and
#' whatnot given some overall info on the dataset.
#'
#' @param proj project name -- used to ID different datasets
#' @param outcome outcome variable as string
#' @param covars character string indicating which set to include.
#'   none is no covariate adjustment.  all is all in the covars_list.
#'   "best" selects the highest correlated covariate.
#' @param covars_list list of covariates to include in the model
#' @param data dataset to use (data frame)
#' @param trt treatment indicator as string
#' @param clust cluster variable as string
#' @param blk block variable as string.  NULL means no blocking.
#' @param ptest Existence of pretest (0/1).  Note, this is just passed
#'   to the final result.  We select the highest correlated covariate,
#'   and do not record whether it was the pretest or not.
#' @param processor A function to apply to the final dataset
#' @param ... Extra arguments to pass to processor function.
data_processor <- function(proj, outcome, covars, covars_list, data, trt,
                           clust, blk, ptest,
                           processor = compare_methods, ... ) {

    # rename outcome and treatment indicator
    data <- data %>%
        mutate(Y = !!sym(outcome),
               Tx = !!sym(trt))

    # EMERGE had multiple treatment arms, so we need to filter appropriately
    if (proj == 'EMERGE') {
        data <- data %>%
            filter(experimentalgroup_desc %in% c(str_remove(trt, 'T_'), 'Control'))
    }

    # build control formula
    if ( is.null(blk) ) {
        form <- as.formula(paste0('Y ~ ', paste(c('Tx', clust), collapse = ' | ')))
    } else {
        form <- as.formula(paste0('Y ~ ', paste(c('Tx', clust, blk), collapse = ' | ')))
    }

    if ( covars != "none" && length( covars_list ) == 0 ) {
        cli::cli_alert( "No covariates to include in {proj}/{outcome}. Changing to 'none'" )
        covars = 'none'
    }

    if (covars == 'none') {
        results <- data %>%
            processor(form,
                      data = ., ... )
    } else if (covars == 'all') {
        cf <- as.formula(paste('~', paste(covars_list, collapse = ' + ')))

        results <- data %>%
            processor(form,
                      data = .,
                      control_formula = cf, ... )
    } else if (covars == 'best') {


        bc <- get_best_covariate( covars_list, outcome = outcome, data = data )
        if ( is.na( bc$best_covariate ) ) {
            cli::cli_alert( "No best covariate found for {proj}/{outcome}" )
            return( NULL )
        }

        cf <- as.formula(paste('~', paste(bc$best_covariate, collapse = ' + ')))

        results <- data %>%
            processor(form,
                      data = .,
                      control_formula = cf, ... )
        results$max_cor = bc$max_cor
    }

    o_desc = attr(data[[outcome]], 'label')
    if (is.null(o_desc)) {
        o_desc = 'No description'
    }

    # spit out results with project and outcome information
    results <- results %>%
        mutate(Project      = proj,
               Pretest      = ptest,
               Covariates   = covars,
               Trt_arm      = ifelse(proj == 'EMERGE', str_remove(trt, 'T_'), NA),
               Outcome      = outcome,
               Outcome_desc = o_desc)

    return(results)
}


#' compare_methods() helper function for data processing
#'
#' Call compare_methods after doing some selection of variables and
#' whatnot given some overall info on the dataset.
#'
#' @param proj project name -- used to ID different datasets
#' @param outcome outcome variable as string
#' @param covars character string indicating which set to include.
#'   none is no covariate adjustment.  all is all in the covars_list.
#'   "best" selects the highest correlated covariate.
#' @param covars_list list of covariates to include in the model
#' @param data dataset to use (data frame)
#' @param trt treatment indicator as string
#' @param clust cluster variable as string
#' @param blk block variable as string
#' @param ptest Existence of pretest (0/1).  Note, this is just passed
#'   to the final result.  We select the highest correlated covariate,
#'   and do not record whether it was the pretest or not.
cm <- function(proj, outcome, covars, covars_list, data, trt, clust, blk, ptest) {

    data_processor(proj=proj, outcome=outcome, covars=covars, covars_list=covars_list,
                   data=data, trt=trt, clust=clust, blk=blk, ptest=ptest,
                   processor = compare_methods,
                   include_method_characteristics = TRUE,
                   include_gee = FALSE,
                   warn_missing = FALSE )

}


#' describe_clusterRCT helper function
#'
#' Call describe_clusterRCT after doing some selection of variables and
#' whatnot given some overall info on the dataset.
#'
#' See cm, above, for description of parameters.
dc <- function(proj, outcome, covars, covars_list, data, trt, clust, blk, ptest) {

    data_processor(proj=proj, outcome=outcome, covars=covars, covars_list=covars_list,
                   data=data, trt=trt, clust=clust, blk=blk, ptest=ptest,
                   processor = describe_clusterRCT,
                   warn_missing = FALSE )
}


# Get indexed dataset and run compare_methods for no, best, all
# covariates
#
# Return result as dataframe.
compare_all_cov <- function(i, size_threshold = NULL, dataset_list = NULL ) {

    dta = grab_data(i, deparse( substitute(i) ),
                    size_threshold = size_threshold,
                    dataset_list = dataset_list )

    if ( is.null(dta) ) {
        return( NULL )
    }

    m1 <- cm(glue::glue("jg_crt_{dta$name}"), "y_std", "none",
             dta$cov_list, dta$data, "treat", "cluster_id", dta$blk, 0)
    m2 <- cm(glue::glue("jg_crt_{dta$name}"), "y_std", "all",
             dta$cov_list, dta$data, "treat", "cluster_id", dta$blk, 0)
    m3 <- cm(glue::glue("jg_crt_{dta$name}"), "y_std", "best",
             dta$cov_list, dta$data, "treat", "cluster_id", dta$blk, 0)

    out <- bind_rows(m1,m2,m3)

    out$blocked = as.numeric( !is.null(dta$blk) )

    out$id = dta$name

    return(out)
}


describe_all_cov <- function(i, dataset_list = NULL ){

    dta = grab_data( i, deparse( substitute(i) ), dataset_list = dataset_list )
    if ( is.null(dta) ) {
        return( NULL )
    }

    m0 <- dc(proj = glue::glue( "jg_crt_{dta$name}"), outcome = "y_std", covars="none",
             covars_list=dta$cov_list, data=dta$data,
             trt="treat", clust="cluster_id", blk=dta$blk, ptest=0)

    m2 <- dc(proj = glue::glue( "jg_crt_{dta$name}"), outcome = "y_std", covars="all",
             covars_list=dta$cov_list, data=dta$data,
             trt="treat", clust="cluster_id", blk=dta$blk, ptest=0)
    #return( m2 )

    # Now 'best' covariate
    m1 <- dc(proj = glue::glue( "jg_crt_{dta$name}"), outcome = "y_std", covars = "best",
             covars_list=dta$cov_list, data=dta$data,
             trt="treat", clust="cluster_id", blk=dta$blk, ptest=0)

    # trt = "treat", "std_baseline", dataset_list[[i]], "treat", "cluster_id", "block_id", 0)

    out <- bind_rows(m0, m1, m2)
    return(out)

}


# Testing and verification code ----


if ( FALSE ) {
    # rename outcome and treatment indicator
    data <- data %>%
        mutate(Y = !!sym(outcome),
               Tx = !!sym(trt))

    # EMERGE had multiple treatment arms, so we need to filter appropriately
    if (proj == 'EMERGE') {
        data <- data %>%
            filter(experimentalgroup_desc %in% c(str_remove(trt, 'T_'), 'Control'))
    }
    # build control formula
    form <- as.formula(paste0('Y ~ ', paste(c('Tx', clust, blk), collapse = ' | ')))

    if (covars == 'none') {
        results <- data %>%
            describe_clusterRCT(form,
                                data = .)
    } else if (covars == 'all') {
        cf <- as.formula(paste('~', paste(covars_list, collapse = ' + ')))

        results <- data %>%
            describe_clusterRCT(form,
                                data = .,
                                control_formula = cf)
    } else if (covars == 'best') {
        # Calculate correlations
        correlations <- sapply(covars_list, function(covariate) {
            cor(data[[outcome]], data[[covariate]], use = "complete.obs")
        })

        # Find the covariate with the highest correlation
        best_covariate <- names(which.max(abs(correlations)))

        cf <- as.formula(paste('~', paste(best_covariate, collapse = ' + ')))

        results <- data %>%
            describe_clusterRCT(form,
                                data = .,
                                control_formula = cf)
    }

    o_desc = attr(data[[outcome]], 'label')
    if (is.null(o_desc)) {
        o_desc = 'No description'
    }

    # spit out results with project and outcome info
    results <- results %>%
        mutate(Project      = proj,
               Pretest      = ptest,
               Covariates   = covars,
               Trt_arm      = ifelse(proj == 'EMERGE', str_remove(trt, 'T_'), NA),
               Outcome      = outcome,
               Outcome_desc = o_desc )

    return(results)
}


if ( FALSE ) {

    compare_all_cov( 5, dataset_list = dfs )

    # Testing (need to load data first)
    describe_all_cov( 1, dataset_list = dfs )
    describe_all_cov( 2, dataset_list = dfs )
}
