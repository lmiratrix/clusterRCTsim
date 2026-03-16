
# Library of functions to get data ready and then compare or describe
# data.  Can take a number of control covariates, and select the
# highest correlation covariate.

# install the clusterRCT package from GitHub
# devtools::install_github("https://github.com/lmiratrix/clusterRCT")

library(tidyverse)
library(haven)
library(tools)
library(writexl)
library(clusterRCT)


#' Get the highest correlated covariate to a given outcome.
#'
#' @return List of two items, the name of the best covariate and the
#'   correlation of that covariate to the outcome.
get_best_covariate <- function( covars_list, outcome, data ) {
  
  # Calculate correlations of covariate to outcome
  correlations <- sapply(covars_list, function(covariate) {
    if ( is.factor(data[[covariate]]) || is.character(data[[covariate]]) ) {
      return(0)
    }
    if ( sd( data[[covariate]], na.rm=TRUE ) == 0 ) {
      # CH - this line is causing an error - I think because "proj" does not exist
      # cli::cli_alert( "Zero variance in covariate {covariate} for {proj}/{outcome}" )
      cli::cli_alert( "Zero variance in covariate {covariate} for {outcome}" )
      return(0)
    }
    cor(data[[outcome]], data[[covariate]], use = "complete.obs")
  })
  
  # Find the covariate with the highest correlation
  max_cor = max(abs(correlations), na.rm=TRUE)
  best_covariate = NA
  if ( is.na( max_cor ) ) {
    cli::cli_alert( "No max correlation found and thus no best covariate found" )
  } else {
    best_covariate <- names(which.max(abs(correlations)))
  }
  return( list( best_covariate = best_covariate, max_cor = max_cor ) )
}





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
      cli::cli_alert( "No best covariate found for {proj}/{outcome}. Skipping dataset." )
      return( NULL )
    }
    
    cf <- as.formula(paste('~', paste(bc$best_covariate, collapse = ' + ')))
    
    # Now call the actual function passed
    results <- data %>%
      processor(form,
                data = .,
                control_formula = cf, ... )
    
    results$best_cov = bc$best_covariate
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
cm <- function(proj, outcome, covars, covars_list, data, trt, clust, blk, ptest, include_gee) {
  
  data_processor(proj=proj, outcome=outcome, covars=covars, covars_list=covars_list, 
                 data=data, trt=trt, clust=clust, blk=blk, ptest=ptest, include_gee=FALSE,
                 processor = compare_methods, 
                 include_method_characteristics = TRUE,
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




