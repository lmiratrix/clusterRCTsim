

# Testing the fingerprint code

# This is a version of master_run_fingerprints that fleshes out the
# sim on the fakeCRT to see if covariates are actually included.
# (They are.)

# load libraries ----
# devtools::install_github("lmiratrix/clusterRCT")
library(tidyverse)
library(haven)
library(clusterRCT)
library(tools)
library(writexl)


source( here::here( "empirical/data_processing_functions.R" ) )

# Where to save result files (one file per outcome will be generated)
RESULT_DIR = here::here( "rerand_sim_results" )


# How many iterations for each outcome?
NUM_SIMULATIONS = 50

# What configuration to we run our simulation over?
COVAR_SETS = c('none', 'best')


# Make place for results, if needed
dir.create( RESULT_DIR, showWarnings = FALSE )




# Wrapper function to run the simulations ----

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
run_sim <- function(proj, outcome, covars, covars_list, data, trt, clust, blk, ptest, 
                    R = NUM_SIMULATIONS,
                    ... ) {
  
  res <- data_processor( proj=proj, outcome=outcome, covars=covars, covars_list=covars_list, 
                         data=data, trt=trt, clust=clust, blk=blk, ptest=ptest, 
                         processor = run_rerandomize_simulation, 
                         include_method_characteristics = FALSE,
                         warn_missing = FALSE,
                         R = R,
                         ... )
  
  res$Outcome_desc = NULL
  res$Pretest = NULL
  
  # This is the same ID as built in the bottom of 'master_analyze_all_data.R'
  UniqueID     = str_remove_all(paste(proj, outcome, trt, sep = '_'), '_NA')
  
  write_rds( res, 
             here::here( glue::glue( "{RESULT_DIR}/{UniqueID}.rds" ) ) )
  res
  
}



# Make some fake data ----


data( fakeCRT )
fakeCRT$Xg = rnorm( nrow(fakeCRT), mean=0, sd=1 )
fakeCRT$Y2 = fakeCRT$Xg + rnorm( nrow(fakeCRT), sd=0.3 )
fakeCRT$XX = rnorm( nrow(fakeCRT) )

outcomes_fake <- c('Yobs', 'Y2')

covars_fake <- c('X.jk', 'C.ijk', "XX", "Xg" )


# Run the simulation ----


# Direct call (for testing)
set.seed( 40404 )
rc <- run_sim( proj = "FakeCRT",
               outcome = "Y2",
               covars = "best",
               covars_list = covars_fake,
               data = fakeCRT,
               trt = "T.x",
               clust = "S.id",
               blk = "D.id",
               ptest = 0,
               R = NUM_SIMULATIONS,
               include_agg = FALSE,
               include_DB = FALSE,
               include_MLM = FALSE )
table(rc$method )
head( rc )

set.seed( 40404 )
rc2 <- run_sim( proj = "FakeCRT",
                outcome = "Y2",
                covars = "none",
                covars_list = covars_fake,
                data = fakeCRT,
                trt = "T.x",
                clust = "S.id",
                blk = "D.id",
                ptest = 0,
                R = NUM_SIMULATIONS,
                include_agg = FALSE,
                include_DB = FALSE,
                include_MLM = FALSE )
head( rc2 )


# Analyze the results ----

summary( rc$ATE_hat - rc2$ATE_hat )

compare_fake = bind_rows( best=rc, none=rc2, .id = "Covariates" )

# making sure we get different answers based on covariates
cc = filter( compare_fake, runID == "empirical", Outcome =="Y2" )
nrow( cc )
ccW = cc %>%
  dplyr::select( method, Covariates, ATE_hat ) %>%
  pivot_wider( names_from = Covariates, values_from = ATE_hat ) 
ccW

compare_fake
ggplot( compare_fake, aes( ATE_hat ) ) +
  geom_dotplot() +
  facet_grid( Covariates ~ method ) +
  theme( legend.position = "bottom" )

compare_fake %>%
  group_by( method, Covariates ) %>%
  summarise( sd = sd( ATE_hat ) ) %>%
  pivot_wider( names_from = Covariates, values_from = sd )



# Just checking what our results look like (unneeded code):
compare_fake <- as_tibble( compare_fake )
compare_fake

compare_fake$weight = get_estimand( compare_fake$method, simple = TRUE )
agg <- summarize_simulation_results( compare_fake, 
                                     summarize_results = "cross" )
agg$real = agg$runID == "empirical"
ggplot( agg, aes( rangeATE, fill=real ) ) +
  facet_wrap( ~ weight ) +
  geom_dotplot()

