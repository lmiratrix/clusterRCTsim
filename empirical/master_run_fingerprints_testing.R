
# This is a verion of the master fingerprint that runs on a fake
# dataset to see if the standardization is broken.

# load libraries ----
# devtools::install_github("lmiratrix/clusterRCT")
library(tidyverse)
library(haven)
library(clusterRCT)
library(tools)
library(writexl)
library( clusterRCT )

source( here::here( "empirical/data_processing_functions.R" ) )

if ( FALSE ) {
  # make fake_data folder and put some fake data in it with a different scale
  dir.create( here::here( "empirical/fake_data/" ), showWarnings = FALSE )
  
  
  studies = read_csv( here::here( "empirical/results/study_descriptions.csv" ),
                      show_col_types = FALSE )
  studies
  
  stud = filter( studies, UniqueID == "CARES_TR_BPI_Incred" )
  clusterRCT:::print.clusterRCTstats( stud )
  stud$ncov.2 = 4
  source( "R/generate_data_from_fingerprint.R" )
  stud$sdY0 = 150
  dat = generate_data_from_fingerprint( structure = NULL,
                                        parameters = stud )  
  
  dat
  sd( dat$Y0 )
  stud$sdY0
  
  saveRDS( dat, here::here( "empirical/fake_data/fake_data.rds" ) )
  
}


# Where the prepped data is
DATA_DIR = here::here( "empirical/fake_data/" )


# Where to save result files (one file per outcome will be generated)
RESULT_DIR = here::here( "empirical/fake_data_results/" )


# How many iterations for each outcome?
NUM_SIMULATIONS = 20

# What configuration to we run our simulation over?
COVAR_SETS = c('none', 'best', 'all')

# Run in parallel?
RUN_PARALLEL = FALSE

# Make place for results, if needed
dir.create( RESULT_DIR, showWarnings = FALSE )



# load cleaned data ----

# these are individual datasets for each non-clean room study, already
# filtered to analytic sample
#
# This code sets a bunch of environment variables to the given
# filenames
for (f in list.files(DATA_DIR, full.names = T)) {
  assign(file_path_sans_ext(basename(f)),
         readRDS(f))
}





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
                         parallel = RUN_PARALLEL,
                         ... )
  
  res$Outcome_desc = NULL
  res$Pretest = NULL
  
  # This is the same ID as built in the bottom of 'master_analyze_all_data.R'
  UniqueID     = str_remove_all(paste(proj, outcome, trt, sep = '_'), '_NA')
  res$UniqueID = UniqueID
  
  write_rds( res, 
             here::here( glue::glue( "{RESULT_DIR}/{UniqueID}.rds" ) ) )
  res
  
}



# OPTIONAL: Set up parallel code

if ( RUN_PARALLEL ) {
  library( furrr )
  oplan = plan(multisession, workers = parallel::detectCores() - 1 )
}



# Run the simulations ----



## First fake data ----

# using W score for WJ-AP outcome - not age-normalized to mean 100/SD 15 used in paper, but best we have
outcomes_mpc <- c('Yobs')

# no pre-tests for math outcomes (see footnote 1, pg. 85)
# other demographic covars not in data 
covars_mpc <- c('X', 'W1', 'W2', 'W3')

# grid for each outcome/covar combo
grid_mpc <- expand_grid(outcome = outcomes_mpc,
                        covars  = COVAR_SETS)
grid_mpc

# compile results
compare_mpc <- pmap(grid_mpc, run_sim, 
                    proj        = 'FAKE',
                    covars_list = covars_mpc,
                    data        = fake_data,
                    trt         = 'Z',
                    clust       = 'clusterID',
                    blk         = 'blockID', 
                    ptest       = 0,
                    .progress   = TRUE) %>%
  reduce(bind_rows)




# Close parallel processing ----

if ( RUN_PARALLEL ) {
  plan(oplan)
}



# Combine and save results ----

compare <- grep('compare_', ls(), value = T) %>%
  map_dfr(~ bind_rows(get(.x))) %>%
  relocate(Project, Outcome, Trt_arm, Covariates) %>%
  # filter to only exemplar CARES outcomes
  mutate(cares_keep = case_when(Project != 'CARES' ~ T,
                                Trt_arm == 'Incredible Years'  & Outcome %in% c('F1_DC_PT', 'F1_TR_BPI', 'F1_TR_WRS')            ~ T,
                                Trt_arm == 'Preschool PATHS'   & Outcome %in% c('F1_DC_EMIT', 'F1_DC_CST_COMPTOT', 'F1_TR_SSRS') ~ T,
                                Trt_arm == 'Tools of the Mind' & Outcome %in% c('F1_DC_PT', 'F1_TR_WRS')                         ~ T,
                                TRUE ~ F)) %>%
  filter(cares_keep == T) %>%
  select(-cares_keep) %>%
  rename(Estimator = method) %>%
  # map estimators to estimands
  mutate(Estimand = case_when(Estimator %in% c('LRi-FE-crve',    'LRi-FE-db',
                                               'LRi-FIpw-crve',  'LRi-FIpw-db',
                                               'LRapw-FE-het',   'LRapw-FE-db',
                                               'LRapw-FIpw-het', 'LRapw-FIpw-db')    ~ 'Person',
                              Estimator %in% c('LRicw-FE-crve',   'LRicw-FE-db',
                                               'LRicw-FIcw-crve', 'LRicw-FIcw-db',
                                               'LRa-FE-het',      'LRa-FE-db',
                                               'LRa-FIcw-het',    'LRa-FIcw-db',
                                               'MLM-FE', 'MLM-RE', 'MLM-FIcw')       ~ 'Cluster'))


as_tibble( compare )

saveRDS(compare,  
        file.path( RESULT_DIR, "master_fake_results.RDS" ) )



# Look at the simulation results ----

compare = readRDS( file.path( RESULT_DIR, "master_fake_results.RDS" ) )

compare <- as_tibble(compare)
compare

table( compare$runID )




