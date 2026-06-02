
# This is a version of master_analyze_all_data.R except it loads the
# data and runs a fingerprint simulation on each dataset.


# load libraries ----
# devtools::install_github("lmiratrix/clusterRCT")
library(tidyverse)
library(haven)
library(clusterRCT)
library(tools)
library(writexl)


source( here::here( "empirical/data_processing_functions.R" ) )

# Where the prepped data is
DATA_DIR = '/data/share/pse/IES CRT/sd7/clean study datasets'


# Where to save result files (one file per outcome will be generated)
RESULT_DIR = here::here( "rerand_sim_results" )


# How many iterations for each outcome?
NUM_SIMULATIONS = 20

# What configuration to we run our simulation over?
COVAR_SETS = c('none', 'best', 'all')

# Run in parallel?
RUN_PARALLEL = TRUE

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



## Fake CRT (demo) ----

if ( FALSE ) {
  data( fakeCRT )
  fakeCRT$Xg = rnorm( nrow(fakeCRT), mean=0, sd=1 )
  fakeCRT$Y2 = fakeCRT$Xg + rnorm( nrow(fakeCRT), sd=1 )
  fakeCRT$XX = rnorm( nrow(fakeCRT) )
  
  outcomes_fake <- c('Yobs', 'Y2')
  
  covars_fake <- c('X.jk', 'C.ijk', "XX", "Xg" )
  
  
  # grid for each outcome/covar combo
  grid_fake <- expand_grid(outcome = outcomes_fake,
                           covars  = COVAR_SETS)
  
  # Run sim across all scenarios
  compare_fake <- pmap(grid_fake, run_sim, 
                       proj        = 'fakeCRT',
                       covars_list = covars_fake,
                       data        = fakeCRT,
                       trt         = 'T.x',
                       clust       = 'S.id',
                       blk         = 'D.id', 
                       ptest       = 0,
                       # R = 2, include_MLM = FALSE,
                       .progress   = TRUE) %>%
    reduce(bind_rows)

  table( compare_fake$runID, compare_fake$Outcome )
}


## Making Pre-K Count ----

# using W score for WJ-AP outcome - not age-normalized to mean 100/SD 15 used in paper, but best we have
outcomes_mpc <- c('eclsb2_scale', 'wj_ap2')

# no pre-tests for math outcomes (see footnote 1, pg. 85)
# other demographic covars not in data 
covars_mpc <- c('age1')

# grid for each outcome/covar combo
grid_mpc <- expand_grid(outcome = outcomes_mpc,
                        covars  = COVAR_SETS)

# compile results
compare_mpc <- pmap(grid_mpc, run_sim, 
                    proj        = 'MPC',
                    covars_list = covars_mpc,
                    data        = zap_labels(mpc),
                    trt         = 'dtx',
                    clust       = 'school1',
                    blk         = 'block', 
                    ptest       = 0,
                    .progress   = T) %>%
  reduce(bind_rows)


## EMERGE ----

outcomes_emerge <- c('team2_scaled', 'eclsb2')

covars_emerge <- c('age1', 'dfemale')

# EMERGE grid also includes treatment arm
grid_emerge <- expand_grid(trt = c('T_BB', 'T_BBSSR'),
                           outcome = outcomes_emerge,
                           covars  = COVAR_SETS)

compare_emerge <- pmap(grid_emerge, run_sim,
                       proj        = 'EMERGE',
                       covars_list = covars_emerge,
                       data        = emerge,
                       clust       = 'teachid',
                       blk         = 'block',
                       ptest       = 0,
                       .progress   = T) %>%
  reduce(bind_rows)

## TRIAD ----

# REMA seems to be on a very different scale than report (see table 3 in report)
# They were converted to a "T score", which is explained on pg. 71 of the 2018 report:
# https://www.mdrc.org/sites/default/files/MPC-High_5s_Impact_FR_0.pdf
# I attempted to recreate the T score and am close to the averages in Table 3, but not quite exact.
outcomes_triad <- c('rema2')

# just pretest and af-am based on table 4 in report
covars_triad <- c('rema1', 'dblack')

grid_triad <- expand_grid(outcome = outcomes_triad,
                          covars  = COVAR_SETS)


compare_triad <- pmap(grid_triad, run_sim,
                      proj        = 'TRIAD',
                      covars_list = covars_triad,
                      data        = zap_labels(triad),
                      trt         = 'dtx',
                      clust       = 'school1',
                      blk         = 'block',
                      ptest       = 1,
                      .progress   = T) %>%
  reduce(bind_rows)

## TRIAD Scale-up ----


# REMA seems to be on a very different scale than report (see table 2 in report)
# maybe T-scores?
outcomes_triadsu <- c('rema2')

# just pretest and af-am based on table 4 in report
covars_triadsu <- c('rema1', 'dblack')

grid_triadsu <- expand_grid(outcome = outcomes_triadsu,
                            covars  = COVAR_SETS)


compare_triadsu <- pmap(grid_triadsu, run_sim, 
                        proj        = 'TRIAD Scale-up',
                        covars_list = covars_triadsu,
                        data        = zap_labels(triadsu),
                        trt         = 'dtx',
                        clust       = 'school1',
                        blk         = 'block',
                        ptest       = 1,
                        .progress   = T) %>%
  reduce(bind_rows)

## Success for All ----

outcomes_sfa <- 'S14WJLWI_SCORE'
covars_sfa <- c('SPECIAL_ED', 'MALE', 'ELL', 'AGE', 'F11PPVT_SCORE', 'F11WJLWI_SCORE')

grid_sfa <- expand_grid(outcome = outcomes_sfa,
                        covars  = COVAR_SETS)


compare_sfa <- pmap(grid_sfa, run_sim, 
                    proj        = 'Success for All',
                    covars_list = covars_sfa,
                    data        = sfa,
                    trt         = 'RA_CODE',
                    clust       = 'SCHOOL_ID',
                    blk         = 'DISTRICT',
                    ptest       = 1,
                    .progress   = T) %>%
  reduce(bind_rows)

## TN Star ----

outcomes_tnstar <- c('gktmathss', 'gktreadss')
covars_tnstar <- c('gender_na', 'female', 'male', 
                   'race_ai', 'race_as_pi', 'race_o', 'race_hi', 
                   'race_na', 'race_b', 'race_w', 
                   'gkfreelunch', 'K_SPECIAL_ED')

grid_tnstar <- expand_grid(outcome = outcomes_tnstar,
                           covars  = COVAR_SETS)


compare_tnstar <- pmap(grid_tnstar, run_sim, 
                       proj        = 'Tennessee STAR',
                       covars_list = covars_tnstar,
                       data        = zap_labels(tnstar),
                       trt         = 'adj_racode',
                       clust       = 'K_TEACHER_ID',
                       blk         = 'K_SCHOOL_ID',
                       ptest       = 0,
                       .progress   = T) %>%
  reduce(bind_rows)

# testing on smaller dataset, one grid row at a time
# issue might be large number of clusters and blocks - try with fewer

# tnstar_test <- tnstar %>%
#   filter(K_SCHOOL_ID %in% unique(.$K_SCHOOL_ID)[1:5])
# 
# compare_tnstar_test <- pmap(grid_tnstar[5,], run_sim, 
#                        proj        = 'Tennessee STAR',
#                        covars_list = covars_tnstar,
#                        data        = zap_labels(tnstar_test),
#                        trt         = 'adj_racode',
#                        clust       = 'K_TEACHER_ID',
#                        blk         = 'K_SCHOOL_ID',
#                        ptest       = 0,
#                        .progress   = T) %>%
#   reduce(bind_rows)

## Beacon Mentor ----

outcomes_beaconmentor <- c('S8TGPA', 'S8TOTCR', 'SRPRST2')
covars_beaconmentor <- c('RADEV', 'S8TEACHID', 'S8FULTMEC', 'ZFemale', 'RATIME', 'ZAGELT24', 'PTTOP75C', 'PTBTM25C')

grid_beaconmentor <- expand_grid(outcome = outcomes_beaconmentor,
                                 covars  = COVAR_SETS)


compare_beaconmentor <- pmap(grid_beaconmentor, run_sim, 
                             proj        = 'Beacon Mentor',
                             covars_list = covars_beaconmentor,
                             data        = beaconmentor,
                             trt         = 'racode_b',
                             clust       = 'S8CLUSTE',
                             blk         = 'RABLOCK',
                             ptest       = 0,
                             .progress   = T) %>%
  reduce(bind_rows)


## PowerTeaching ----

outcomes_powerteaching <- c('TSTSSTAND')
covars_powerteaching <- c('COHORT1', 'COHORT2', 'MALE_BL',
                          'HISPANIC_BL', 'WHITE_BL', 'BLACK_BL',
                          'FRLUNCH_BL', 'ELL_BL', 'SPECIAL_ED_BL', 
                          'TSTSSTAND_BL', 'TSTSS_BL', 'TSTPLVL_BL')

grid_powerteaching <- expand_grid(outcome = outcomes_powerteaching,
                                  covars  = COVAR_SETS)


compare_powerteaching <- pmap(grid_powerteaching, run_sim, 
                              proj        = 'PowerTeaching',
                              covars_list = covars_powerteaching,
                              data        = powerteaching,
                              trt         = 'RA_CODE',
                              clust       = 'SCHOOL_ID',
                              blk         = 'DIST',
                              ptest       = 1,
                              .progress   = T) %>%
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



write_xlsx(compare,  '/data/share/pse/IES CRT/sd7/master_run_fingerprints.xlsx')
saveRDS(compare,    '/data/share/pse/IES CRT/sd7/master_run_fingerprints.RDS')


