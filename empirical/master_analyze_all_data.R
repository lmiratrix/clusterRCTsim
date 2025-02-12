
# Core Analysis Script
# 
# This script performs the main empirical analysis by running 
# `compare_methods()` and `describe_clusterRCT()` on all non-clean 
# room datasets. The results are then combined with the clean room 
# output, stacked, and saved as `empirical/results/master241204_combined.RDS`.

#
# Files saved:
# - Three excel sheets
# - master240528_combined.RDS (or whatever date)

# was getting error running devtools::install_github(), so used this approach instead
# download.file("https://github.com/lmiratrix/clusterRCT/archive/refs/heads/main.zip", destfile = "/data/share/pse/IES CRT/clusterRCTsim/empirical/clusterRCT.zip")
# unzip("/data/share/pse/IES CRT/clusterRCTsim/empirical/clusterRCT.zip")
# install.packages("/data/share/pse/IES CRT/clusterRCTsim/clusterRCT-main", repos = NULL, type = "source")

# load libraries ----
# devtools::install_github("lmiratrix/clusterRCT")
library(tidyverse)
library(haven)
library(clusterRCT)
library(tools)
library(writexl)



# load cleaned data ----

# these are individual datasets for each non-clean room study, already filtered to analytic sample
for (f in list.files('/data/share/pse/IES CRT/sd7/clean study datasets', full.names = T)) {
  assign(file_path_sans_ext(basename(f)),
         readRDS(f))
}

# Helper functions for compiling results ----
# the functions below run compare_methods() and describe_clusterRCT() for a given outcome and set of covariates.
# writing the functions in this generalized way allows the code below to be streamlined and avoid repetition

# compare_methods() helper function
cm <- function(proj, outcome, covars, covars_list, data, trt, clust, blk, ptest) {
  
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
  
  # run compare_methods() three times: no covars, "best" covar, and all covars
  if (covars == 'none') {
    results <- data %>%
      compare_methods(form,
                      data = .,
                      include_method_characteristics = F)
  } else if (covars == 'all') {
    cf <- as.formula(paste('~', paste(covars_list, collapse = ' + ')))
    
    results <- data %>%
      compare_methods(form,
                      data = .,
                      control_formula = cf,
                      include_method_characteristics = F)
  } else if (covars == 'best') {
    # Calculate correlations
    correlations <- sapply(covars_list, function(covariate) {
      cor(data[[outcome]], data[[covariate]], use = "complete.obs")
    })
    
    # Find the covariate with the highest correlation
    best_covariate <- names(which.max(abs(correlations)))
    
    cf <- as.formula(paste('~', paste(best_covariate, collapse = ' + ')))
    
    results <- data %>%
      compare_methods(form,
                      data = .,
                      control_formula = cf,
                      include_method_characteristics = F)
  }
  
  # spit out results with project and outcome information
  results <- results %>%
    mutate(Project      = proj,
           Pretest      = ptest,
           Covariates   = covars,
           Trt_arm      = ifelse(proj == 'EMERGE', str_remove(trt, 'T_'), NA),
           Outcome      = outcome,
           Outcome_desc = attr(data[[outcome]], 'label'))
  
  return(results)
}

# describe_clusterRCT helper function
dc <- function(proj, outcome, covars, covars_list, data, trt, clust, blk, ptest) {
  
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
  
  # run describe_clusterRCT() three times: no covars, "best" covar, and all covars
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
  
  # spit out results with project and outcome info
  results <- results %>%
    mutate(Project      = proj,
           Pretest      = ptest,
           Covariates   = covars,
           Trt_arm      = ifelse(proj == 'EMERGE', str_remove(trt, 'T_'), NA),
           Outcome      = outcome,
           Outcome_desc = attr(data[[outcome]], 'label'))
  
  return(results)
}


# Run compare_methods() ----

## Making Pre-K Count ----

# using W score for WJ-AP outcome - not age-normalized to mean 100/SD 15 used in paper, but best we have
outcomes_mpc <- c('eclsb2_scale', 'wj_ap2')

# no pre-tests for math outcomes (see footnote 1, pg. 85)
# other demographic covars not in data 
covars_mpc <- c('age1')

# grid for each outcome/covar combo
grid_mpc <- expand_grid(outcome = outcomes_mpc,
                        covars  = c('none', 'best', 'all'))

# compile results
compare_mpc <- pmap(grid_mpc, cm, 
                    proj        = 'MPC',
                    covars_list = covars_mpc,
                    data        = mpc,
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
                           covars  = c('none', 'best', 'all'))

compare_emerge <- pmap(grid_emerge, cm,
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
                          covars  = c('none', 'best', 'all'))


compare_triad <- pmap(grid_triad, cm,
                      proj        = 'TRIAD',
                      covars_list = covars_triad,
                      data        = triad,
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
                            covars  = c('none', 'best', 'all'))


compare_triadsu <- pmap(grid_triadsu, cm, 
                        proj        = 'TRIAD Scale-up',
                        covars_list = covars_triadsu,
                        data        = triadsu,
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
                        covars  = c('none', 'best', 'all'))


compare_sfa <- pmap(grid_sfa, cm, 
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
                           covars  = c('none', 'best', 'all'))


compare_tnstar <- pmap(grid_tnstar, cm, 
                       proj        = 'Tennessee STAR',
                       covars_list = covars_tnstar,
                       data        = tnstar,
                       trt         = 'adj_racode',
                       clust       = 'K_TEACHER_ID',
                       blk         = 'K_SCHOOL_ID',
                       ptest       = 0,
                       .progress   = T) %>%
  reduce(bind_rows)

## Beacon Mentor ----

outcomes_beaconmentor <- c('S8TGPA', 'S8TOTCR', 'SRPRST2')
covars_beaconmentor <- c('RADEV', 'S8TEACHID', 'S8FULTMEC', 'ZFemale', 'RATIME', 'ZAGELT24', 'PTTOP75C', 'PTBTM25C')

grid_beaconmentor <- expand_grid(outcome = outcomes_beaconmentor,
                                 covars  = c('none', 'best', 'all'))


compare_beaconmentor <- pmap(grid_beaconmentor, cm, 
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
                                  covars  = c('none', 'best', 'all'))


compare_powerteaching <- pmap(grid_powerteaching, cm, 
                              proj        = 'PowerTeaching',
                              covars_list = covars_powerteaching,
                              data        = powerteaching,
                              trt         = 'RA_CODE',
                              clust       = 'SCHOOL_ID',
                              blk         = 'DIST',
                              ptest       = 1,
                              .progress   = T) %>%
  reduce(bind_rows)

# Run describe_clusterRCT() ----

## Making Pre-K Count ----

describe_mpc <- pmap(grid_mpc, dc, 
                    proj        = 'MPC',
                    covars_list = covars_mpc,
                    data        = mpc,
                    trt         = 'dtx',
                    clust       = 'school1',
                    blk         = 'block', 
                    ptest       = 0,
                    .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

## EMERGE ----

describe_emerge <- pmap(grid_emerge, dc,
                        proj        = 'EMERGE',
                        covars_list = covars_emerge,
                        data        = emerge,
                        clust       = 'teachid',
                        blk         = 'block',
                        ptest       = 0,
                        .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

## TRIAD ----

describe_triad <- pmap(grid_triad,dc,
                      proj        = 'TRIAD',
                      covars_list = covars_triad,
                      data        = triad,
                      trt         = 'dtx',
                      clust       = 'school1',
                      blk         = 'block',
                      ptest       = 1,
                      .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

## TRIAD Scale-up ---- 

describe_triadsu <- pmap(grid_triadsu, dc, 
                        proj        = 'TRIAD Scale-up',
                        covars_list = covars_triadsu,
                        data        = triadsu,
                        trt         = 'dtx',
                        clust       = 'school1',
                        blk         = 'block',
                        ptest       = 1,
                        .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

## Success for All ----

describe_sfa <- pmap(grid_sfa, dc, 
                    proj        = 'Success for All',
                    covars_list = covars_sfa,
                    data        = sfa,
                    trt         = 'RA_CODE',
                    clust       = 'SCHOOL_ID',
                    blk         = 'DISTRICT',
                    ptest       = 1,
                    .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()


## TN Star ----

describe_tnstar <- pmap(grid_tnstar, dc, 
                       proj        = 'Tennessee STAR',
                       covars_list = covars_tnstar,
                       data        = tnstar,
                       trt         = 'adj_racode',
                       clust       = 'K_TEACHER_ID',
                       blk         = 'K_SCHOOL_ID',
                       ptest       = 0,
                       .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

## Beacon Mentor ----

describe_beaconmentor <- pmap(grid_beaconmentor, dc, 
                             proj        = 'Beacon Mentor',
                             covars_list = covars_beaconmentor,
                             data        = beaconmentor,
                             trt         = 'racode_b',
                             clust       = 'S8CLUSTE',
                             blk         = 'RABLOCK',
                             ptest       = 0,
                             .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

## PowerTeaching ----

describe_powerteaching <- pmap(grid_powerteaching, dc, 
                              proj        = 'PowerTeaching',
                              covars_list = covars_powerteaching,
                              data        = powerteaching,
                              trt         = 'RA_CODE',
                              clust       = 'SCHOOL_ID',
                              blk         = 'DIST',
                              ptest       = 1,
                              .progress   = T) %>%
  reduce(bind_rows) %>%
  as.data.frame()

# Compile and save results ----

## Read in output from clean room datasets ----
# the same helper functions were run on the three clean room studies, and the results were compiled

compare_ies <- readRDS("/data/share/pse/IES CRT/sd7/iesdata_comparemethods_v2.RDS")
describe_ies <- readRDS("/data/share/pse/IES CRT/sd7/iesdata_describeclusterrct_v2.RDS")

## Combine and save compare_methods() results ----

compare <- grep('compare_', ls(), value = T) %>%
  map_dfr(~ bind_rows(get(.x))) %>%
  relocate(Project, Outcome, Outcome_desc, Trt_arm, Pretest, Covariates) %>%
  # filter to only exemplar CARES outcomes
  mutate(cares_keep = case_when(Project != 'CARES' ~ T,
                                Trt_arm == 'Incredible Years'  & Outcome %in% c('F1_DC_PT', 'F1_TR_BPI', 'F1_TR_WRS')            ~ T,
                                Trt_arm == 'Preschool PATHS'   & Outcome %in% c('F1_DC_EMIT', 'F1_DC_CST_COMPTOT', 'F1_TR_SSRS') ~ T,
                                Trt_arm == 'Tools of the Mind' & Outcome %in% c('F1_DC_PT', 'F1_TR_WRS')                         ~ T,
                                TRUE ~ F),
         Outcome_desc = str_trim(Outcome_desc)) %>%
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
                                               'MLM-FE', 'MLM-RE', 'MLM-FIcw')       ~ 'Cluster'),
         
         Outcome_desc = str_trim(Outcome_desc))

## Combine and save describe_clusterRCT() results ----

describe <- grep('describe_', ls(), value = T) %>%
  map_dfr(~ bind_rows(get(.x))) %>%
  relocate(Project, Outcome, Outcome_desc, Trt_arm, Pretest) %>%
  relocate(SER2.1, .after = R2.1) %>%
  relocate(SER2.2, .after = R2.2) %>%
  # filter to only exemplar CARES outcomes
  mutate(cares_keep = case_when(Project != 'CARES' ~ T,
                                Trt_arm == 'Incredible Years'  & Outcome %in% c('F1_DC_PT', 'F1_TR_BPI', 'F1_TR_WRS')            ~ T,
                                Trt_arm == 'Preschool PATHS'   & Outcome %in% c('F1_DC_EMIT', 'F1_DC_CST_COMPTOT', 'F1_TR_SSRS') ~ T,
                                Trt_arm == 'Tools of the Mind' & Outcome %in% c('F1_DC_PT', 'F1_TR_WRS')                         ~ T,
                                TRUE ~ F),
         Outcome_desc = str_trim(Outcome_desc)) %>%
  filter(cares_keep == T) %>%
  select(-cares_keep)

## Create and save merged master dataset ----

combined <- inner_join(describe %>%
                         select(-Pretest),
                       compare) %>%
  # create a project-outcome-arm level identifier and standardized (ES units) ATE and SE
  mutate(UniqueID     = str_remove_all(paste(Project, Outcome, Trt_arm, sep = '_'), '_NA'),
         ATE_hat_std  = ATE_hat / sdY0,
         SE_hat_std   = SE_hat / sdY0) %>%
  relocate(UniqueID) %>%
  relocate(ATE_hat_std, SE_hat_std, .after = SE_hat)

write_xlsx(compare,  '/data/share/pse/IES CRT/sd7/master241204_compare.xlsx')
write_xlsx(describe, '/data/share/pse/IES CRT/sd7/master241204_describe.xlsx') 
write_xlsx(combined, '/data/share/pse/IES CRT/sd7/master241204_combined.xlsx') 
saveRDS(combined,    '/data/share/pse/IES CRT/sd7/master241204_combined.RDS')
saveRDS(combined,    '/data/share/pse/IES CRT/clusterRCTsim/empirical/results/master241204_combined.RDS')

