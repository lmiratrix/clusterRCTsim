

getwd()

library( tidyverse )

dat = read_csv( "empirical/jg_results.csv" ) 

head( dat ) 

stud = filter( dat, method == "LR_FE_CRVE", Covariates == "all" ) %>%
  dplyr::select( -c( method:df ), -Pretest, -Covariates, -Trt_arm, -Outcome )

stud %>%
  dplyr::select( Project, n, J, K, p.tx, ncov.1, ncov.2 )

