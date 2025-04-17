


library( tidyverse )


res <- read_rds( "fingerprint/results_master_finger/master_run_fingerprints.RDS" ) %>%
  as_tibble()
res
table( res$runID, res$Outcome )

table( res$Outcome, res$Project )

rem = filter( res, Outcome == "rema2" )
table( rem$Project )
