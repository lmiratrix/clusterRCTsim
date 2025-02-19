

# Make table of features for the different estimators



# Make list of what estimators are of what type ----

if ( FALSE ) {
  
  
  source( here::here( "empirical/load_empirical_results.R" ) )
  source( here::here( "empirical/makeDotplot.R" ) )
  
  
  ee = unique( combined$Estimator )
  ee
  tt = tibble( Estimator = ee,
               type = case_when( grepl( "MLM", ee ) ~ "MLM",
                                 grepl( "LRa", ee ) ~ "LRa",
                                 grepl( "LRi", ee ) ~ "LRi" ),
               FI = 0+ grepl( "-FI", ee ),
               LR = 0+ grepl( "LR", ee ),
               agg = 0+ grepl( "LRa", ee ),
               db = 0+ grepl( "db", ee ),
               crve = 0+ grepl( "crve", ee ),
               het = 0+ grepl( "het", ee ) )
  
  
  table( combined$UniqueID )
  
  rm( tt )
  rm( ee )
  
  est_profile <- est_profile %>%
    mutate( est_class = paste0( ifelse( LR, "LR", "MLM" ),
                                ifelse( FI, "-FI", "-FE" ) ) )
  
  
  st <- filter( combined, 
                Project == "PowerTeaching",
                Covariates == "best" ) %>%
    dplyr::select( Estimator, ATE_hat, SE_hat ) %>%
    mutate( dup = 0 + duplicated( round( ATE_hat, digits = 7 ) ) )
  
  est_profile <- left_join( est_profile, 
                            dplyr::select( st, Estimator, dup ),
                            by ="Estimator")
  
  
  rm( st )
  
  est_profile
  datapasta::tribble_paste( est_profile )

}

est_profile <- tibble::tribble(
          ~Estimator, ~type, ~FI, ~LR, ~agg, ~db, ~crve, ~het, ~est_class, ~dup,
         "LRa-FE-db", "LRa",   0,   1,    1,   1,     0,    0,    "LR-FE",    0,
        "LRa-FE-het", "LRa",   0,   1,    1,   0,     0,    1,    "LR-FE",    1,
       "LRa-FIcw-db", "LRa",   1,   1,    1,   1,     0,    0,    "LR-FI",    0,
      "LRa-FIcw-het", "LRa",   1,   1,    1,   0,     0,    1,    "LR-FI",    1,
       "LRapw-FE-db", "LRa",   0,   1,    1,   1,     0,    0,    "LR-FE",    0,
      "LRapw-FE-het", "LRa",   0,   1,    1,   0,     0,    1,    "LR-FE",    1,
     "LRapw-FIpw-db", "LRa",   1,   1,    1,   1,     0,    0,    "LR-FI",    0,
    "LRapw-FIpw-het", "LRa",   1,   1,    1,   0,     0,    1,    "LR-FI",    1,
       "LRi-FE-crve", "LRi",   0,   1,    0,   0,     1,    0,    "LR-FE",    0,
         "LRi-FE-db", "LRi",   0,   1,    0,   1,     0,    0,    "LR-FE",    1,
     "LRi-FIpw-crve", "LRi",   1,   1,    0,   0,     1,    0,    "LR-FI",    0,
       "LRi-FIpw-db", "LRi",   1,   1,    0,   1,     0,    0,    "LR-FI",    1,
     "LRicw-FE-crve", "LRi",   0,   1,    0,   0,     1,    0,    "LR-FE",    0,
       "LRicw-FE-db", "LRi",   0,   1,    0,   1,     0,    0,    "LR-FE",    1,
   "LRicw-FIcw-crve", "LRi",   1,   1,    0,   0,     1,    0,    "LR-FI",    0,
     "LRicw-FIcw-db", "LRi",   1,   1,    0,   1,     0,    0,    "LR-FI",    1,
            "MLM-FE", "MLM",   0,   0,    0,   0,     0,    0,   "MLM-FE",    0,
          "MLM-FIcw", "MLM",   1,   0,    0,   0,     0,    0,   "MLM-FI",    0,
            "MLM-RE", "MLM",   0,   0,    0,   0,     0,    0,   "MLM-FE",    0
   )

