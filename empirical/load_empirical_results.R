
# Load libraries and set of results
#
# Also set paths for storing generated results.


# Load libraries and data ----

library(tidyverse)
library(readxl)
library(gridExtra)
library(openxlsx)

options(list(dplyr.summarise.inform = FALSE))
theme_set( theme_classic() )



PLOT_PATH = '/data/share/pse/IES CRT/tables/'
PLOT_PATH = paste0( here::here( "empirical", "tables" ), "/" )

plot_path <- function( filename ) {
  here::here( PLOT_PATH, filename )
}
if ( !dir.exists( PLOT_PATH ) ) {
  cli::cli_alert("Making directory to store tables and figures at {PLOT_PATH}" )
  dir.create( PLOT_PATH, recursive = T )
}



# Load the results from compare() ----

combined <-  read_rds( here::here( "empirical/results/master241204_combined.RDS" ) ) %>%
  as_tibble()

# add in literature info for curated studies 
curated_litinfo <- read_excel(here::here("empirical/results/study_research_info.xlsx")) %>%
  filter(Collection == 'Curated')

combined <- left_join(combined, curated_litinfo %>%
                        select(Project, Outcome, Trt_arm, lead_author:article_url))

combined



# Add in the Josh datasets ----

if ( TRUE ) {
  # Add in josh results
  josh = read_csv( here::here( "empirical", "results", "jg_results.csv" ) )
  
  #josh <- josh %>%
  #  mutate(Name_Year = str_c(str_extract(lead_author, "^[^,]+"), "-", year))
  
  
  names(josh)
  head( josh )
  josh$UniqueID = josh$Project
  
  josh$Project = paste0( str_extract(josh$lead_author, "^[^,]+"), "/set_", josh$study_id )
  
  
  josh %>%
    dplyr::filter( blocked == 1 ) %>%
    dplyr::select( Project, Outcome, sample, dataset_id, study_id ) %>%
    unique() %>%
    print( n = 25 )
  
  
  tb = sort( table( josh$study_id ) )
  tb2 = sort( table( josh$Project ) )
  stopifnot( all( tb == tb2 ) )
  rm( tb, tb2 )
  
  table( josh$UniqueID )
  n_distinct(josh$UniqueID)
  skimr::skim( josh )
  
  setdiff( names(combined), names(josh) )
  setdiff( names(josh), names(combined) )
  
  josh <- josh %>%
    rename( Estimator = method,
            Estimand = weight,
            ATE_hat_std = ATE_hat,
            SE_hat_std = SE_hat )
  
  #%>%
  # mutate( UniqueID = paste0( UniqueID, "/", id ) ) %>%
  #  dplyr::select( -id )
  
  length( unique( josh$UniqueID ) )
  
  # All original datasets are blocked
  combined$blocked = 1
  
  combined = bind_rows( `Curated` = combined, `Public` = josh, .id="Collection" )
  names(combined)
  rm( josh )
}



# Check all estimators have same estimand for all estimates
table( combined$Estimator, combined$Estimand )

source( here::here( "empirical/R/clean_result_table.R" ) )
combined <- combined %>%
  clean_result_table( ) %>%
  make_df_flags()


# Divide the datasets into blocked and non-blocked ----


unblocked = filter( combined, blocked == 0 ) %>%
  dplyr::select( -blocked )

combined = filter( combined, blocked == 1 ) %>%
  dplyr::select( -blocked )


# Drop CRVE FI SEs for singleton cases ----

table( combined$Estimator )

drp <- combined$Estimator %in% c( "LRicw-FIcw-crve",
                                  "LRi-FIpw-crve" ) &
  combined$num_singletons > 0
table( drp )

combined$SE_hat
combined$SE_hat[ drp ] <- NA

combined$SE_hat_std[ drp ] <- NA

if ( FALSE ) {
filter( combined, UniqueID == "jg_crt_6", Covariates == "best") %>%
  dplyr::select( Estimator, num_singletons, ATE_hat_std, SE_hat, SE_hat_std )
}

rm( drp )





# Make table of all the outcomes ----
outs <- combined %>% 
  filter( Covariates == "best" ) %>%
  filter( !duplicated( UniqueID ) ) %>%
  dplyr::select( -c( Estimator:df ), -Covariates, -c( Estimand:disfavored ) )

outs

outs %>%
  group_by( Collection ) %>%
  summarise( n_study = n_distinct(Project),
             n_outcome = n_distinct(UniqueID) )


# Make studies table ----


setdiff( names( combined ),
         names( unblocked ) )
setdiff(   names( unblocked ),
           names( combined ))

studies <- bind_rows( Blocked = combined, Unblocked = unblocked, .id="blocked" ) %>%
  dplyr::filter( Covariates == "best" )

n_defined <- studies %>% 
  group_by( UniqueID ) %>%
  summarise( n_ATE = sum( !is.na( ATE_hat ) ),
             n_SE = sum( !is.na( SE_hat ) ) )

studies <- studies %>%
  dplyr::filter( Estimator %in% c( "LRi-crve", "LRi-FE-crve" ) ) %>%
  dplyr::select( -c( Estimator:disfavored ) )

studies <- left_join( studies, n_defined, by="UniqueID" )

rm( n_defined )

cat( "Data loaded\n" )
