
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

cat( "Data loaded\n" )
