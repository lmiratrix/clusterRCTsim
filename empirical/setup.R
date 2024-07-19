
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

combined <-  read_rds( here::here( "empirical/results/master240528_combined.RDS" ) ) %>%
  as_tibble()

combined

if ( TRUE ) {
  # Add in josh results
  josh = read_csv( here::here( "empirical", "results", "jg_results.csv" ) )

  head( josh )
  josh$UniqueID = josh$Project
  
  table( josh$UniqueID )
  skimr::skim( josh )
  
  setdiff( names(combined), names(josh) )
  setdiff( names(josh), names(combined) )
  
  josh <- josh %>%
    rename( Estimator = method,
            ATE_hat_std = ATE_hat,
            SE_hat_std = SE_hat ) %>%
    mutate( UniqueID = paste0( UniqueID, "/", id ) ) %>%
    dplyr::select( -id )
  
  length( unique( josh$UniqueID ) )
  
  # All original datasets are blocked
  combined$blocked = 1

  josh_blocked = josh   %>%
    filter( blocked == 1 )
  
  josh_unblocked = josh %>%
    filter( blocked == 0 )
  josh_unblocked$Collection ="JG"
  length( unique( josh_unblocked$UniqueID ) )
  
  combined = bind_rows( `Curated` = combined, `JG` = josh_blocked, .id="Collection" )
  names(combined)
  rm( josh_blocked )
  rm( josh )

}

#head( combined$UniqueID)
#table( combined$Trt_arm)


# Remove "of the" from all the names, etc
combined$UniqueID <- gsub("Preschool", "", combined$UniqueID)
combined$UniqueID <- gsub("Incredible Years", "Incred", combined$UniqueID)
combined$UniqueID <- gsub("_F1", "", combined$UniqueID)
combined$UniqueID <- gsub("of the", "", combined$UniqueID)
combined$UniqueID <- gsub("PATHS", "PTH", combined$UniqueID)
combined$UniqueID <- gsub(" Mentor", "", combined$UniqueID)
combined$UniqueID <- gsub("Success for All", "Success", combined$UniqueID)
combined$UniqueID <- gsub("Tennessee STAR", "TennStar", combined$UniqueID)
combined$UniqueID <- gsub("scaled", "", combined$UniqueID)
combined$UniqueID <- gsub("  ", "_", combined$UniqueID)
combined$UniqueID <- gsub(" ", "", combined$UniqueID)

table( combined$UniqueID, useNA = "always" )
