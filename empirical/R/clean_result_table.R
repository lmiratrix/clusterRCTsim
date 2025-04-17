
# Function to clean the combined results and make some useful descriptive covariates



clean_result_table <- function(combined ) {
  
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
  
  
  if ( FALSE ) {
    cc <- combined %>%
      dplyr::select( UniqueID, Estimator, SE_hat, n, J, K, df, Estimand, Covariates ) %>%
      arrange( UniqueID, Estimator ) %>%
      filter( Estimator == "DB_FI_Cluster_Cluster" ) %>%
      print( n = 20 )
    
    cc %>%
      dplyr::select( -SE_hat ) %>%
      group_by( UniqueID ) %>% 
      pivot_wider( names_from = Covariates, values_from = c( df ) ) %>%
      mutate( del = none - best ) %>%
      dplyr::filter( del != 1 ) %>%
      print( n = 30 )
  }
  
  combined
}



make_df_flags <- function( combined ) {
  
  # Calculate some data set descriptive such as low DF datasets ----
  
  combined <- combined %>% 
    group_by( UniqueID, Covariates ) %>%
    mutate( low_df = ifelse( min( df ) < 10, "low", "high" ),
            cons_df = J - 2*K,
            has_single = num_singletons > 0,
            class = case_when(
              low_df == "low" & has_single ~ "Low df & singletons",
              low_df == "low" & !has_single ~ "Low df only",
              low_df != "low" & has_single ~ "Singletons",
              TRUE ~ "Large"
            ),
            has_single = ifelse( has_single, "singetons", "no singletons" ) ) %>%
    ungroup()
  
  combined
}

