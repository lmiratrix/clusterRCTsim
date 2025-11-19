
# Load and prep empirical results and calculate various summary measures and flags


#j_studies = read_csv( here::here( "empirical/results/jg_describe.csv" ) )
source( here::here( "empirical/load_empirical_results.R" ) )


source( here::here( "empirical/calculate_all_ranges.R" ) )

setdiff( summary_stats$UniqueID, summary_stats$UniqueID  )

summary_stats$Estimand = "Overall"
summary_stats_exclFI$Estimand = "Overall"
summary_stats_exclFI$Estimand = paste0( summary_stats_exclFI$Estimand, "_noFI" )
summary_stats_estimand_exclFI$Estimand = paste0( summary_stats_estimand_exclFI$Estimand, "_noFI" )

ranges = bind_rows( all = summary_stats, 
                    estimand = summary_stats_estimand,
                    noFI = summary_stats_exclFI,
                    noFI_estimand = summary_stats_estimand_exclFI,
                    .id = "group" )

table( ranges$group )


rangesW = ranges
rangesW  <- left_join( rangesW %>%
                         dplyr::select( group, Covariates, UniqueID, ATE_N:Estimand ),
                       studies, by=c( "UniqueID", "Covariates") )
names(rangesW)

# Massage range to big dataframe ----

ranges <- ranges %>%
  dplyr::select( group, UniqueID, Covariates, Estimand, Range, Ratio ) %>%
  pivot_longer( cols=c( Range, Ratio ),
                names_to = "measure", values_to = "Range" ) %>%
  mutate( measure = fct_recode( measure,
                                ate = "Range",
                                se = "Ratio" ) ) %>%
  filter( Covariates == "best" ) 


table( ranges$Estimand, ranges$measure )

filter( summary_stats,  UniqueID == "TRIAD_rema2" )
print( filter( ranges, UniqueID == "TRIAD_rema2" ), n = 30 )

rm( summary_stats, summary_stats_estimand, summary_stats_estimand_exclFI, summary_stats_exclFI )

ranges <- mutate( ranges,
                  measure = paste0( measure, "_", Estimand ),
                  measure = str_replace( measure, "_est", "" ),
                  measure = str_replace( measure, "_Overall", "" ) )
table( ranges$measure )

ranges <- mutate( ranges,
                  ate = str_detect(measure, "ate"),
                  big = ifelse( Range > ifelse( ate, 0.10, 2 ), "***", "" ) ) 

table( ranges$UniqueID )
table( ranges$measure )


ranges

# Save study descriptions to file ----

rW <- ranges %>%
  dplyr::select( UniqueID, measure, Range, big ) %>%
  pivot_wider( names_from = measure, values_from = c(Range, big), names_vary = "slowest" )
rW

ranges <- left_join( ranges, studies, by="UniqueID" )

stopifnot( all( ranges$Covariates.x == ranges$Covariates.y ) )
ranges <- ranges %>%
  rename( Covariates = Covariates.x ) %>%
  dplyr::select( -Covariates.y )

# We only have blocked studies
table( ranges$blocked )
setdiff( studies$UniqueID, ranges$UniqueID )
setdiff( ranges$UniqueID, studies$UniqueID )


studies <- left_join( studies, rW, by="UniqueID" )

names(studies)


write_csv( studies, here::here( "empirical/results/study_descriptions.csv" ) )
names( studies )
nrow( studies )

colnames(studies) <- str_replace( colnames(studies), "_Cluster", "_Clu" )
colnames(studies) <- str_replace( colnames(studies), "_Person", "_Per" )
colnames(studies) <- str_replace( colnames(studies), "_exclFI", "_noFI" )
colnames(studies) <- str_replace( colnames(studies), "Range_", "R_" )


