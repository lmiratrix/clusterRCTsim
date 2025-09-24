

# Make all the figures for the main empirical paper (for the blocked
# datasets)
#
# Colin, edited by Luke
# 2024


# load data ----

# load all data
source( here::here( "empirical/load_empirical_results.R" ) )
source( here::here( "empirical/makeDotplot.R" ) )


if ( FALSE ) {
  # filter to curated studies
  combined <- combined %>%
    filter(Collection == 'Curated')
}


source( here::here( "empirical/calculate_all_ranges.R" ) )



# Which studies are troublemakers? ----

summary_stats

summary_stats %>% 
  ungroup() %>%
  filter( Ratio > 2 ) %>%
  dplyr::select( UniqueID, Covariates, Collection, low_df, has_single, Ratio ) %>%
  arrange( Covariates, Ratio ) %>%
  knitr::kable( digits = 1 )

if ( FALSE ) {
  filter( combined, UniqueID == "CARES_DC_CST_COMPTOT_PTH", Covariates == "best") %>%
    dplyr::select( Estimator, num_singletons, ATE_hat_std, SE_hat_std ) %>%
    arrange( SE_hat_std )
  
  summary_stats %>%
    filter( UniqueID == "CARES_DC_CST_COMPTOT_PTH" ) %>%
    dplyr::select( Covariates, Range, Ratio, low_df, has_single )
}



# Make Paper Figures ----


## Figure 1: The Range (across 11 Estimators) of Estimates of β for {n_outcomes} Outcomes ----
table( summary_stats$class )
nrow( summary_stats )

range( summary_stats$Range )
n_outcomes <- n_distinct( summary_stats$UniqueID )
n_estimators <- n_distinct( combined$Estimator )
n_estimators
nrow( summary_stats )
table( table( summary_stats$UniqueID ) )

range_x_range <- 0.30

summary_stats %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β', 
              ylab     = 'Number of Outcomes', 
              title    = 'Figure 1', 
              filename   = 'final paper figures/fig1_ATErange',
              subtitle = glue::glue( 'The Range (across {n_estimators} Estimators) of Estimates of β for {n_outcomes} Outcomes' ),
              binw = 0.01,
              coordratio = 0.1,
              rmlegend = FALSE) +
  scale_x_continuous( limits=c(0, range_x_range ) )



## Figure 2: The Range of Estimates of β for all Estimators Targeting a Given Estimand  ----
summary_stats_estimand %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β',
              ylab     = 'Number of Outcomes',
              filename = 'final paper figures/fig2_ATErange_byestimand',
              title = 'Figure 2',
              subtitle = glue::glue( 'The Range of Estimates of β for all Estimators\nTargeting a Given Estimand for {n_outcomes} Outcomes' ),
              coordratio = .15,
              axisfs = 14,
              binw = 0.01,
              rmlegend = F,
              byest    = TRUE,
              xscale = c(0, 0.3)) 
# scale_x_continuous( limits=c(0, range_x_range ) )



## Figure 3: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) for {n_outcomes} outcomes. ----

summary_stats %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 3',
              subtitle = glue::glue( 'The Ratio (across {n_estimators} estimators) of the Largest to Smallest\nEstimates of SE(β) for {n_outcomes} outcomes' ),
              filename = 'final paper figures/fig3_SEratio',
              rmlegend = FALSE,
              binw     = 0.1,
              coordratio = 1.25,
              axisfs = 14,
              tcratio  = T,
              tcbreaks = c(1:4))



## Figure 4: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) among all estimators targeting a given estimand for {n_outcomes} outcomes ----

summary_stats_estimand %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 4',
              subtitle = glue::glue( 'The Ratio (across {n_estimators} estimators) of the Largest to Smallest Estimates of SE(β)\namong all estimators targeting a given estimand for {n_outcomes} outcomes' ),
              filename = 'final paper figures/fig4_SEratio_byestimand',
              byest = T,
              rmlegend = T,
              binw     = 0.1,
              coordratio = 2,
              tcratio  = T,
              tcbreaks = c(1:4))





## SE plots excluding fully interacted models ----

n_est_noFI = combined %>% 
  filter( !grepl('-FI', Estimator) ) %>%
  pull( Estimator ) %>%
  n_distinct()
n_est_noFI

summary_stats_exclFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 3a',
              subtitle = glue::glue( 'The Ratio (across {n_est_noFI} estimators) of the Largest to Smallest\nEstimates of SE(β) for {n_outcomes} outcomes\n(excluding fully interacted models)' ),
              filename = 'final paper figures/fig3a_SEratio_exclFI',
              rmlegend = T,
              binw     = 0.05,
              coordratio = 0.8)

### Figure 5: Figure 4 without fully-interactive models ----
summary_stats_estimand_exclFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 5',
              subtitle = glue::glue( 'The Ratio (across {n_est_noFI} estimators) of the Largest to Smallest Estimates of SE(β)\namong all estimators targeting a given estimand for {n_outcomes} outcomes\n(excluding fully interacted models)' ),
              filename = 'final paper figures/fig5_SEratio_byestimand_exclFI',
              byest = T,
              rmlegend = T,
              binw     = 0.1,
              coordratio = 2.2,
              tcratio  = T,
              tcbreaks = c(1:5))





