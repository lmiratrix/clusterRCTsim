

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



# Summarize results for aggregate plots  ----

## Summarize ATE estimates to calculate ranges ----
# summarize ATE by study-outcome
ate_summary <- combined %>%
  group_by(UniqueID, class, low_df, has_single, Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = "drop") %>%
  mutate(Range = Max - Min)

# summarize ATE by study-outcome-estimand
ate_summary_estimand <- combined %>%
  group_by(Collection, class, low_df, has_single, Project, Trt_arm, Outcome, Covariates, Estimand) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = 'drop') %>%
  mutate(Range = Max - Min)

## Summarize SEhat estimates to calculate ranges ----

se_summary <- combined %>%
  group_by(Collection, class, low_df, has_single, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

se_summary_estimand <- combined %>%
  group_by(Collection, class, low_df, has_single, Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

## SE ranges (excluding fully interacted models) ----

se_summary_exclFI <- combined %>%
  filter(!grepl('-FI', Estimator)) %>%
  group_by(Collection, class, low_df, has_single, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

se_summary_estimand_exclFI <- combined %>%
  filter(!grepl('-FI', Estimator)) %>%
  group_by(Collection, class, low_df, has_single, Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )





# Make Paper Figures ----
## Figure 1: The Range (across 11 Estimators) of Estimates of β for {n_outcomes} Outcomes ----
table( ate_summary$class )
nrow( ate_summary )

range( ate_summary$Range )
n_outcomes <- n_distinct( ate_summary$UniqueID )
n_estimators <- n_distinct( combined$Estimator )
n_estimators
nrow( ate_summary )
table( table( ate_summary$UniqueID ) )

range_x_range <- 0.30

ate_summary %>%
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
ate_summary_estimand %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β',
              ylab     = 'Number of Outcomes',
              filename = 'final paper figures/fig2_ATErange_byestimand',
              title = 'Figure 2',
              subtitle = glue::glue( 'The Range of Estimates of β for all Estimators\nTargeting a Given Estimand for {n_outcomes} Outcomes' ),
              coordratio = .10,
              binw = 0.01,
              rmlegend = F,
              byest    = TRUE ) +
  scale_x_continuous( limits=c(0, range_x_range ) )



## Figure 3: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) for {n_outcomes} outcomes. ----
se_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 3',
              subtitle = glue::glue( 'The Ratio (across {n_estimators} estimators) of the Largest to Smallest\nEstimates of SE(β) for {n_outcomes} outcomes' ),
              filename = 'final paper figures/fig3_SEratio',
              rmlegend = FALSE,
              binw     = 0.1,
              coordratio = 1,
              tcratio  = T,
              tcbreaks = c(1:5))



## Figure 4: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) among all estimators targeting a given estimand for {n_outcomes} outcomes ----
se_summary_estimand %>%
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
              coordratio = 1.5,
              tcratio  = T,
              tcbreaks = c(1:5))

## SE plots excluding fully interacted models ----

n_est_noFI = combined %>% 
                           filter( !grepl('-FI', Estimator) ) %>%
                           pull( Estimator ) %>%
                           n_distinct()
n_est_noFI

se_summary_exclFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 3a',
              subtitle = glue::glue( 'The Ratio (across {n_est_noFI} estimators) of the Largest to Smallest\nEstimates of SE(β) for {n_outcomes} outcomes\n(excluding fully interacted models)' ),
              filename = 'fig3a_SEratio_exclFI',
              rmlegend = T,
              binw     = 0.05,
              coordratio = 0.8)

se_summary_estimand_exclFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 4a',
              subtitle = glue::glue( 'The Ratio (across {n_est_noFI} estimators) of the Largest to Smallest Estimates of SE(β)\namong all estimators targeting a given estimand for {n_outcomes} outcomes\n(excluding fully interacted models)' ),
              filename = 'fig4a_SEratio_byestimand_exclFI',
              byest = T,
              rmlegend = T,
              binw     = 0.05,
              coordratio = 0.8)




# And figures for the non-blocked datasets ----


