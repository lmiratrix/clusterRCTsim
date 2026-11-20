
# NOTE: This version (stripped) makes figures without the titles and
# captions, for the final paper layout.



# Make all the figures for the main empirical paper (for the blocked
# datasets)
#
# Colin, edited by Luke
# 2024


# load data ----

# load all data
source( here::here( "empirical/load_empirical_results.R" ) )
source( here::here( "empirical/R/makeDotplot.R" ) )


if ( FALSE ) {
  # filter to curated studies
  combined <- combined %>%
    filter(Collection == 'Curated')
}


source( here::here( "empirical/calculate_all_ranges.R" ) )



# Make Paper Figures ----


## Figure 1: The Range (across 11 Estimators) of Estimates of β for {n_outcomes} Outcomes ----
n_outcomes <- n_distinct( summary_stats$UniqueID )
n_outcomes
n_estimators <- n_distinct( combined$Estimator )
n_estimators

range_x_range <- 0.30

ss_best <- summary_stats %>%
  filter(Covariates == "best")
ss_best
table( ss_best$ATE_N )
ss_best$ATE_under = ifelse( ss_best$ATE_N < max(ss_best$ATE_N ), "yes", "no" )

plt <- ss_best %>%
  makeDotplot(var      = 'Range',
              x_mark_var = "ATE_under",
              xlab     = 'Range of Point Estimates of β', 
              ylab     = '', 
              binw = 0.01,
              coordratio = 0.1,
              rmlegend = FALSE) +
  scale_x_continuous( limits=c(0, range_x_range ), breaks = seq( 0, range_x_range, 0.05 ) ) +
  legend_inside()
plt

my_ggsave( plt,               
           filename   = 'final paper figures/fig1_ATErange',
           dim = c( 7, 3 ) )




## Figure 2: The Range of Estimates of β for all Estimators Targeting a Given Estimand  ----
ss_best_estimand = summary_stats_estimand %>%
  filter(Covariates == "best") %>%
  group_by( Estimand ) %>%
  mutate( ATE_under = ifelse( ATE_N < max(ATE_N ), "yes", "no" ) ) %>%
  ungroup()

plt <- ss_best_estimand %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β',
              ylab     = '',
              x_mark_var = "ATE_under",
              coordratio = .15,
              binw = 0.01,
              rmlegend = F,
              byest    = TRUE,
              xscale = c(0, 0.3))  +
  scale_x_continuous( limits=c(0, range_x_range ), breaks = seq( 0, range_x_range, 0.1 ) )+
  legend_inside()
plt

my_ggsave( plt,               
           filename = 'final paper figures/fig2_ATErange_byestimand',
           dim = c( 7, 3 ) )





## Figure 3: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) for {n_outcomes} outcomes. ----

table( ss_best$SE_N )

ss_best$SE_under = ifelse( ss_best$SE_N < max(ss_best$SE_N ), "yes", "no" )


plt <- ss_best %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = '',
              x_mark_var = "SE_under",
              rmlegend = FALSE,
              binw     = 0.1,
              coordratio = 1.25,
              tcratio  = T,
              tcbreaks = c(1:4))+
  legend_inside()
plt

my_ggsave( plt,               
           filename = 'final paper figures/fig3_SEratio',
           dim = c( 7, 3 ) )



## Figure 4: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) among all estimators targeting a given estimand for {n_outcomes} outcomes ----

ss_best_estimand = summary_stats_estimand %>%
  filter(Covariates == "best") %>%
  group_by( Estimand ) %>%
  mutate( SE_under = ifelse( SE_N < max(SE_N), "yes", "no" ) ) %>%
  ungroup()
table( ss_best_estimand$SE_under )

plt <- ss_best_estimand %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = '',
              x_mark_var = "SE_under",
              byest = T,
              rmlegend = T,
              binw     = 0.1,
              coordratio = 2,
              tcratio  = T,
              tcbreaks = c(1:4))+
  legend_inside()
plt

my_ggsave( plt,               
           filename = 'final paper figures/fig4_SEratio_byestimand',
           dim = c( 7, 3 ) )






