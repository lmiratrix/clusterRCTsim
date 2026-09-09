
# "Fingerprint" simulations are simulations calibrated to the
# structure of a given dataset, or set of descriptive statistics
# regarding a specific dataset.
#
# They usually have no treatment effect, and are more about much
# variation we would see across estimators when there is in fact no
# variation in the treatment effect, and thus in the possible
# estimands one might consider.

# This script demonstrates running such simulations for a series of
# examples.


source( "R/calibrated_simulation_code.R" )


# ---- First simulation ----

# This makes data purely based on the descriptive statistics.

studies = read_csv( here::here( "empirical/results/study_descriptions.csv" ),
                    show_col_types = FALSE )
studies

stud = filter( studies, UniqueID == "CARES_TR_BPI_Incred" )
clusterRCT:::print.clusterRCTstats( stud )

# Demo of DGP
dat = generate_data_from_fingerprint( structure = NULL,
                                      parameters = stud )  
dat

# Running simulation itself
res <- run_fingerprint_simulation( structure = NULL,
                                   parameters = stud, 
                                   R = 100,
                                   include_MLM = FALSE,
                                   summarize_results = "cross", parallel = FALSE )
res$ratioSE[ res$ratioSE > 5 ] = 5

res %>% group_by( weight ) %>%
  summarise( mn_rangeATE = mean(rangeATE),
             sd_rangeATE = sd(rangeATE),
             mn_ratioSE = mean( ratioSE ),
             sd_ratioSE = sd( ratioSE ),
             per_max = mean( ratioSE == 5 ) )
             



# Second simulation (simulation code check) ----

# Run permutation simulation on fixed dataset

# TODO: See if this compares to just generating data over and over?
# We should get the same pattern of results.

set.seed( 1039 )
library( tidyverse )
library( clusterRCT )


ATE = 0.125
tx_het = TRUE

model.params.list <- list(
  M = 1                            # number of outcomes
  , J = 30                          # number of schools
  , K = 10                          # number of districts
  , nbar = 10                       # number of individuals per school
  , S.id = NULL                     # N-length vector of school assignments
  , D.id = NULL                     # N-length vector of district assignments
  , Xi0 = 0                         # scalar grand mean outcome under no treatment
  , MDES = ATE            # minimum detectable effect size
  , R2.3 = 0.1              # percent of district variation
  , ICC.3 = 0.2             # district intraclass correlation
  , omega.3 = 0.2*tx_het           # ratio of district effect size variability
  , R2.2 = 0.1              # percent of school variation
  , ICC.2 = 0.2             # school intraclass correlation
  , omega.2 = 0.0*tx_het          # ratio of school effect size variability
  , R2.1 = 0.1    # percent of indiv variation explained
)


sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
head( sim.data )


raw_res2 <- run_rerandomize_simulation( Yobs ~ T.x | S.id | D.id, 
                                        data=sim.data,
                                        R = 100,
                                        include_MLM = FALSE )
head( raw_res2 )

sim_res2 <- summarize_simulation_results( raw_res2 )
head( sim_res2 )

ggplot( sim_res2, aes( method ) ) +
  geom_point( aes( y = ESEhat, color = "SE_hat") ) +
  geom_errorbar( aes( ymin = Q10, ymax = Q90, color = "SE_hat" ), width = 0 ) +
  geom_point( aes( y = SE, color = "true SE" ), size = 2 ) +
  coord_flip()

sim_res2
ggplot( sim_res2, aes( method ) ) +
  geom_point( aes( y = EATE, color = "ATE_hat") ) +
  geom_errorbar( aes( ymin = EATE+2*SE, ymax = EATE-2*SE, color = "ATE_hat" ), 
                 width = 0 ) +
  geom_point( aes( y = 0, color = "true ATE" ), size = 2 ) +
  coord_flip()



sL <- sim_res2 %>%
  pivot_longer( cols = c(EATE, ESEhat, sdSEhat, SE), 
                names_to="metric", values_to="value" )
sL

ggplot( sL, aes( method, value, pch=weight )) +
  facet_wrap( ~ metric, scales="free_x" ) +
  geom_hline( yintercept = 0 ) +
  geom_point() +
  coord_flip()




