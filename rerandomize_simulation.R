

library( tidyverse )
source( "R/calibrated_simulation_code.R" )

data( fakeCRT )
head( fakeCRT )

#formula = Yobs ~ T.x | S.id | D.id
#data = fakeCRT

sim_res <- run_rerandomize_simulation( Yobs ~ T.x | S.id | D.id, 
                                       data=fakeCRT,
                                       R = 3,
                                       summarize_results = FALSE,
                                       include_MLM = FALSE )

# Raw results
sim_res

# Aggregate results if desired
# sum_res <- process_results( sim_res )
sum_res

ggplot( sim_res, aes( method, SE ) ) +
  geom_point( aes( y = ESEhat, color = "SE_hat") ) +
  geom_errorbar( aes( ymin = Q10, ymax = Q90, color = "SE_hat" ), width = 0 ) +
  geom_point( aes( color = "true SE" ), size = 2 ) +
  coord_flip()

sL <- sim_res %>%
  pivot_longer( cols = c(EATE, ESEhat, sdSEhat, SE), names_to="metric", values_to="value" )
sL

ggplot( sL, aes( method, value, pch=weight )) +
  facet_wrap( ~ metric, scales="free_x" ) +
  geom_hline( yintercept = 0 ) +
  geom_point() +
  coord_flip()

sim_res



# Second simulation ----

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


sim_res2 <- run_rerandomize_simulation( Yobs ~ T.x | S.id | D.id, 
                                       data=sim.data,
                                       R = 100,
                                       include_MLM = FALSE )
sim_res2

ggplot( sim_res2, aes( method, SE ) ) +
  geom_point( aes( y = ESEhat, color = "SE_hat") ) +
  geom_errorbar( aes( ymin = Q10, ymax = Q90, color = "SE_hat" ), width = 0 ) +
  geom_point( aes( color = "true SE" ), size = 2 ) +
  coord_flip()

sL <- sim_res2 %>%
  pivot_longer( cols = c(EATE, ESEhat, sdSEhat, SE), names_to="metric", values_to="value" )
sL

ggplot( sL, aes( method, value, pch=weight )) +
  facet_wrap( ~ metric, scales="free_x" ) +
  geom_hline( yintercept = 0 ) +
  geom_point() +
  coord_flip()

sim_res2





