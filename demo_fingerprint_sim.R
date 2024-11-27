
# Illustrate making and using the "fingerprint" DGP process in our
# simulation framework

# In particular, this script shows how we can load the descriptive
# statistics and the skeleton of a dataset and run a simulation based
# on those measured statistics.




source( "R/data_generation_code.R" )
source( "R/simulation_support_code_fixed_district.R" )

library( lme4 )
library( lmerTest )
library( future )
library( clusterRCT )

set.seed( 404404 )

# Get some original data ----
data( fakeCRT )
head( fakeCRT )
fakeCRT$X2 = rnorm( nrow(fakeCRT) )



# Analyze our original data to get DGP parameters ----

# Describe it to get our DGP stats
dd = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                          control_formula = ~ X.jk + C.ijk + X2 )
dd

# Get the structure of it as well (cluster sizes, etc)
struc <- get_structure( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
struc




# We can generate fake data now ----
fakeStar <- generate_data_from_fingerprint( struc, dd, tau = 0.20,
                                            names = c( "Yobs", "Z","sid", "did" ) )
head( fakeStar )

star_stat <- describe_clusterRCT( Yobs ~ Z | sid | did, 
                                  control_formula = ~ X1 + X2 + W,
                                  data=fakeStar )
star_stat

# Compare to original!
bind_rows( dd, star_stat )

# Our treatment effect is a constant 0.20.  No heterogeniety in
# treatment effect.
compare_methods( Yobs ~ Z | sid | did, 
                 control_formula = ~ X1 + X2 + W,
                 data=fakeStar, include_disfavored = FALSE, 
                 include_method_characteristics = FALSE )


# Here we generate data from a file ----

# We need the skeleton and the describe data results stored in some
# result file

# load skeleton
skeleton <- readRDS( "empirical/results/beaconmentor_J-K-T.RDS" )
names(skeleton) <- c( "clusterID", "blockID", "Z" )

# load stats
stats <- readRDS( "empirical/results/master240528_combined.RDS" )
head( stats )
ss <- filter( stats, UniqueID == "Beacon Mentor_S8TGPA",
              Estimator == "Agg_FE_Cluster", 
              Covariates == "best" 
) %>%
  as_tibble()
ss
if ( FALSE ) {
  ss$block_ICC = 0
  ss$cluster_ICC = 0.2
  ss$R2.1 = 0
  ss$R2.2 = 0.20
}
fakeStar <- generate_data_from_fingerprint( skeleton, ss, tau = 0.20,
                                            names = c( "Yobs", "Z","sid", "did" ) )
fakeStar
#nrow( fakeStar )
#nrow( skeleton )

star_stat <- describe_clusterRCT( Yobs ~ Z | sid | did, 
                                  control_formula = ~ X1 + X2,
                                  data=fakeStar )

bind_rows( ss, star_stat ) %>%
  dplyr::select( -c( UniqueID:Trt_arm, Pretest:Estimand ) ) %>%
  as.data.frame()



# Run the simulation ----

if ( FALSE ) {
  sing <- single.CRT.fixed.district( dgp_function = generate_data_from_fingerprint,
                                     dgp_args = list( structure = struc, parameters = dd,
                                                      names = c( "Yobs", "Z","sid", "did" ) ),
                                     n.runs = 3 )
}

sim <- run.scenario.fixed.district( 
  dgp_function = generate_data_from_fingerprint,
  dgp_args = list( structure = skeleton, parameters = ss, tau = 0.20,
                   names = c( "Yobs", "Z","sid", "did" ) ),
  #just.describe.data = TRUE,
  R = 10, n.runs = 1, parallel = FALSE )

sim

rngs <- sim %>% group_by( subrun ) %>%
  summarise( range = diff( range( ATE.hat ) ) )
rngs

ggplot( rngs, aes( range ) ) +
  geom_dotplot(  )

summary( rngs$range )


# Run the simulation without the skeleton ----

# We can generate data on the fly if we don't pass the skeleton

# load stats
stats <- readRDS( "empirical/results/master240528_combined.RDS" )
head( stats )
ss <- filter( stats, UniqueID == "Beacon Mentor_S8TGPA",
              Estimator == "Agg_FE_Cluster", 
              Covariates == "best" 
) %>%
  as_tibble()

fakeStar <- generate_data_from_fingerprint( parameters = ss, tau = 0.20,
                                            names = c( "Yobs", "Z","sid", "did" ) )
head( fakeStar )

# Number of observations now only approximate
nrow( fakeStar )
nrow( skeleton )

star_stat <- describe_clusterRCT( Yobs ~ Z | sid | did, 
                                  control_formula = ~ X1 + X2,
                                  data=fakeStar )

bind_rows( ss, star_stat ) %>%
  dplyr::select( -c( UniqueID:Trt_arm, Pretest:Estimand ) ) %>%
  as.data.frame()


# And now the simulation

sim <- run.scenario.fixed.district( 
  dgp_function = generate_data_from_fingerprint,
  dgp_args = list( structure = NULL, parameters = ss, tau = 0.20,
                   names = c( "Yobs", "Z","sid", "did" ) ),
  #just.describe.data = TRUE,
  R = 10, n.runs = 1, parallel = FALSE )

sim

rngs <- sim %>% group_by( subrun ) %>%
  summarise( range = diff( range( ATE.hat ) ) )
rngs

ggplot( rngs, aes( range ) ) +
  geom_dotplot(  )

summary( rngs$range )


