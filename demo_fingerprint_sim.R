
# Illustrate making and using grids in our simulation framework




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



# Run the simulation ----

if ( FALSE ) {
  sing <- single.CRT.fixed.district( dgp_function = generate_data_from_fingerprint,
                                     dgp_args = list( structure = struc, parameters = dd,
                                                      names = c( "Yobs", "Z","sid", "did" ) ),
                                     n.runs = 3 )
}

sim <- run.scenario.fixed.district( 
  dgp_function = generate_data_from_fingerprint,
  dgp_args = list( structure = struc, parameters = dd, tau = 0.20,
                   names = c( "Yobs", "Z","sid", "did" ) ),
  #just.describe.data = TRUE,
  R = 10, n.runs = 3, parallel = FALSE )

sim

rngs <- sim %>% group_by( subrun ) %>%
  summarise( range = diff( range( ATE.hat ) ) )
rngs

ggplot( rngs, aes( range ) ) +
  geom_dotplot( binwidth = 0.0025 )

summary( rngs$range )

