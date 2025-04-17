#
# Run a single simulation scenario
#
# This runs a nested simulation with R superpopulation (at the cluster
# level) replications.
#
# Output is a tibble, one row per method per simulation run.

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_district.R" ) )


# This will run a single scenario!

grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                         ATE = 0.2, sd.ATE = 0.1,
                         ICC = 0.2,
                         variable.J = TRUE, variable.p = TRUE, variable.nbar = TRUE,
                         num.cov.school = 3, size.intercept.cor = -0.8,
                         size.impact.cor = 0.75, proptx.impact.cor = 0.75, nbar.J.cor = 0.9, 
                         empirical = TRUE )

grid
describe_district_grid( grid )


rr = run.scenario.fixed.district( grid,
                   R = 2,
                   .progress="text",
                   include_method_characteristics=FALSE)


rr = as_tibble( rr )
rr
summary( rr$ATE.hat )
table( rr$run, rr$subrun )

head( rr )  
