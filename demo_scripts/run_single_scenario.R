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
                         ATE = 0.2, sd.ATE = 0.4,
                         ICC.district = 0.2,
                         variable.p = TRUE,
                         cv.district = 0.6,
                         cv.nbar = 0.6,
                         size.intercept.cor = -0.8,
                         R2 = NULL,
                         num.cov.district = 0,
                         sd.school = 0.4,
                         ICC.school = 0.2, 
                         R2.school = NULL,
                         num.cov.school = 0,
                         cv.school = 0.6,
                         size.impact.cor.school = 0.5,
                         R2.indiv = NULL,
                         num.cov.indiv = 0,
                         empirical = FALSE )


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
