#
# Run a single simulation scenario
#
# This runs a nested simulation with R superpopulation (at the cluster
# level) replications.
#
# Output is a tibble, one row per method per simulation run.

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_block.R" ) )


# This will run a single scenario!

grid <- build_blocks( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                         ATE = 0.2,
                         ICC.block = 0.2,
                         variable.p = TRUE,
                         cv.block = 0.6,
                         cv.nbar = 0.6,
                         size.intercept.cor = -0.8,
                         R2 = NULL,
                         num.cov.block = 0,
                         sd.cluster = 0.4,
                         ICC.cluster = 0.2, 
                         R2.cluster = NULL,
                         num.cov.cluster = 0,
                         cv.cluster = 0.6,
                         size.impact.cor.cluster = 0.5,
                         R2.indiv = NULL,
                         num.cov.indiv = 0,
                         empirical = FALSE )


grid
describe_block_grid( grid )


rr = run.scenario.fixed.block( grid,
                   R = 2,
                   .progress="text",
                   include_method_characteristics=FALSE)

rr

summary( rr$ATE.hat )
table( rr$run, rr$subrun )

head( rr )  
