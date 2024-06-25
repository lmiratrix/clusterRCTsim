#
# Run a single scenario
#
# This will provide a nested simulation with R superpopulation
# replications.
#
# Output is a tibble, one row per method per simulation run.

source( "data_generation_code.R" )
source( "simulation_support_code_crt.R" )


# This will run a single scenario!
rr = run.scenario( K = 10, 
                   J.bar = 10,
                   p.bar = 0.5,
                   ATE.district = 0.2, 
                   sd.district = 0.1,
                   ICC.district = 0.2,
                   variable.J = F,
                   variable.p = T,
                   
                   # cluster-level params
                   n.bar = 50, 
                   sd.school = 0.1, 
                   ICC.school = 0.2,
                   variable.n = T,
                   
                   # size and correlation params
                   size.impact.cor.district = 0,
                   size.intercept.cor.district = 0.5,
                   size.impact.cor.school = 0,
                   proptx.impact.cor.district = 0,
                   
                   R = 2,
                   .progress="text",
                   include_method_characteristics=FALSE)


rr = as_tibble( rr )
rr
table( rr$run, rr$subrun )

head( rr )  
