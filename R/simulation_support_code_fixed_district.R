

# A version of simulation_support_code_crt.R that allows for a fixed
# set of districts
#
# Has extra methods for taking a 'grid' of district level
# characteristics and generating data.


library( tidyverse )
library( purrr )
library( multiwayvcov )

library( clusterRCT )


source( here::here( "R/util.R" ) )
source( here::here( "R/data_generation_code.R" ) )





#' Run a single round of simulation.
#'
#' Generate a cluster-randomized trial dataset, and then repeatedly
#' rerandomize the fixed dataset n.runs times, each time estimating
#' the ATE using all our methods. It also calculates the true
#' finite-sample ATE (person and site weighted).
#'
#' There are two ways to generate data: via the grid or via
#' dgp_function (with possible additional dgp_args).  
#' 
#' @inheritParams run.scenario.fixed.district
#'
#' @return Dataframe with the estimates along with the true baseline
#'   values.
single.CRT.fixed.district = function(
    grid = NULL,
    
    # method params
    n.runs = 3,
    just.describe.data = FALSE, 
    just.return.data = FALSE,
    dgp_function = NULL,
    dgp_args = NULL,
    ...,
    .progress = FALSE ) {
  
  # Make the dataset
  stopifnot( is.null( grid ) || is.null( dgp_function ) )
  if ( !is.null( dgp_function ) ) {
    df <- do.call( dgp_function, dgp_args )
  } else if ( !is.null( grid ) ) {
    df <- generate_data_for_districts( grid = grid )
  } else {
    stop( "Need to provide either a grid or a data generating function." )
  }
  
  if ( just.return.data ) {
    return( df )
  }
  
  if ( just.describe.data ) {
    return( describe.data( df, print = FALSE ) )
  }
  
  # Calculate true finite-sample estimands
  beta_jks <- df %>% 
    group_by(did, sid) %>% 
    summarize(n_jk = n(), 
              beta_jk = mean(Y1-Y0), .groups="drop")
  
  ATE.student <- with( df, mean( Y1 - Y0 ) )
  ATE.school <- mean(beta_jks$beta_jk)
  
  ATE.dist.student <- beta_jks %>% 
    summarize(n_k = sum(n_jk),
              s = sum(n_jk*beta_jk)) %>% 
    summarize(ATE.dist.student = sum(s/n_k) / n()) %>% 
    pull(ATE.dist.student)
  ATE.dist.school <- beta_jks %>% 
    summarize(J_k = n(),
              s = sum(beta_jk)) %>% 
    summarize(ATE.dist.school = sum(s/J_k) / n()) %>% 
    pull(ATE.dist.school)
  
  # Do finite-sample simulation, rerandomizing n.runs times.
  ests = map_dfr( 1:n.runs, function( x ) {
    df = rerandomize.data( df )
    get.estimates( df, ... )
  }, .id="subrun", .progress=.progress )
  
  
  # Gather results and ship!
  ests = mutate( ests,
                 ATE.finite.student = ATE.student,
                 ATE.finite.school = ATE.school,
                 ATE.finite.dist.student = ATE.dist.student,
                 ATE.finite.dist.school = ATE.dist.school ) %>%
    as_tibble()
  
  if (FALSE) {
    ests %>% filter(subrun==1)
    
    # dists of tx and co
    df %>% 
      group_by(did, sid) %>% 
      summarize(Z = mean(Z),
                mn = mean(Yobs)) %>% 
      ggplot(aes(x=mn, group=Z, color=as.factor(Z))) +
      geom_density()
    
    df %>% 
      group_by(Z) %>% 
      summarize(mn = mean(Yobs))
  }
  
  ests    
}



if ( FALSE ) {
  
  grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.ATE = 0.3,
                           ICC = 0.2,
                           variable.J = TRUE, variable.p = TRUE, variable.nbar = TRUE,
                           num.covariates = 3, size.intercept.cor = -0.8,
                           size.impact.cor = 0.75, proptx.impact.cor = 0.75, nbar.J.cor = -0.9, 
                           empirical = TRUE ) %>%
    arrange( J )
  
  grid
  grid$J[1] = 4
  grid$n.bar[1] = 100
  grid$district.tau[1] = -2
  #ggplot( grid, aes( n.bar, J ) ) +
  #  geom_point()
  cor( grid$J, grid$n.bar )  
  describe_district_grid( grid )
  
  ac = single.CRT.fixed.district( grid = grid, 
                                  n.runs = 50, 
                                  just.describe.data = TRUE )
  ac
  
  
  
  ac = single.CRT.fixed.district( grid = grid, 
                                  n.runs = 3 )
  ac
  
  ac = single.CRT.fixed.district( grid = grid, 
                                  n.runs = 3, just.return.data = TRUE )
  ac
  
  get.estimates( ac )
  
  
  a = single.CRT.fixed.district( grid = grid, 
                                 n.runs = 100, 
                                 .progress = TRUE )
  table( a$subrun )
  a  
  
  a %>% 
    filter( disfavored != 1 ) %>%
    group_by( method, weight ) %>%
    summarise( SE.b = sd( ATE.hat ) / sqrt( n() ),
               bias.stu = mean( ATE.hat - ATE.finite.student ),
               bias.sch = mean( ATE.hat - ATE.finite.school ) ) %>%
    arrange( weight ) %>%
    knitr::kable( digits = 2 )
  
  
}





#' Run a simulation with R trials and return all the simulation runs
#' as a large dataframe.
#'
#' This function takes a given grid (or, alternatively, a data
#' generating process) and runs a simulation using it.
#'
#' @param grid Dataframe, one row per district, of DGP parameters.
#' @param dgp_function Function to generate dataset.  Will use if grid
#'   is NULL.
#' @param dgp_args List of arguments to pass to dgp_function, if there
#'   is one.
#' @param R Number of datasets to generate.
#' @param n.runs Number of "finite sample" replicates for each
#'   dataset.  So total of R x n.runs simulation runs will be done.
#' @param parallel Should the simulation runs be parallelized using
#'   the purrr package?
#'
#' @return Stacked dataframe of all the results of the
#'   compare_methods() calls.
run.scenario.fixed.district = function( 
    grid = NULL,
    #ICC.district = sd( grid$district.alpha ),
    
    # method params
    n.runs = 3,
    R = 10,
    parallel = FALSE, 
    dgp_function = NULL,
    dgp_args = NULL,
    ID=NULL, .progress="none", ... ) {
  
  stopifnot( is.character(.progress) & length(.progress) == 1 )
  
  ptm = proc.time()
  
  if ( is.null( ID ) ) {
    ID = -999
  }
  
  scat( "Running scenario %s\n", ID )
  
  #describe_district_grid( grid )
  
  # scat( "Assumed ICC (district) = %.2f\nR = %d X %d\n", 
  #       ICC.district,
  #       R, n.runs )
  scat( "R = %d X %d\n", R, n.runs )
  
  # Get all the ... arguments in list
  dots = list(...)
  if ( length( dots ) == 0 ) {
    cat( "\tNo extra arguments passed.\n" )
  } else {
    scat( "\tExtra arguments: %s\n\n", paste(names(dots), collapse=", ") )
  }
  
  single_shot <- function( i, ... ) { 
    single.CRT.fixed.district( 
      grid = grid,
      dgp_function = dgp_function,
      dgp_args = dgp_args,
      n.runs=n.runs, ... )
  }
  
  if ( parallel ) {
    library( furrr )
    plan(multisession, workers = parallel::detectCores() - 1 )
    
    rps = future_map_dfr( 1:R, single_shot, ..., 
                          .id="run", 
                          .options = furrr_options(seed = TRUE),
                          .progress = TRUE )
  } else {
    rps = map_dfr( 1:R, single_shot, ..., 
                   .id="run", .progress=.progress )
  }
  
  if ( "subrun" %in% colnames( rps ) ) {
    rps$subrun = paste0( rps$run, ".", rps$subrun )
  }

  scat("**\n**\tTotal time elapsed:\n")
  tot.time = proc.time() - ptm
  print(tot.time)
  sim.per.min = R / (tot.time["elapsed"] / 60)
  scat("Simulations per minute = %.2f\n", sim.per.min )
  rps$sim.per.min = sim.per.min
  rps
}






generate_clusterCRT <- function(J, n, ICC, tau, cor.ICC.tau, p) {
  
  # Set mean vector and covariance matrix for bivariate normal distribution
  mean_vec <- c(0.2, 0)
  cov_matrix <- matrix(c(tau, cor.ICC.tau * sqrt(tau * ICC), 
                         cor.ICC.tau * sqrt(tau * ICC), ICC), 
                       nrow = 2)
  
  # Generate cluster-level random effects and cluster-average outcomes
  cluster_random_effects <- MASS::mvrnorm(J, mu = mean_vec, Sigma = cov_matrix)
  
  # Extract cluster effects and average outcomes
  cluster_effects <- cluster_random_effects[, 1]
  cluster_avg_outcomes <- cluster_random_effects[, 2]
  
  # Generate individual-level data for each cluster
  data <- tibble(
    sid = rep(1:J, each = n),
    Z = rep(0 + (sample(J) <= p*J), each = n),
    Y0 = rep(cluster_avg_outcomes, each = n) +
      rnorm(J * n, mean = 0, sd = sqrt(1-ICC)),
    Y1 = Y0 + rep(cluster_effects, each = n),
    Yobs = ifelse( Z, Y1, Y0 )
  )
  data$did = 1
  data$ATE.finite.student = data$ATE.finite.school = mean( cluster_effects )
  data
}



if ( FALSE ) {
  
  generate_clusterCRT( 10, 20, 0.2, 0.3^2, 0.75, 0.5 )
  
  
  dd_check =   generate_clusterCRT( J = 9, n = 500, ICC = 0, tau = 0.3^2, cor.ICC.tau = -1, p = 0.5 )
  
  
  undebug( single.CRT.fixed.district )
  single.CRT.fixed.district( dgp_function = generate_clusterCRT, 
                             dgp_args = list( J = 10, n = 20, ICC = 0.2, tau = 0.3^2, cor.ICC.tau = 0.75, p = 0.5 ) )
  
  run.scenario.fixed.district( dgp_function = generate_clusterCRT, 
                    dgp_args = list( J = 10, n = 20, ICC = 0.2, tau = 0.3^2, cor.ICC.tau = 0.75, p = 0.5 ),
                    n.runs = 3, R = 10 )
}




# Function to run and save simulation and simulation grid ----


#' Run a simulation and save results to "results" directory
#'
#'
#' @inheritParams run.scenario.fixed.district
#' @param name Name of the simulation (will cast this to character),
#'   used in filenames.
#'   
#' @return The dataframe of simulation results.
run_sim <- function( name = "", 
                     R, n.runs, 
                     grid = NULL, 
                     dgp_function = NULL, dgp_args = NULL, 
                     parallel = TRUE ) {
  
  name = as.character(name)
  
  if ( !is.null( grid ) ) {
    saveRDS( grid, file=paste0( "results/superpop_SE_simulation_grid", 
                                name, ".rds" ) )
  }
  
  cat( "\n\n" )
  cli::cli_alert_info( "Running simulation {name}\n" )
  #options( warn = 2 )
  sim <- run.scenario.fixed.district( 
    grid,
    dgp_function = dgp_function,
    dgp_args = dgp_args,
    ID = name,
    # method params
    n.runs = n.runs,
    R = R,
    .progress="Rs",
    parallel = parallel )
  
  if ( !is.null(grid) ) {
    attr( sim, "grid" ) <- grid
  }
  
  saveRDS( sim, file=glue::glue( "results/superpop_SE_simulation{name}.rds" ) )
  cli::cli_alert_info( "Simulation {name} complete and saved.\n" )
  
  invisible( sim )
}



# Test code -----


if ( FALSE ) {
  
  grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.ATE = 0.1,
                           ICC = 0.2,
                           variable.J = TRUE, variable.p = TRUE, variable.nbar = TRUE,
                           num.covariates = 3, size.intercept.cor = -0.8,
                           size.impact.cor = 0.75, proptx.impact.cor = 0.75, nbar.J.cor = 0.9, 
                           empirical = TRUE )
  
  grid
  describe_district_grid( grid )
  
  a = single.CRT.fixed.district( grid = grid )
  a
  
  df = run.scenario.fixed.district( grid, n.bar = 50, sd.school = 0.1, 
                                    size.impact.cor.school = 0.75,
                                    variable.n = TRUE,
                                    n.runs = 3 )  
  
  df
}

