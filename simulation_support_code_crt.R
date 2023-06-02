
# temporary copy of simulation_support_code.R
#  so that I can easily see the original simulation_support_code 
#  as I adapt it to multisite trials

# mapping:
#  - single.MLM.trial <=> single.CRT
# updated:
#  - rerandomize.data: now re-cluster-randomizes data

# This script has all the support code to run a large simulation to look at
# stability of the different SE estimators and true performances of the ATE
# estimators.
#
# The simulation script runs a bunch of scenarios and saves the massive amount
# of results to a file to be processed.
#
# This script has the code to run a single trial, etc.

library( dplyr )
library( blkvar )
library( purrr )
library( multiwayvcov )

library( clusterRCT )


#### Utility functions #####

if ( FALSE ) {
  # for testing crashing in the simulator.  This wrapper will make things crash sometimes.
  my_compare_methods = function( ... ){
    if ( sample(5,1) == 1 ) {
      asdfadfga
    } else {
      compare_methods( ... )
    }
  }
}


scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}


get.estimates = function( df, ... ) {
  safe_comp = safely( compare_methods )
  ests = safe_comp( Yobs, Z, sid, data=df, include_block=FALSE, ... )
  if ( is.null( ests[[2]] ) ) {
    ests = ests[[1]]
    ests = rename( ests, 
                   ATE.hat = ATE_hat,
                   SE.hat = SE )
    
  } else {
    cat( "\nCaught error:\n" )
    print( ests[[2]] )
    ests = data.frame( ATE.hat = c(NA),
                       SE.hat = c(NA) )    
  }
  ests
}



# Reassign treatment within a district (to do finite sample inference)
rerandomize.data = function( dat, p ) {
  stopifnot(length(unique(did)) == 1)
  dat %>%
    mutate(
      # sid = droplevels(sid),   # avoid empty factor level issues with clusters=sid
      Z = randomizr::cluster_ra(clusters = sid, prob = p),
      Yobs = ifelse( Z == 1, Y1, Y0 ) )
}
# Reassign treatment to all districts
rerandomize.data.full = function( dat ) {
  grid <- attr(dat, "district_key")
  
  dat <- dat %>% 
    split(~did) %>% 
    purrr::map2(grid$p, ~rerandomize.data(.x, .y)) %>% 
    purrr::list_rbind()
  
  attr(dat, "district_key") <- grid
  dat
}


##### Code to explore scenarios and sample datasets #####

# Utility function to describe a given data set to explore the data generation
# process
# TODO
describe.data = function( df ) {
  
  # we want:
  #  - true ATEs
  #  - true correlations from district_key
  #  - etc...
  
  df = rerandomize.data.full( df )
  
  
  # district params
  sstat_distr <- attr(df, "district_key") %>% 
    summarize(cor.J.tau_dist = cor(J, district.tau),
              cor.p.tau_dist = cor(p, district.tau),
              sd.tau_dist = sd(district.tau),
              sd.J_dist = sd(J),
              sd.p_dist = sd(p))
  
  # compute true ATE values
  ATEs <- df %>% 
    group_by(did, sid) %>% 
    summarize(Njk = n(),
              Bjk = mean(Y1-Y0),
              .groups="drop_last") %>% 
    summarize(Nk = sum(Njk),
              Jk = n(),
              ATE.distr = sum(Njk*Bjk) / Nk,
              ATE.clust = sum(Bjk),
              ATE.indiv = sum(Njk*Bjk),
              .groups="drop_last") %>% 
    summarize(ATE.distr = mean(ATE.distr),
              ATE.clust = sum(ATE.clust) / sum(Jk),
              ATE.indiv = sum(ATE.indiv) / sum(Nk))
  
  
  # compute design-based point ATE estimates
  #  - district-weighted, cluster-weighted, and individual-weighted
  ATEhats <- df %>% 
    group_by(did, sid) %>% 
    summarize(n.clust = n(),
              Z.clust = mean(Z),
              Yhat.clust = mean(Yobs),
              .groups="drop_last") %>% 
    summarize(J.k = length(unique(sid)),
              n.k = sum(n.clust),
              Bhat.k.clust = sum(Z.clust*Yhat.clust) / sum(Z.clust) -
                sum((1-Z.clust)*Yhat.clust) / sum(1-Z.clust),
              Bhat.k.indiv = sum(n.clust*Z.clust*Yhat.clust) / sum(Z.clust*n.clust) -
                sum(n.clust*(1-Z.clust)*Yhat.clust) / sum((1-Z.clust)*n.clust)) %>% 
    summarize(ATE.distr = mean(Bhat.k.clust),
              ATE.clust = sum(J.k*Bhat.k.clust) / sum(J.k),
              ATE.indiv = sum(n.k*Bhat.k.indiv) / sum(n.k))
  
  res <- bind_rows(ATEs, ATEhats) %>% 
    mutate(ATE = c("True ATEs", "Estimated ATEs"), .before=everything())
  
  print(sstat_distr)
  print(res)
}



# Testing code
if ( FALSE ) {
  library( tidyverse )
  
  df = generate_multilevel_data_no_cov( n.bar=200, J=50,
                                        tau.11.star = 0.1^2,
                                        ICC = 0.20,
                                        p = 0.70,
                                        variable.n = TRUE,
                                        variable.p = TRUE,
                                        size.impact.correlate = TRUE,
                                        proptx.impact.correlate = TRUE,
                                        finite.model = FALSE )
  
  # empirical sd of true treatment effects across clusters
  tau.S = attr( df, "tau.S")
  tau.S    
  
  sites = df %>% group_by( sid ) %>% 
    summarise( n = n(),
               p.Z = mean( Z ),
               ATE.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ),
               .groups = "drop")
  
  head( sites )
  
  qplot( sites$n, binwidth=10 )
  summary( sites$n )
  qplot( sites$p.Z, binwidth=0.025 )
  summary( sites$p.Z )
  qplot( sites$n, sites$ATE.hat )
  cor( sites$n, sites$ATE.hat )
  
  qplot( sites$p.Z, sites$ATE.hat )
  cor( sites$p.Z, sites$ATE.hat )
  
  
  head( df )  
  df = rerandomize_data(df)
  head( df )
  
  describe.data(df)
  describe_data(df)  
  
  
  get.estimates(df)
}


##### Data-generating process #####

# generate variable cluster sizes
#  - analogue to blkvar::block_distn()
cluster_distn <- function(K, J.bar, size.ratio) {
  stopifnot(J.bar > 4)
  nk <- blkvar::block_distn(J=K, n.bar=J.bar-4, size.ratio)
  nk[nk<4] <- 4
  
  nk
}

# Generate a cluster-randomized trial, with multiple districts
# 
# District average effects may be correlated with
#  district size (i.e., # clusters) and proportion of clusters treated
generate_districts <- function(
    # district-level params
    K, J.bar,
    ATE.district, sd.district,
    variable.J = FALSE,

    # cluster-level params
    n.bar, tau.11, ICC = 0.20,
    variable.n = FALSE,
    
    # size and correlation params
    size.ratio.clust = 0.6, 
    size.ratio.indiv = 0.6,
    size.impact.correlate.clust = 0,
    size.impact.correlate.indiv = 0,
    proptx.impact.correlate.clust = F,
    correlate.strength.clust = 0.75,
    correlate.strength.indiv = 0.75,
    ...) {
  # idea: induce noisy correlation by...
  #  1. sort (noised) random district effects
  #  2. order num. clusters & prop. tx clusters
  
  stopifnot(size.impact.correlate.clust %in% c(-1,0,1))
  stopifnot(size.impact.correlate.indiv %in% c(-1,0,1))
  
  # values for J
  if (variable.J) {
    num_clusters <- cluster_distn(K, J.bar, size.ratio=size.ratio.clust)
  } else {
    num_clusters <- rep(J.bar, K)
  }
  # values for proportion treated
  prop_tx_clusters <- runif(K, min=0.25, max=0.75)
  
  # district-level parameters
  if (size.impact.correlate.clust == 0) {
    grid <- tibble(
      district.tau = rnorm(K, mean=ATE.district, sd=sd.district),
      J = num_clusters,
      p = prop_tx_clusters
    ) %>% 
      mutate(did = 1:n(), .before=everything())
  } else {
    grid <- tibble(
      district.tau = rnorm(K, mean=ATE.district, sd=sd.district),
    ) %>%
      arrange(size.impact.correlate.clust * 
                # (gamma.10 + rnorm(K, mean=0, sd=eps))
                (correlate.strength.clust^2 * district.tau +
                   (1-correlate.strength.clust^2) * rnorm(K, sd=sd.district)) ) %>% 
      mutate(did = 1:n(), .before=everything(),
             J = sort(num_clusters),
             p = sort(prop_tx_clusters))
  }
  
  # check correlation
  if (F) {
    grid %>% 
      summarize(corJ = cor(district.tau, J, method="spearman"),
                corp = cor(district.tau, p, method="spearman"))
  }
  
  # make student-level dataset
  sim.data <- grid %>%
    rowwise() %>%
    mutate(sim_data = list(
      generate_multilevel_data_no_cov( 
        n.bar = n.bar, 
        J = J,
        p = 0.5,   # doesn't matter, since we rerandomize tx to clusters anyways
        tau.11 = tau.11,
        gamma.10 = district.tau,
        ICC = ICC,
        
        variable.n = variable.n,
        finite.model = FALSE,
        size.impact.correlate = size.impact.correlate.indiv,
        proptx.impact.correlate = FALSE,
        correlate.strength = correlate.strength.indiv,
        size.ratio = size.ratio.indiv
      ))) %>%
    tidyr::unnest(cols=sim_data) %>% 
    select(-J, -p, -district.tau)
  
  # ensure that sids are unique
  sid_key <- sim.data %>%
    select(did, sid) %>%
    unique() %>%
    mutate(sid_new = 1:n())
  
  sim.data <- sim.data %>%
    left_join(sid_key, by=c("did", "sid")) %>%
    mutate(sid = sid_new) %>%
    select(-sid_new)
  
  # record J, p, district.tau for each district
  attr(sim.data, "district_key") <- grid
  
  # assign cluster-randomized treatment
  sim.data <- rerandomize.data.full(sim.data)
  
  sim.data
}

if (F) {
  df <- generate_districts(
    K = 5, 
    J.bar = 20,
    ATE.district = 0.2, 
    sd.district = 0.1,
    variable.J = T,

    # cluster-level params
    n.bar = 50, 
    tau.11 = 0.1, 
    ICC = 0.20,
    variable.n = T,
    
    # size and correlation params
    size.ratio.clust = 0.6, 
    size.ratio.indiv = 0.6,
    size.impact.correlate.clust = 1,
    size.impact.correlate.indiv = 1,
    proptx.impact.correlate.clust = F,
    correlate.strength.clust = 0.5,
    correlate.strength.indiv = 0.5)
  
  describe.data(df)
  
}


##### Simulation code #####

# Run a single round of simulation.
#
# Generate a cluster-randomized trial, and then repeatedly rerandomize the fixed 
# dataset n.runs times, each time estimating the  ATE using all our methods.  
# It also calculates the true finite-sample ATE (person and site weighted).
#
# @return Dataframe with the estimates along with the true baseline values.
single.CRT = function( K, J.bar,
                       ATE.district, sd.district,
                       variable.J = FALSE,

                       # cluster-level params
                       n.bar, tau.11, ICC = 0.20,
                       variable.n = FALSE,
                       
                       # size and correlation params
                       size.ratio.clust = 0.6, 
                       size.ratio.indiv = 0.6,
                       size.impact.correlate.clust = 0,
                       size.impact.correlate.indiv = 0,
                       proptx.impact.correlate.clust = F,
                       correlate.strength.clust = 0.75,
                       correlate.strength.indiv = 0.75,
                       
                       # method params
                       n.runs = 3,
                       just.describe.data = FALSE, 
                       just.return.data = FALSE,
                       
                       ... ) {
  
  df <- generate_districts(
    K = K, 
    J.bar = J.bar,
    ATE.district = ATE.district, 
    sd.district = sd.district,
    variable.J = variable.J,
    
    # cluster-level params
    n.bar = n.bar, 
    tau.11 = tau.11, 
    ICC = ICC,
    variable.n = variable.n,
    
    # size and correlation params
    size.ratio.clust = size.ratio.clust, 
    size.ratio.indiv = size.ratio.indiv,
    size.impact.correlate.clust = size.impact.correlate.clust,
    size.impact.correlate.indiv = size.impact.correlate.indiv,
    proptx.impact.correlate.clust = proptx.impact.correlate.clust,
    correlate.strength.clust = correlate.strength.clust,
    correlate.strength.indiv = correlate.strength.indiv)
  
  if ( just.return.data ) {
    return( df )
  }
  
  if ( just.describe.data ) {
    return( describe.data( df ) )
  }
  
  # Calculate true finite sample estimands
  params = df %>% group_by( sid ) %>%
    summarise( ATE = mean( Y1 ) - mean( Y0 ),
               n = n(), .groups="drop" )
  ATE.person = with( df, mean( Y1 - Y0 ) )
  ATE.site = mean( params$ATE )
  
  # Do finite-sample simulation, rerandomizing n.runs times.  
  ests = plyr::rdply( n.runs, {
    df = rerandomize.data( df )
    get.estimates( df, ... ) 
  }, .id="subrun" )
  
  # Gather results and ship!
  ests = mutate( ests,
                 ATE.finite.person=ATE.person, 
                 ATE.finite.site=ATE.site )
  ests    
}



# Test the simulation code
if ( FALSE ) {
  
  rst = single.MLM.trial( 20, 10, 0.2^2, variable.n = TRUE, variable.p = TRUE,
                          dependence = TRUE, proptx.dependence = TRUE )
  head( rst )
  table( rst$subrun )
  
}


# Run a simulation with R trials and return all the simulation runs as a large dataframe.
run.scenario = function( J, n.bar, tau, dependence, proptx.dependence, variable.n, variable.p, 
                         ATE = 0.20, ICC = 0.20, p.tx = 0.65, R = 10, n.runs=3, 
                         ID=NULL, .progress="none", ... ) {
  ptm = proc.time()
  
  #  n = n.bar * J
  #  e.time = 1.13^(J/5) + 1.0183^(n/100)
  
  #  R.adj = round( R * pmin( 1, 3/e.time ) )
  #if ( J > 40 ) {
  #  R.adj = R / 2
  #} else {
  #  R.adj = R
  #}
  #R.adj = R
  
  if ( is.null( ID ) ) {
    ID = -999
  }
  
  scat( "Running scenario %d\n\tJ=%d\tn.bar=%d\ttau=%.2f\tATE=%.2f\tICC=%.2f\tprop tx=%.2f\tdependence=%s\tprop dep=%s\tR=%d\n", 
        ID, J, n.bar,tau, ATE, ICC, p.tx, dependence, proptx.dependence, R)
  
  rps = plyr::rdply( R,  single.MLM.trial( n.bar=n.bar, J=J, tau.11.star=tau, 
                                           dependence = dependence, proptx.dependence = proptx.dependence,
                                           variable.n = variable.n, variable.p = variable.p,
                                           ATE.superpop=ATE, ICC=ICC, p.tx=p.tx, n.runs=n.runs, ... ),
                     .id="run", .progress="text" )
  
  rps$subrun = paste0( rps$run, ".", rps$subrun )
  
  scat("**\n**\tTotal time elapsed:\n")
  tot.time = proc.time() - ptm
  print(tot.time)
  sim.per.min = R / (tot.time["elapsed"] / 60)
  scat("Simulations per minute = %.2f\n", sim.per.min )
  rps$sim.per.min = sim.per.min
  rps
}


# Testing run.scenario()
if ( FALSE ) {
  rr = run.scenario( J = 4, n.bar = 16, tau = 0.3^2, dependence=TRUE, 
                     proptx.dependence = TRUE, 
                     variable.n = TRUE, variable.p = TRUE,
                     ATE = 0.2, R = 10 )
  head( rr )  
  table( rr$run )
  table(  table( rr$run ) )
  table(  table( rr$subrun ) )
  tt = sort( table( rr$subrun ) )
  tt
  
  filter( rr, subrun %in% c( "1.2", "4.1", "4.2", "2.1" ) )
  
}



#### Make scenario configurations #####

make.scenario.list = function( group = "main" ) {
  
  if ( group == "main" ) {
    scenarios = expand.grid( J = c( 80, 40, 20, 10 ),
                             n.bar = c( 8000, 4000, 2000, 1000 ), # put in totals here
                             dependence = c( 1, 0 ),
                             proptx.dependence = c( 1, 0, -1 ), 
                             variable.n = c( TRUE, FALSE ),
                             variable.p = c( TRUE, FALSE ),
                             ATE = 0.20,
                             tau = c( 0, 0.1, 0.20 )^2,
                             p.tx = 0.65,
                             ICC = 0.20 )
  } else if ( simstudy == "initial" ) {
    scenarios = expand.grid( J = c( 80, 40, 20 ),
                             n.bar = c( 8000, 4000, 2000 ), # put in totals here
                             dependence = c( 1, 0 ),
                             proptx.dependence = c( 1, 0, -1 ), 
                             variable.n = c( TRUE, FALSE ),
                             variable.p = c( TRUE, FALSE ),
                             ATE = 0.20,
                             tau = c( 0, 0.1, 0.20 )^2,
                             p.tx = 0.70,
                             ICC = 0.60 )
  } else {
    error( "No defined group in make.scenario.list" )
  }
  
  nrow( scenarios )
  
  # drop redundant scenarios (dependence irrelevant if quantity not varying)
  scenarios = filter( scenarios, (variable.n == TRUE) | (dependence == 0) )
  scenarios = filter( scenarios, (variable.p == TRUE) | (proptx.dependence == 0) )
  nrow( scenarios )
  
  # drop redundant/not well defined scenarios (dependence irrelevant if no cross site variation )
  scenarios = filter( scenarios, !( (tau == 0) & (dependence != 0) ) )
  scenarios = filter( scenarios, !( (tau == 0) & (proptx.dependence != 0) ) )
  nrow( scenarios )
  
  table( scenarios$variable.n, scenarios$dependence )
  table( scenarios$variable.p, scenarios$proptx.dependence )
  
  table( scenarios$variable.n, scenarios$tau )
  
  # convert n.bar from totals to people per site
  scenarios = mutate( scenarios, n.bar = round( n.bar / J ) )
  scenarios = as_tibble( scenarios )
  scenarios
  nrow( scenarios )
  
  scenarios$ID = 1:nrow( scenarios )
  
  scenarios  
}


# Given an Id 'index' return the charactaristics of the given generated scenario.
# 
# This is so we don't have to store all the factors of the experiment in all the dataframes.
get.scenario.by.id = function( index, group = "main" ) {
  scenarios = make.scenario.list( group=group )
  scenarios[ index, ]
}


# Testing code
if ( FALSE ) {
  sc = make.scenario.list()
  nrow( sc )
  head( sc )
  table( var_n = sc$variable.n, var_p = sc$variable.p )
  
  get.scenario.by.id(12)
}


