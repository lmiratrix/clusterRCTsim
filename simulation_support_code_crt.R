
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
  ests = safe_comp( Yobs ~ Z | sid | did, data=df, ... )
  if ( is.null( ests[[2]] ) ) {
    ests = ests[[1]]
    ests = rename( ests, 
                   ATE.hat = ATE_hat,
                   SE.hat = SE_hat )
    
  } else {
    cat( "\nCaught error:\n" )
    print( ests[[2]] )
    ests = data.frame( ATE.hat = c(NA),
                       SE.hat = c(NA) )    
  }
  ests
}



# Reassign treatment within a district (to do finite sample inference)
rerandomize.data.district = function( dat, p ) {
  stopifnot(length(unique(dat$did)) == 1)
  dat %>%
    mutate(
      # sid = droplevels(sid),   # avoid empty factor level issues with clusters=sid
      Z = randomizr::cluster_ra(clusters = sid, prob = p),
      Yobs = ifelse( Z == 1, Y1, Y0 ) )
}
# Reassign treatment to all districts
rerandomize.data = function( dat ) {
  grid <- attr(dat, "district_key")
  
  dat <- dat %>% 
    split(~did) %>% 
    purrr::map2(grid$p, ~rerandomize.data.district(.x, .y)) %>% 
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
  
  df = rerandomize.data( df )
  
  
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
  stopifnot(J.bar > 5)
  nk <- blkvar::block_distn(J=K, n.bar=J.bar-4, size.ratio)
  nk[nk<4] <- 4
  
  nk
}

# Generate a cluster-randomized trial, with multiple districts
# 
# District average effects may be correlated with
#  district size (i.e., # clusters) and proportion of clusters treated
# School average effects may have all correlations
#  encoded in blkvar package
generate_districts <- function(
    # district-level params
    K, J.bar, p.bar=0.5,
    ATE.district, sd.district,
    variable.J = FALSE,
    variable.p = FALSE,
    
    # correlation params at cluster level
    size.ratio.clust = 0.6, 
    size.impact.correlate.clust = 0,
    proptx.impact.correlate.clust = 0,
    
    # school-level params
    n.bar, tau.11, ICC = 0.20,
    variable.n = FALSE,
    
    # correlation params at individual level
    size.ratio.indiv = 0.6,
    size.impact.correlate.indiv = 0,
    correlate.strength.indiv = 0.75) {
  
  stopifnot("size.impact.correlate.clust must be in (-1,1)" = 
              size.impact.correlate.clust >= -1 & size.impact.correlate.clust <= 1)
  stopifnot("proptx.impact.correlate.clust must be in (-1,1)" = 
              proptx.impact.correlate.clust >= -1 & proptx.impact.correlate.clust <= 1)
  stopifnot(size.impact.correlate.indiv %in% c(-1,0,1))
  
  # district J values
  if (variable.J) {
    num_clusters <- cluster_distn(K, J.bar, size.ratio=size.ratio.clust)
  } else {
    num_clusters <- rep(J.bar, K)
  }
  # district proportions treated
  if (variable.p) {
    pp <- min(p.bar, 1-p.bar)
    prop_tx_clusters <- runif(K, min=p.bar-pp/2, max=p.bar+pp/2)
  } else {
    prop_tx_clusters <- rep(p.bar, K)
  }
  
  # district-level parameters
  grid <- tibble(
    did = 1:K,
    district.tau = rnorm(K, mean=ATE.district, sd=sd.district))
  
  if (size.impact.correlate.clust != 0) {
    stopifnot("Cannot correlate district sizes and average impacts if sizes do not vary" = variable.J)
    
    # fuzzy sort num_clusters
    num_clusters_temp <- num_clusters[order(
      size.impact.correlate.clust^2 * num_clusters + 
        (1-size.impact.correlate.clust^2) * rnorm(K, sd=sd(num_clusters)) )]
    num_clusters <- num_clusters_temp[ order(order(sign(size.impact.correlate.clust)*grid$district.tau)) ]
  }
  if (proptx.impact.correlate.clust != 0) {
    stopifnot("Cannot correlate district proportions treated and average impacts if proportions treated do not vary" = variable.p)
    # fuzzy sort prop_tx_clusters
    p_temp <- prop_tx_clusters[order(
      proptx.impact.correlate.clust^2 * prop_tx_clusters + 
        (1-proptx.impact.correlate.clust^2) * rnorm(K, sd=sd(prop_tx_clusters)) )]
    prop_tx_clusters <- p_temp[ order(order(sign(proptx.impact.correlate.clust)*grid$district.tau)) ]
  }
  
  grid$J <- num_clusters
  grid$p <- prop_tx_clusters
  
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
  sim.data <- rerandomize.data(sim.data)
  
  sim.data
}

if (F) {
  df <- generate_districts(
    K = 100,
    J.bar = 20, 
    p.bar = 0.5,
    ATE.district = 0.2, 
    sd.district = 0.1,
    variable.J = T,
    variable.p = T,
    
    # correlation params at cluster level
    size.ratio.clust = 0.6, 
    size.impact.correlate.clust = 0,
    proptx.impact.correlate.clust = 0,
    
    # school-level params
    n.bar = 10, 
    tau.11 = 0, 
    ICC = 0.20,
    variable.n = T,
    
    # correlation params at individual level
    size.ratio.indiv = 0.6,
    size.impact.correlate.indiv = 0,
    correlate.strength.indiv = 0.75
  )
  
  df %>% 
    group_by(did, sid) %>% 
    summarize(n = n(),
              Z = mean(Z),
              tauhat = mean(Y1-Y0)) %>% 
    summarize(J = n(),
              p = mean(Z),
              tauhat = mean(tauhat)) %>% 
    summarize(corJ = cor(tauhat, J, method="spearman"),
              corp = cor(tauhat, p, method="spearman"))
  
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
single.CRT = function( # district-level params
  K, J.bar, p.bar=0.5,
  ATE.district, sd.district,
  variable.J = FALSE,
  variable.p = FALSE,
  
  # correlation params at cluster level
  size.ratio.clust = 0.6, 
  size.impact.correlate.clust = 0,
  proptx.impact.correlate.clust = 0,
  
  # school-level params
  n.bar, tau.11, ICC = 0.20,
  variable.n = FALSE,
  
  # correlation params at individual level
  size.ratio.indiv = 0.6,
  size.impact.correlate.indiv = 0,
  correlate.strength.indiv = 0.75,
  
  # method params
  n.runs = 3,
  just.describe.data = FALSE, 
  just.return.data = FALSE,
  ... ) {
  
  df <- generate_districts(
    K = K, 
    J.bar = J.bar,
    p.bar = p.bar,
    ATE.district = ATE.district, 
    sd.district = sd.district,
    variable.J = variable.J,
    variable.p = variable.p,
    
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
  
  # TODO: are we still interested in the same finite stuffs?
  
  # Gather results and ship!
  ests = mutate( ests,
                 ATE.finite.person=ATE.person, 
                 ATE.finite.site=ATE.site )
  ests    
}



# Test the simulation code
if ( F ) {
  
  res <- single.CRT(
    K=10, J.bar=10, p.bar=0.5,
    variable.J=T, variable.p=T,
    ATE.district=0.2, sd.district=0.1, 
    size.impact.correlate.clust = 0.5,
    proptx.impact.correlate.clust = -0.5,
    n.bar=50, tau.11=0.1^2,
    variable.n=T,
    include_method_characteristics=F)
  
  res
}


# Run a simulation with R trials and return all the simulation runs as a large dataframe.
run.scenario = function( 
    K, J.bar, p.bar,
    ATE.district, sd.district,
    variable.J,
    variable.p,
    
    # correlation params at cluster level
    size.ratio.clust = 0.6,
    size.impact.correlate.clust,
    proptx.impact.correlate.clust,
    
    # school-level params
    n.bar, tau.11, ICC = 0.20,
    variable.n,
    
    # correlation params at individual level
    size.ratio.indiv = 0.6,
    size.impact.correlate.indiv,
    correlate.strength.indiv = 0.75,
    
    # method params
    n.runs = 3,
    R = 10,

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
  
  scat( "Running scenario %d\n\tK=%d\tJ.bar=%d\tp.bar=%.2f\tdistrict ATE=%.2f\tdistrict sd=%.2f
        district size correlate=%.2f\tdistrict proptx correlate=%.2f
        n.bar=%d\tschool sd=%.2f\tICC=%.2f
        school size correlate=%s
        R=%d\n", 
        ID, K, J.bar, p.bar,
        ATE.district, sd.district, 
        size.impact.correlate.clust, proptx.impact.correlate.clust,
        n.bar, tau.11, ICC,
        size.impact.correlate.indiv,
        R)
  
  rps = plyr::rdply( R, single.CRT( 
    K = K, 
    J.bar = J.bar,
    p.bar = p.bar,
    ATE.district = ATE.district, 
    sd.district = sd.district,
    variable.J = variable.J,
    variable.p = variable.p,
    
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
    correlate.strength.indiv = correlate.strength.indiv,
    n.runs=n.runs, ... ),
    .id="run", .progress=.progress )
  
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
  rr = run.scenario( K = 10, 
                     J.bar = 10,
                     p.bar = 0.5,
                     ATE.district = 0.2, 
                     sd.district = 0.1,
                     variable.J = T,
                     variable.p = T,
                     
                     # cluster-level params
                     n.bar = 50, 
                     tau.11 = 0.1, 
                     ICC = 0.2,
                     variable.n = T,
                     
                     # size and correlation params
                     size.impact.correlate.clust = 0,
                     size.impact.correlate.indiv = 0,
                     proptx.impact.correlate.clust = 0,

                     R = 2,
                     .progress="text",
                     include_method_characteristics=F)
  head( rr )  
}



#### Make scenario configurations #####

make.scenario.list = function( group = "main" ) {
  
  if ( group == "main" ) {
    scenarios = expand.grid( 
      K = c(5,10,20), 
      J.bar = c(10,20,40,80),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      variable.J = c(T,F),
      variable.p = c(T,F),
      
      # cluster-level params
      n.bar = c(25,50,100), 
      tau.11 = c(0.1, 0.2), 
      ICC = 0.2,
      variable.n = c(T,F),
      
      # size and correlation params
      size.impact.correlate.clust = c(-0.75,0,0.75),
      size.impact.correlate.indiv = c(-0.75,0,0.75),
      proptx.impact.correlate.clust = c(-1,0,1))
  } else if ( group == "initial" ) {
    scenarios = expand.grid( 
      K = c(10), 
      J.bar = c(20,80),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      variable.J = c(T,F),
      variable.p = c(T,F),
      
      # cluster-level params
      n.bar = c(25,100), 
      tau.11 = c(0.1), 
      ICC = 0.2,
      variable.n = c(T,F),
      
      # size and correlation params
      size.impact.correlate.clust = c(0,0.75),
      size.impact.correlate.indiv = c(0,0.75),
      proptx.impact.correlate.clust = c(0,1))
    
  } else {
    error( "No defined group in make.scenario.list" )
  }
  
  nrow( scenarios )
  
  scenarios <- scenarios %>% 
    filter(
      variable.J | size.impact.correlate.clust==0,
      variable.p | proptx.impact.correlate.clust==0,
      sd.district!=0 | (size.impact.correlate.clust == 0 & proptx.impact.correlate.clust == 0),
      variable.n | size.impact.correlate.indiv == 0,
      tau.11 != 0 | size.impact.correlate.indiv == 0
    )
  
  # convert n.bar from totals to people per site
  
  scenarios$ID = 1:nrow( scenarios )
  
  scenarios  
}


# Given an Id 'index' return the characteristics of the given generated scenario.
# 
# This is so we don't have to store all the factors of the experiment in all the dataframes.
get.scenario.by.id = function( index, group = "main" ) {
  scenarios = make.scenario.list( group=group )
  scenarios[ index, ]
}


# Testing code
if ( FALSE ) {
  sc = make.scenario.list(group="initial")
  nrow( sc )
  head( sc )
  table( var_n = sc$variable.n, var_p = sc$variable.p )
  
  get.scenario.by.id(12, group="initial")
}


