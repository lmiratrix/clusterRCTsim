

# This script has all the support code to run large simulations to
# look at stability of the different SE estimators and true
# performances of the ATE estimators.
#
# This script has functions to run a single trial, generate a list of
# all scenarios, and a few other things.

library( tidyverse )
library( purrr )
library( multiwayvcov )

library( clusterRCT )


#### Utility functions #####

if ( FALSE ) {
  # for testing crashing in the simulator.  This wrapper will make things crash sometimes.
  my_compare_methods = function( ... ) {
    if ( sample(5,1) == 1 ) {
      asdfadfga
    } else {
      compare_methods( ... )
    }
  }
}


# Print things to console using sprintf.
scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}



# Given dataframe df, estimate impacts.
# 
# df has to be in canonical form of Yobs ~ Z | sid | did.
#
# @return Dataframe with estimates of ATE and SE, one row per method, 
# from compare_methods() function
get.estimates = function( df, ... ) {
  safe_comp = safely( clusterRCT::compare_methods )
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


# Reassign treatment to all clusters within blocks.
rerandomize.data = function( dat ) {
  
  grid <- attr( dat, "district_key" )
  
  if ( "did" %in% names(dat) ) {
    dat <- dat %>% group_by( did ) %>% 
      mutate( p = mean( Z[ !duplicated( sid ) ] ),
              Z = randomizr::block_and_cluster_ra(
                blocks=did,
                clusters=sid, 
                prob_unit=p),
              Yobs = ifelse(Z, Y1, Y0)) %>% 
      dplyr::select(-p)
  } else {
    dat <- dat %>% 
      mutate( p = mean( Z[ !duplicated( sid ) ] ),
              Z = randomizr::cluster_ra(
                clusters=sid, 
                prob_unit=p),
              Yobs = ifelse(Z, Y1, Y0)) %>%
      dplyr::select(-p)
  }
  
  attr( dat, "district_key" ) <- grid
  
  return( dat )
}


#### Code to explore scenarios and sample datasets ####

# Utility function to describe a given data set to explore the data
# generation process
#
# TODO
describe.data = function( df ) {
  
  # we want:
  #  - true ATEs
  #  - true correlations from district_key
  #  - etc...
  
  # Why this needed?
  df = rerandomize.data( df )
  
  
  # district params
  d_key = attr(df, "district_key")
  if ( !is.null( d_key ) ) {
    sstat_distr <-d_key %>% 
      summarize(cor.J.tau_dist = cor(J, district.tau),
                cor.p.tau_dist = cor(p, district.tau),
                sd.tau_dist = sd(district.tau),
                sd.J_dist = sd(J),
                sd.p_dist = sd(p))
    print(sstat_distr)
    
  } else {
    df$did = "single"
  }
  
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
  
  print(res)
}

# use describe_clusterRCT()


##### Testing code #####

# Testing code of blkvar's generate_multilevel_data_no_cov
# and describe_data method
if ( FALSE ) {
  library( tidyverse )
  
  df = blkvar::generate_multilevel_data_no_cov( n.bar=200, J=50,
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
  df = rerandomize.data(df)
  head( df )
  
  describe.data(df)
  describe_data(df)  
  
  
  get.estimates(df)
}



#### Simulation code ####

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
  ICC.district,
  variable.J = FALSE,
  variable.p = FALSE,
  
  # correlation params at cluster level
  size.ratio.district = 0.6, 
  size.impact.cor.district = 0,
  proptx.impact.cor.district = 0,
  size.intercept.cor.district = 0,
  
  # school-level params
  n.bar, sd.school, ICC.school = 0.20,
  variable.n = FALSE,
  
  # correlation params at individual level
  size.ratio.school = 0.6,
  size.impact.cor.school = 0,
  cor.strength.school = 0.75,
  
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
    ICC.district = ICC.district,
    variable.J = variable.J,
    variable.p = variable.p,
    
    # cluster-level params
    n.bar = n.bar, 
    sd.school = sd.school, 
    ICC.school = ICC.school,
    variable.n = variable.n,
    
    # size and correlation params
    size.ratio.district = size.ratio.district, 
    size.ratio.school = size.ratio.school,
    size.impact.cor.district = size.impact.cor.district,
    size.intercept.cor.district = size.intercept.cor.district,
    size.impact.cor.school = size.impact.cor.school,
    proptx.impact.cor.district = proptx.impact.cor.district,
    cor.strength.school = cor.strength.school)
  
  if ( just.return.data ) {
    return( df )
  }
  
  if ( just.describe.data ) {
    return( describe.data( df ) )
  }
  
  # Calculate true finite-sample estimands
  beta_jks <- df %>% 
    group_by(did, sid) %>% 
    summarize(n_jk = n(), 
              beta_jk = mean(Y1-Y0), .groups="drop")
  
  ATE.student <- with( df, mean( Y1 - Y0 ) )
  # ATE.student = beta_jks %>% 
  #   summarize(n_k = sum(n_jk),
  #             s = sum(n_jk*beta_jk)) %>% 
  #   summarize(ATE.student = sum(s) / sum(n_k))
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
  }, .id="subrun" )
  
  # Gather results and ship!
  ests = mutate( ests,
                 ATE.finite.student = ATE.student,
                 ATE.finite.school = ATE.school,
                 ATE.finite.dist.student = ATE.dist.student,
                 ATE.finite.dist.school = ATE.dist.school)
  
  if (F) {
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



# Test the simulation code
if ( FALSE ) {
  
  # NOTE: DB estimators will return infinite SEs if there are
  #  fewer than two tx and co clusters in any district
  res <- single.CRT(
    K=10, J.bar=20, p.bar=0.5,
    variable.J=F, variable.p=F,
    ATE.district=0.2, sd.district=0.1, 
    ICC.district=0.2,
    
    size.impact.cor.district = 0,
    proptx.impact.cor.district = -0,
    
    n.bar=50, sd.school=0.1, ICC.school=0.3,
    variable.n=TRUE,
    include_method_characteristics=F)
  
  head(res)
}


# Run a simulation with R trials and return all the simulation runs as a large dataframe.
run.scenario = function( 
    K, J.bar, p.bar,
    ATE.district, sd.district,
    ICC.district,
    variable.J,
    variable.p,
    
    # correlation params at cluster level
    size.ratio.district = 0.6,
    size.impact.cor.district,
    size.intercept.cor.district = 0,
    proptx.impact.cor.district,
    
    # school-level params
    n.bar, sd.school, ICC.school = 0.2,
    variable.n,
    
    # correlation params at individual level
    size.ratio.school = 0.6,
    size.impact.cor.school,
    cor.strength.school = 0.75,
    
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
  
  scat( "Running scenario %d\n\tK=%d\tJ.bar=%d\tp.bar=%.2f\t
        district ATE=%.2f\tdistrict sd=%.2fdistrict ICC=%.2f\t
        district size correlate=%.2f\tdistrict proptx correlate=%.2f
        n.bar=%d\tschool sd=%.2f\tICC=%.2f
        school size correlate=%s
        R=%d\n", 
        ID, K, J.bar, p.bar,
        ATE.district, sd.district, 
        ICC.district,
        size.impact.cor.district, proptx.impact.cor.district,
        n.bar, sd.school, ICC.school,
        size.impact.cor.school,
        R)
  
  # Get all the ... arguments in list
  dots = list(...)
  scat( "\tExtra arguments: %s\n", paste(names(dots), collapse=", ") )
  
  rps = map_dfr( 1:R, \(i) single.CRT( 
    K = K, 
    J.bar = J.bar,
    p.bar = p.bar,
    ATE.district = ATE.district, 
    sd.district = sd.district,
    ICC.district = ICC.district,
    variable.J = variable.J,
    variable.p = variable.p,
    
    # cluster-level params
    n.bar = n.bar, 
    sd.school = sd.school, 
    ICC.school = ICC.school,
    variable.n = variable.n,
    
    # size and correlation params
    size.ratio.district = size.ratio.district, 
    size.ratio.school = size.ratio.school,
    size.impact.cor.district = size.impact.cor.district,
    size.intercept.cor.district = size.intercept.cor.district,
    size.impact.cor.school = size.impact.cor.school,
    proptx.impact.cor.district = proptx.impact.cor.district,
    cor.strength.school = cor.strength.school,
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
                     include_method_characteristics=F)
  head( rr )  
}



#### Make scenario configurations #####
#base.values
#function to replicate base valuse across the rows and only change one variable
expand.fixed.one<-function(manip.var.name,manip.var.vals, base.values){
  base.comp<-base.values %>% select(!sim.id)
  scenarios<- base.comp %>% slice(rep(1:n(),each=length(manip.var.vals))) %>%
    select(-all_of(manip.var.name)) %>% mutate(!!manip.var.name:=manip.var.vals) %>% 
    filter(!pmap_lgl(.,~all(c(...)==base.comp))) %>% mutate (sim.id=manip.var.name)}

#set base.values
base.values<-data.frame(K=10, 
                        J.bar = 10,
                        p.bar = 0.5,
                        ATE.district = 0.2, 
                        sd.district = 0.1,
                        ICC.district = 0.2,
                        variable.J = T,
                        variable.p = F,
                        n.bar=25,
                        sd.school = 0.1, 
                        ICC.school = 0.2,
                        variable.n = T,
                        size.impact.cor.district = 0,
                        size.impact.cor.school = 0,
                        proptx.impact.cor.district = 0,
                        size.intercept.cor.district = 0,
                        cor.strength.school=0,  
                        sim.id="base.values")

base.values.variable.p<- base.values %>% mutate(variable.p = T,sim.id="base.values.variable.p")
base.values.sch.cor<- base.values %>% mutate( size.impact.cor.school = 1,sim.id="base.values.sch.cor")
base.values.sch.cor.neg<- base.values %>% mutate( size.impact.cor.school = -1,sim.id="base.values.sch.cor.neg")  
base.values.estimand<- base.values %>% mutate( size.impact.cor.school = 1,cor.strength.school=.75,sim.id="base.values.estimand")
base.values.estimand.variable.p<- base.values.estimand %>% mutate(variable.p = T,sim.id="base.values.estimand.variable.p")
base.values.estimand.2<-base.values.estimand %>% mutate(sd.district=.3,sd.school=.3,sim.id="base.values.estimand.2")
base.values.estimand.variable.p.2<- base.values.estimand.2 %>% mutate(variable.p = T,sim.id="base.values.estimand.variable.p.2")
base.values.performance<- base.values %>% mutate(size.impact.cor.district = .2, 
                                                 size.impact.cor.school =1, cor.strength.school=.2,sim.id="base.values.performance")
base.values.performance.variable.p<- base.values.performance %>% mutate(variable.p = T,sim.id="base.values.performance.variable.p")
base.values.performance.sch.cor<- base.values.performance %>% mutate( size.impact.cor.school = 1,sim.id="base.values.performance.sch.cor")
base.values.performance.sch.cor.neg<- base.values.performance %>% mutate( size.impact.cor.school = -1,sim.id="base.values.performance.sch.cor.neg")  




make.scenario.list = function( group = "main" ) {
  
  if ( group == "main" ) {
    scenarios = expand.grid( 
      K = c(5,10,20), 
      J.bar = c(10,20,40,80),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = 0.2,
      variable.J = c(T,F),
      variable.p = c(T,F),
      
      # cluster-level params
      n.bar = c(25,50,100), 
      sd.school = c(0.1, 0.2), 
      ICC.school = 0.2,
      variable.n = c(T,F),
      
      # size and correlation params
      size.impact.cor.district = c(0,0.75),
      size.impact.cor.school = c(0,1),
      proptx.impact.cor.district = c(0, 0.75))
  } else if ( group == "sim_base_1" ) {
    scenarios<-bind_rows(expand.fixed.one(manip.var.name="size.impact.cor.district",
                                          manip.var.vals<-seq(-1,1,.1),base.values),
                         expand.fixed.one(manip.var.name="size.impact.cor.school",
                                          manip.var.vals<-c(-1,0,1),base.values),
                         expand.fixed.one(manip.var.name="proptx.impact.cor.district",
                                          manip.var.vals<-seq(-1,1,.1),base.values.variable.p),
                         expand.fixed.one(manip.var.name="size.intercept.cor.district",
                                          manip.var.vals<-seq(-1,1,.1),base.values),
                         expand.fixed.one(manip.var.name="ICC.district",
                                          manip.var.vals<-seq(.05,.7,.05),base.values),
                         expand.fixed.one(manip.var.name="ICC.school",
                                          manip.var.vals<-seq(.05,.7,.05),base.values), 
                         expand.fixed.one(manip.var.name="sd.district",
                                          manip.var.vals<-seq(0,.5,.02),base.values), 
                         expand.fixed.one(manip.var.name="sd.school",
                                          manip.var.vals<-seq(0,.5,.02),base.values),
                         expand.fixed.one(manip.var.name="cor.strength.school",
                                          manip.var.vals<-seq(0,1,.1),base.values.sch.cor),
                         expand.fixed.one(manip.var.name="cor.strength.school",
                                          manip.var.vals<-seq(0,1,.1),base.values.sch.cor.neg),
                         base.values,base.values.variable.p,base.values.sch.cor,base.values.sch.cor.neg)
    
  } else if ( group == "initial" ) {
    scenarios = expand.grid( 
      K = c(10), 
      J.bar = c(20,50),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = 0.2,
      variable.J = c(T,F),
      variable.p = c(T,F),
      
      # cluster-level params
      n.bar = c(50), 
      sd.school = c(0.1), 
      ICC.school = 0.2,
      variable.n = c(T,F),
      
      # size and correlation params
      size.impact.cor.district = c(0,0.75),
      size.impact.cor.school = c(0,1),
      proptx.impact.cor.district = c(0,0.75))
    
  } else if ( group == "sim2" ) {
    scenarios = expand.grid( 
      K = c(10), 
      J.bar = c(25),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = 0.2,
      variable.J = c(T),
      variable.p = c(F),
      
      # cluster-level params
      n.bar = c(25,100), 
      sd.school = c(0.1), 
      ICC.school = 0.2,
      variable.n = c(T),
      
      # size and correlation params
      size.impact.cor.district = c(-0.75,0,0.75),
      size.impact.cor.school = c(-1,0,1),
      proptx.impact.cor.district = c(0))
    
  } else if ( group == "sim3" ) {
    scenarios = expand.grid( 
      K = c(5,10,20), 
      J.bar = c(10,25,80),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = 0.2,
      variable.J = c(T),
      variable.p = c(T),
      
      # cluster-level params
      n.bar = c(25,50,100), 
      sd.school = c(0.1), 
      ICC.school = 0.2,
      variable.n = c(T),
      
      # size and correlation params
      size.impact.cor.district = c(0,0.75),
      size.impact.cor.school = c(0,1),
      proptx.impact.cor.district = c(0))
    
  } else if ( group == "sim_icc" ) {
    scenarios = expand.grid( 
      K = c(10), 
      J.bar = c(20),
      p.bar = c(0.5),
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = c(0.1, 0.25, 0.4),
      variable.J = c(T),
      variable.p = c(T),
      
      # cluster-level params
      n.bar = c(50), 
      sd.school = c(0.1), 
      ICC.school = c(0.1, 0.25, 0.4),
      variable.n = c(T),
      
      # size and correlation params
      size.impact.cor.district = c(0,0.75),
      size.impact.cor.school = c(0,1),
      proptx.impact.cor.district = c(0))
    
  } else if ( group == "sim_interceptcor" ) {
    scenarios = expand.grid( 
      K = c(15), 
      J.bar = c(15),
      p.bar = c(0.75),       # higher pbar
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = c(0.2),
      variable.J = c(T),
      variable.p = c(T),
      
      # cluster-level params
      n.bar = c(50), 
      sd.school = c(0.1), 
      ICC.school = c(0.2),
      variable.n = c(T),
      
      # size and correlation params
      size.impact.cor.district = c(0, 0.8),
      size.impact.cor.school = c(0, 1),
      proptx.impact.cor.district = c(0, 0.8),
      size.intercept.cor.district = c(0, 0.8))
    
  } else if ( group == "sim_small" ) {
    scenarios = expand.grid( 
      K = c(4), 
      J.bar = c(6),
      p.bar = c(0.75),
      ATE.district = 0.2, 
      sd.district = 0.1,
      ICC.district = c(0.2),
      variable.J = c(T),
      variable.p = c(T),
      
      # cluster-level params
      n.bar = c(50), 
      sd.school = c(0.1), 
      ICC.school = c(0.2),
      variable.n = c(T),
      
      # size and correlation params
      size.impact.cor.district = c(0, 0.75),
      size.impact.cor.school = c(0, 1),
      proptx.impact.cor.district = c(0, 0.75),
      size.intercept.cor.district = c(0, 0.75))
    
  } else {
    error( "No defined group in make.scenario.list" )
  }
  
  nrow( scenarios )
  
  # filter out replicated correlation rows (if no room for variation)
  scenarios <- scenarios %>% 
    filter(
      variable.J | size.impact.cor.district==0,
      variable.J | size.intercept.cor.district==0,
      variable.p | proptx.impact.cor.district==0,
      sd.district!=0 | (size.impact.cor.district == 0 & proptx.impact.cor.district == 0),
      variable.n | size.impact.cor.school == 0,
      sd.school != 0 | size.impact.cor.school == 0
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


##### Testing code #####

# Testing code
if ( FALSE ) {
  sc = make.scenario.list(group="main")
  nrow( sc )
  head( sc )
  
  get.scenario.by.id(12, group="main")
}

