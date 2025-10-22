

# This script has all the support code to run large simulations to
# look at stability of the different SE estimators and true
# performances of the ATE estimators.
#
# This script has functions to run a single trial, generate a list of
# all scenarios, and a few other things.
#
# NOTE: This regenerates blocks (blocks) and so operates in a very
# super-population framework.  We now usually fix the blocks, and
# use the companion "fixed_block" version of this script.

library( tidyverse )
library( purrr )
library( multiwayvcov )

library( clusterRCT )

source( here::here( "R/util.R" ) )


#### Simulation code ####


# Run a single round of simulation.
#
# Generate a cluster-randomized trial, and then repeatedly rerandomize the fixed 
# dataset n.runs times, each time estimating the ATE using all our methods.  
# It also calculates the true finite-sample ATE (person and site weighted).
#
# @return Dataframe with the estimates along with the true baseline values.
single.CRT = function( # block-level params
  K, J.bar, p.bar=0.5,
  ATE.block, sd.block,
  ICC.block,
  variable.p = FALSE,
  
  # correlation params at cluster level
  cv.block = 0.6, 
  size.impact.cor.block = 0,
  proptx.impact.cor.block = 0,
  size.intercept.cor.block = 0,
  
  # cluster-level params
  n.bar, sd.cluster, ICC.cluster = 0.20,
  
  # correlation params at individual level
  cv.cluster = 0.6,
  size.impact.cor.cluster = 0,
  cor.strength.cluster = 0.75,
  
  # method params
  n.runs = 3,
  just.describe.data = FALSE, 
  just.return.data = FALSE,
  ... ) {
  
  df <- generate_blocks(
    K = K, 
    J.bar = J.bar,
    p.bar = p.bar,
    ATE.block = ATE.block, 
    sd.block = sd.block,
    ICC.block = ICC.block,
    variable.p = variable.p,
    
    # cluster-level params
    n.bar = n.bar, 
    sd.cluster = sd.cluster, 
    ICC.cluster = ICC.cluster,
    
    # size and correlation params
    cv.block = cv.block, 
    cv.cluster = cv.cluster,
    size.impact.cor.block = size.impact.cor.block,
    size.intercept.cor.block = size.intercept.cor.block,
    size.impact.cor.cluster = size.impact.cor.cluster,
    proptx.impact.cor.block = proptx.impact.cor.block)
  
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
  
  # Calculate the estimands
  ATE.student <- with( df, mean( Y1 - Y0 ) )
  # ATE.student = beta_jks %>% 
  #   summarize(n_k = sum(n_jk),
  #             s = sum(n_jk*beta_jk)) %>% 
  #   summarize(ATE.student = sum(s) / sum(n_k))
  ATE.cluster <- mean(beta_jks$beta_jk)
  ATE.dist.student <- beta_jks %>% 
    summarize(n_k = sum(n_jk),
              s = sum(n_jk*beta_jk)) %>% 
    summarize(ATE.dist.student = sum(s/n_k) / n()) %>% 
    pull(ATE.dist.student)
  ATE.dist.cluster <- beta_jks %>% 
    summarize(J_k = n(),
              s = sum(beta_jk)) %>% 
    summarize(ATE.dist.cluster = sum(s/J_k) / n()) %>% 
    pull(ATE.dist.cluster)
  
  # Do finite-sample simulation, rerandomizing n.runs times.  
  ests = map_dfr( 1:n.runs, function( x ) {
    df = rerandomize.data( df )
    get.estimates( df, ... )
  }, .id="subrun" )
  
  # Gather results and ship!
  ests = mutate( ests,
                 ATE.finite.student = ATE.student,
                 ATE.finite.cluster = ATE.cluster,
                 ATE.finite.dist.student = ATE.dist.student,
                 ATE.finite.dist.cluster = ATE.dist.cluster)
  
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










# Test the single.CRT simulation code
if ( FALSE ) {
  
  # NOTE: DB estimators will return infinite SEs if there are
  #  fewer than two tx and co clusters in any block
  res <- single.CRT(
    K=10, J.bar=20, p.bar=0.5,
    variable.p=F,
    ATE.block=0.2, sd.block=0.1, 
    ICC.block=0.2,
    
    size.impact.cor.block = 0,
    proptx.impact.cor.block = -0,
    
    n.bar=50, sd.cluster=0.1, ICC.cluster=0.3,
    include_method_characteristics=F)
  
  head(res)
  table( res$subrun )
}


# Method to run full nested simulation ----

#' Run a simulation with R trials and return all the simulation runs
#' as a large dataframe.
#'
#' @param sd.cluster Treatment variation at the cluster level.
#' 
#' @param n.runs Number of "finite sample" replicates for each
#'   generated dataset.
#' @param R Number of datasets.
run.scenario = function( 
  # block-level params
  K, J.bar, p.bar,
  ATE.block, sd.block,
  ICC.block,
  variable.p,
  
  # correlation params at cluster level
  cv.block = 0.6,
  size.impact.cor.block,
  size.intercept.cor.block = 0,
  proptx.impact.cor.block,
  
  # cluster-level params
  n.bar, sd.cluster, ICC.cluster = 0.2,
  
  # correlation params at individual level
  cv.cluster = 0.6,
  size.impact.cor.cluster,
  cor.strength.cluster = 0.75,
  
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
        block ATE=%.2f\tblock sd=%.2fblock ICC=%.2f\t
        block size correlate=%.2f\tblock proptx correlate=%.2f
        n.bar=%d\tcluster sd=%.2f\tICC=%.2f
        cluster size correlate=%s
        R=%d\n", 
        ID, K, J.bar, p.bar,
        ATE.block, sd.block, 
        ICC.block,
        size.impact.cor.block, proptx.impact.cor.block,
        n.bar, sd.cluster, ICC.cluster,
        size.impact.cor.cluster,
        R)
  
  # Get all the ... arguments in list
  dots = list(...)
  scat( "\tExtra arguments: %s\n", paste(names(dots), collapse=", ") )
  
  rps = map_dfr( 1:R, \(i) single.CRT( 
    K = K, 
    J.bar = J.bar,
    p.bar = p.bar,
    ATE.block = ATE.block, 
    sd.block = sd.block,
    ICC.block = ICC.block,
    variable.p = variable.p,
    
    # cluster-level params
    n.bar = n.bar, 
    sd.cluster = sd.cluster, 
    ICC.cluster = ICC.cluster,
    
    # size and correlation params
    cv.block = cv.block, 
    cv.cluster = cv.cluster,
    size.impact.cor.block = size.impact.cor.block,
    size.intercept.cor.block = size.intercept.cor.block,
    size.impact.cor.cluster = size.impact.cor.cluster,
    proptx.impact.cor.block = proptx.impact.cor.block,
    cor.strength.cluster = cor.strength.cluster,
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
                     ATE.block = 0.2, 
                     sd.block = 0.1,
                     ICC.block = 0.2,
                     variable.p = T,
                     
                     # cluster-level params
                     n.bar = 50, 
                     sd.cluster = 0.1, 
                     ICC.cluster = 0.2,
                     
                     # size and correlation params
                     size.impact.cor.block = 0,
                     size.intercept.cor.block = 0.5,
                     size.impact.cor.cluster = 0,
                     proptx.impact.cor.block = 0,
                     
                     R = 2,
                     .progress="text",
                     include_method_characteristics=F)
  head( rr )  
}



