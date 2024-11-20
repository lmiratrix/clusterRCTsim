
# Utility functions for the simulation code



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
  if ( is.null( df$did ) || length( unique( df$did ) ) == 1 ) {
    ests = safe_comp( Yobs ~ Z | sid, data=df, ... )
  } else {
    ests = safe_comp( Yobs ~ Z | sid | did, data=df, ... )
  }
  if ( is.null( ests$error ) ) {
    ests = ests$result
    ests = rename( ests, 
                   ATE.hat = ATE_hat,
                   SE.hat = SE_hat )
    
  } else {
    cat( "\nCaught error in compare_methods:\n" )
    print( ests$error )
    browser()
    saveRDS( df, file=here::here( "../clusterRCTsim/data_triggers_warning.rds" ) )
    ests = data.frame( ATE.hat = c(NA),
                       SE.hat = c(NA) )    
  }
  ests
}


# Reassign treatment to all clusters within blocks.
rerandomize.data = function( dat ) {
  
  stopifnot( all( c("sid", "Z", "Y1", "Y0" ) %in% colnames(dat) ) )
  
  grid <- attr( dat, "district_key" )
  
  if ( "did" %in% names(dat) ) {
    dat <- dat %>% group_by( did ) %>% 
      mutate( p = mean( Z[ !duplicated( sid ) ] ) ) %>%
      ungroup() %>%
      mutate( Z = randomizr::block_and_cluster_ra( blocks=did,
                                                   clusters=sid, 
                                                   prob_unit=p),
              Yobs = ifelse(Z, Y1, Y0) ) %>% 
      dplyr::select(-p)
  } else {
    # No blocking
    dat <- dat %>% 
      mutate( p = mean( Z[ !duplicated( sid ) ] ),
              Z = randomizr::cluster_ra( clusters=sid, 
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
describe.data = function( df, 
                          # formula = NULL, 
                          print=TRUE ) {
  
  # we want:
  #  - true ATEs
  #  - true correlations from district_key
  #  - etc...
  
  # district params
  d_key = attr(df, "district_key")
  
  #if ( !is.null( formula ) ) {
  #  df = clusterRCT:::make_canonical_data(formula, data=df )
  #}
  
  # Why this needed?
  #df = rerandomize.data( df )
  
  if ( !( "did" %in% colnames(df) ) ) {
    df$did = "single"
  }
  
  if ( !is.null( d_key ) ) {
    sstat_distr <- d_key %>% 
      summarize(cor.J.tau_dist = cor(J, district.tau),
                cor.p.tau_dist = cor(p, district.tau),
                sd.tau_dist = sd(district.tau),
                sd.J_dist = sd(J),
                sd.p_dist = sd(p))
    
    if ( print ) {
      print(sstat_distr)
    }
  } 
  
  
  # compute true ATE values
  ATEs <- df %>% 
    group_by(did, sid) %>% 
    summarize(Njk = n(),
              Bjk = mean(Y1-Y0),
              .groups="drop_last") %>% 
    summarize(Nk = sum(Njk),
              Jk = n(),
              ATE.clust = mean(Bjk),
              ATE.indiv = weighted.mean(Bjk, w=Njk),
              .groups="drop") %>% 
    summarize(ATE.distr = mean(ATE.clust),
              ATE.clust = weighted.mean(ATE.clust, w=Jk),
              ATE.indiv = weighted.mean(ATE.indiv, w=Nk) )
  
  
  # compute design-based point ATE estimates
  #  - district-weighted, cluster-weighted, and individual-weighted
  ATEhats <- df %>% 
    group_by(did, sid) %>% 
    summarize(Njk = n(),
              Z = mean(Z),
              Yhat.clust = mean(Yobs),
              .groups="drop_last")
  
  stopifnot( all( ATEhats$Z %in% c(0,1) ) )
  
  ATEhats <- ATEhats %>%
    summarize( Jk = length(unique(sid)),
               Nk = sum(Njk),
               Bhat.k.clust = sum(Z*Yhat.clust) / sum(Z) -
                 sum((1-Z)*Yhat.clust) / sum(1-Z),
               Bhat.k.indiv = sum(Njk*Z*Yhat.clust) / sum(Z*Njk) -
                 sum(Njk*(1-Z)*Yhat.clust) / sum((1-Z)*Njk)) %>% 
    summarize(ATE.distr = mean(Bhat.k.clust),
              ATE.clust = weighted.mean( Bhat.k.clust, w = Jk ),
              ATE.indiv = weighted.mean( Bhat.k.indiv, w = Nk ) )
  
  res <- bind_rows(ATEs, ATEhats) %>% 
    mutate(ATE = c("True ATEs", "Estimated ATEs"), .before=everything())
  
  if ( print ) {
    print(res)
    invisible( res ) 
  } else {
    res
  }
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



#' Taken from the clusterRCT::describe_clusterRCT() method, this takes
#' a list of parameter values in x and prints them out in a nicely
#' formatted block.
#' 
#' @param x List of parameter values (or estimates)
print_stats <- function( x ) {
  
  if (  is.null( x$K ) || x$K == 1 ) {
    clusterRCT:::scat( "Cluster RCT: %d units in %d clusters\n", x$n, x$J )
  } else {
    clusterRCT:::scat( "Cluster RCT: %d units in %d clusters across %d blocks\n",
                       x$n, x$J, x$K )
  }
  
  
  if ( !is.null( x$K ) && x$K != 1 ) {
    clusterRCT:::scat( "Block Statistics:\n\tAvg units/block: %.2f (coef var %.2f)\n\tAvg clusters/block: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f-%.1f w/ IQR = %.1f\n\tICC: %.2f\n\tAvg tx: %.2f (coef var %.2f)\n\t25-75 Quantile tx: %.2f-%.2f w/ IQR %.2f\n",
                       x$n_block, x$n_block_cv,
                       x$Jbar, x$Jcv,
                       x$J.25, x$J.75, x$J.IQR,
                       x$block_ICC,
                       x$tx.avg, x$tx.cv,
                       x$tx.25, x$tx.75, x$tx.IQR)
    if ( hasName( x, "R2.3" ) ) {
      clusterRCT:::scat( "\tR2.3: %.2f (%d fixed effects)\n", x$R2.3, x$K)
    }
  }
  
  clusterRCT:::scat( "Cluster Statistics:\n\tAvg size: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f-%.1f w/ IQR = %.1f\n\tICC: %.2f\n\tProp treated: %.2f\n",
                     x$nbar, x$ncv,
                     x$n.25, x$n.75, x$n.IQR,
                     x$cluster_ICC,
                     x$p.tx )
  if ( hasName( x, "R2.2" ) ) {
    clusterRCT:::scat( "\tR2.2: %.2f (%d covariates)\n", x$R2.2, x$ncov.2)
  }
  
  
  clusterRCT:::scat( "Unit Statistics:\n\tstddev( Y0 ): %.2f\n",
                     x$sdY0 )
  if ( hasName( x, "R2.1" ) ) {
    clusterRCT:::scat( "\tR2.1: %.2f (%d covariates)\n", x$R2.1, x$ncov.1)
  }
  
  if ( ".missing" %in% names(x) ) {
    clusterRCT:::scat( "missing data counts:\n" )
    print( x$.missing[[1]] )
  }
  
  notes = attr( x, "notes" )
  if ( !is.null(notes) && ( length( notes ) > 0 ) ) {
    cat( "Notes:\n" )
    for ( n in notes ) {
      cat( "\t", n, "\n" )
    }
  }
  
  cat( glue::glue( "{x$num_singletons} singletons, {x$num_doubletons} doubletons." ) )
  cat( "\n" )
  if ( x$matched_pairs ) {
    cat( "Matched pairs experiment\n" )
  }
  
  invisible( x )
  
}


