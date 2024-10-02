
# was getting error running devtools::install_github(), so used this approach instead
# download.file("https://github.com/lmiratrix/clusterRCT/archive/refs/heads/main.zip", destfile = "clusterRCT.zip")
# unzip("clusterRCT.zip")
# install.packages("clusterRCT-main", repos = NULL, type = "source")

# load libraries ----
# devtools::install_github("lmiratrix/clusterRCT")
library(tidyverse)
library(haven)
library(clusterRCT)
library(tools)
library(writexl)



# analyze cleaned real data ----
if ( FALSE ) {
  f <- "/data/share/pse/IES CRT/sd7/clean study datasets/beaconmentor.RDS"
  assign(file_path_sans_ext(basename(f)),
         readRDS(f))
  
  
  # Run compare_methods() ----
  
  ## Beacon Mentor ----
  
  # outcomes_beaconmentor <- c('S8TGPA', 'S8TOTCR', 'SRPRST2')
  # covars_beaconmentor <- c('RADEV', 'S8TEACHID', 'S8FULTMEC', 'ZFemale', 'RATIME', 'ZAGELT24', 'PTTOP75C', 'PTBTM25C')
  
  
  # example with no covariates
  cm_bm_nocovar <- compare_methods(S8TOTCR ~ racode_b | S8CLUSTE | RABLOCK,
                                   data = beaconmentor,
                                   include_method_characteristics = F)
  
  # example with one covariate
  cm_bm_onecovar <- compare_methods(S8TOTCR ~ racode_b | S8CLUSTE | RABLOCK,
                                    data = beaconmentor,
                                    control_formula = ~ ZFemale,
                                    include_method_characteristics = F)
  
  
  cc = bind_rows( best = cm_bm_onecovar, none = cm_bm_nocovar, .id = "Covariates" )
  
  cc %>%
    filter( method %in% c( "Agg_FI_Cluster_Cluster", "DB_FI_Cluster_Cluster" ) ) %>%
    dplyr::select( Covariates, method, df, ATE_hat, SE_hat ) %>%
    pivot_wider( names_from = method, values_from = c( ATE_hat, SE_hat ) ) %>%
    rename_with( ~ str_replace(., "_Cluster_Cluster", "" ) ) %>%
    mutate( discATE = is.na( ATE_hat_Agg_FI ) != is.na( ATE_hat_DB_FI ) ) %>%
    mutate( discSE = is.na( SE_hat_Agg_FI ) != is.na( SE_hat_DB_FI ) )
}



# Make fake data with beaconmentor structure ----

structure <- read_rds( here::here( "empirical/results/beaconmentor_J-K-T.RDS" ) )
structure
names(structure) <- c( "clusterID", "blockID", "Z" )
table( structure$blockID, structure$Z )
structure$N = 20
length( unique( structure$blockID ) )
length( unique( structure$clusterID ) )

# Get description of data
mstr <- read_rds(here::here( "empirical/results/master240528_combined.RDS" ) )

parameters <- filter( mstr, Outcome == "S8TGPA", Covariates=="best" ) %>%
  slice_head(n=1) %>%
  as.list()



# Some helper methods ----

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



generate_data_from_fingerprint <- function( structure, parameters, 
                                            tau = 0.20 ) {
  
  K = length( unique( structure$blockID ) )
  J = length( unique( structure$clusterID ) )  
  n = nrow( structure )
  structure <- mutate( structure,
                       blockID = as.character(blockID),
                       clusterID = as.character(clusterID ) )
  
  block_int = rnorm( K, mean = 0, sd = sqrt( parameters$block_ICC ) )
  names( block_int ) = unique( structure$blockID )
  
  cluster_int = rnorm( J, mean = 0, sd = sqrt( parameters$cluster_ICC ) )
  names( cluster_int ) = unique( structure$clusterID )
  
  W = rnorm( J )
  names(W) = names( cluster_int )
  X = rnorm( n )
  
  beta_W = sqrt( parameters$R2.2 * parameters$cluster_ICC )
  beta_X = sqrt( (1-parameters$block_ICC-parameters$cluster_ICC) * parameters$R2.1 )
  
  sigma_e = sqrt( ( 1 - parameters$block_ICC - parameters$cluster_ICC ) * (1-parameters$R2.1 ) )
  
  
  # TODO: Possibly group mean center individual level covariate.  Good
  # idea?
  dat <- structure %>% 
    mutate( e = rnorm( n(), mean = 0, sd = sigma_e ),
            X = X ) %>%
    group_by( clusterID ) %>%
    mutate( X = X - mean(X) ) %>%
    ungroup() %>%
    mutate(
      Y = block_int[ blockID ] + cluster_int[ clusterID ] + Z * tau + e + W[clusterID] * beta_W + X * beta_X )
  
  # Make it so we end up with right number of covariates
  if ( parameters$ncov.2 > 0 ) {
    dat$X = dat$X + W[ dat$clusterID ]
  }
  
  X = blkvar:::make_covariate_set( dat$X, k = parameters$ncov.1, 
                                   var_prefix = "X" )
  dat <- bind_cols( dat, X )
  
  if ( parameters$ncov.2 > parameters$ncov.1 ) {
    W = blkvar:::make_covariate_set( rnorm(J), 
                                     k = parameters$ncov.2 - parameters$ncov.1, 
                                     var_prefix = "W" )
    W$clusterID = names( cluster_int )
    dat <- left_join( dat, W, by="clusterID" )
  }
  
  dat
  
}


if ( FALSE ) {
  
  structure = tibble( clusterID = 1:8,
                      Z = rep( c(0,1), 4 ) )
  structure = bind_rows( structure, structure, structure )
  s2 = structure
  s2$clusterID = rep( 1:4, nrow(s2)/4 )
  structure = bind_rows( structure, s2 )
  nrow( structure )
  structure$blockID <- round( 1/3 + structure$clusterID / 3 )
  
  table( structure$blockID, structure$clusterID )
  undebug( generate_data_from_fingerprint )
  dat <- generate_data_from_fingerprint( structure,
                                         parameters = list( ncov.1 = 2,
                                                            ncov.2 = 1,
                                                            block_ICC = 0.8,
                                                            cluster_ICC = 0.0,
                                                            R2.2 = 0.0,
                                                            R2.1 = 0.9 ) )
  
  ggplot( dat, aes( X1, Y, col=clusterID ) ) +
    facet_wrap( ~blockID ) +
    geom_point()
}


# Check methods on fake beacon dataset ----



beac = generate_data_from_fingerprint( structure,
                                       parameters )
beac

# Note--the individual covariates give us more cluster covariates.
# Adjust DGP?
dd <- describe_clusterRCT( Y ~ Z | clusterID | blockID, data = beac,
                     control_formula = ~ X1 + X2 )
dd

# Compare to the listed parameters
print_stats( parameters )


cc <- compare_methods( Y ~ Z | clusterID | blockID, data = beac,
                           control_formula = ~ X1 + X2,
                       include_method_characteristics = F)
cc

cc2 <- compare_methods( Y ~ Z | clusterID | blockID, data = beac,
                       include_method_characteristics = F)

cc = bind_rows( best = cc, none = cc2, .id = "Covariates" )

cc %>%
  filter( method %in% c( "Agg_FI_Cluster_Cluster", "DB_FI_Cluster_Cluster" ) ) %>%
  dplyr::select( Covariates, method, df, ATE_hat, SE_hat ) %>%
  pivot_wider( names_from = method, values_from = c( ATE_hat, SE_hat ) ) %>%
  rename_with( ~ str_replace(., "_Cluster_Cluster", "" ) ) %>%
  mutate( discATE = is.na( ATE_hat_Agg_FI ) != is.na( ATE_hat_DB_FI ) ) %>%
  mutate( discSE = is.na( SE_hat_Agg_FI ) != is.na( SE_hat_DB_FI ) )

