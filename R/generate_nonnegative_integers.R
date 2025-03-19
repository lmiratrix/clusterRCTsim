

# Generate cluster sizes or other sets of nonnegative integers.
#
# This script compares the coefficient of variation to the mean, and
# then generates from one of three distributions:
#
# A) sigma2 < mean: binomial, selecting n and p to make the mean and
#                   sigma2 the desired targets.
# B) sigma2 = mean: poisson 
# C) sigma2 > mean: negative binomial
#
# For (A) the variance is initially approximate.  In order to get the
# exact variance, we need to spike in a two-size distribution
# surrounding the target mean.  This spike reduces the variance to the
# desired level.  (It appears as if the approximation from the
# binomial always gives a variance a bit too high, making this work.)
#
# There are some combinations of very low variance with non-integer
# means that appear to be impossible. E.g., mean = 0.5 and variance of
# 0.01, since even a coin flip on 0 or 1 has a variance of at least
# 0.5^2 = 0.25.  In these cases, the function will give the spike,
# only, and also produce a warning.





#' Generate n nonnegative integers from a distribution with given mean
#' and standard deviation.
#' 
#' @return List of integers 
generate_nonnegative_integers <- function(n, mean, sd, verbose = FALSE) {
  stopifnot( mean >= 0 )
  stopifnot( sd >= 0 )
  
  if ( sd == 0 ) {
    if ( round(mean) != mean ) {
      warning( glue::glue( "Cannot have 0 variation sizes that average out to {round(mean,digits=3)}." ) )
    }
    return( rep( round(mean), n ) )
  }
  
  var <- sd^2
  
  # Calculate the variance of a two-size spike
  mu_f = floor( mean )
  mu_c = ceiling( mean )
  p_spike = mean - mu_f
  var_spike = p_spike * (1 - p_spike)
  if ( var_spike > var ) {
    if ( mu_f != mean ) {
      warning( glue::glue( "Cannot generate a distribution with this mean ({round(mean,digits=3)}) and sd ({round(sqrt(var),digit=3)}).  Achieved sd={round(sqrt(var_spike),digits=3)})." ) )
      
      return( mu_f + rbinom( n, 1, p_spike ) )
    }
  }
  
  
  if ( var < mean*0.9999 ) {
    # Use binomial with some n
    p <- (mean - var)/mean
    bin_n <- mean / p
    
    # Round n to an integer (ideally quite large)
    bin_n <- ceiling(bin_n)
    
    # Recalculate p to maintain the exact mean
    p <- mean/bin_n
    
    pvar = bin_n * p * (1-p)
    if ( pvar < var ) {
      stop( "Error in binomial calculation" )
    }
    if ( pvar > var * 1.0001 ) {
      # We need to reduce our variance by spiking in our spike
      
      p_use_spike = (var - pvar) / (var_spike - pvar)
      if ( p_use_spike < 0 ) {
        stop( "Error in spike calculation" )
      }
      
      draws = rbinom( n, bin_n, p )
      
      Z = rbinom( n, 1, p_use_spike )
      n_spk = sum(Z)
      draws[Z == 1] = mu_f + rbinom( n_spk, 1, p_spike )
      
      if ( verbose ) {
        attr( draws, "spike" ) = p_use_spike
      }
      return( draws )
      
    } else {
      return( rbinom( n, bin_n, p ) )
    }    
  } else if ( var < mean*1.0001 )  {
    return( rpois( n, lambda = mean ) )
  } else {
    # negative binomial when var > mean; this one is exact
    
    size <- mean^2 / (var - mean)
    prob <- mean / var
    rr <- rnbinom(n, size = size, prob = prob)
    
    return( rr )
  }
}

# SOME NOTES:

# For poisson the CV is the following
# sqrt( lambda ) / lambda = 1 / sqrt( lambda )


# Negative binomial N( r, p )
# mean = r( 1 - p ) / p
# variance = r( 1 - p ) / p^2
# CV = r(1-p)/p / sqrt( r(1-p) / p^2 ) 
#    = sqrt(r(1-p))
#
# mean = var / p
# p = var / mean
# and r = (1 - var/mean) / (var/mean) = (mean - var) / var 
#.    = mean / var - 1
# So r = p * mean / ( 1 - p )
#    r = p^2 / (1-p)
# so p*mean/(1-p) = p^2 / (1-p)
# 



if ( FALSE ) {
  
  # Test the function
  
  check_function <- function( mean, cv, R = 10000 ) {
    
    rps = generate_cluster_sizes( R, mean, cv*mean )
    
    tibble( EY = mean(rps),
            sdY = sd(rps),
            cvY = sdY / EY,
            min = min(rps),
            max = max(rps) )
  }
  
  check_function( 0.5, 0, 10000 )
  check_function( 10, 0.5, 10000 )
  check_function( 100, 1, 10000 )
  check_function( 100, 0.1, 10000 )
  check_function( 100, 0.02, 10000 )
  
  gg <- generate_cluster_sizes( 10000, 100, 0.02, verbose=TRUE )
  attributes( gg )
  qplot( gg )
  
  library( tidyverse )
  grid = expand_grid( 
    mean = c( 0.5, 1, 1.5, 3, 5, 10, 30, 50, 100, 200 ),
    cv = seq( 0, 1.5, by = 0.1 ) 
  ) %>% 
    mutate( sigma2 = (mean * cv)^2,
            small_sigma = sigma2 < mean )
  grid
  
  grid <- grid %>%
    mutate( res = map2( mean, cv, check_function, R=10000 ) ) %>%
    unnest( cols=res )
  
  
  grid
  summary( grid$EY - grid$mean )
  qplot( grid$EY - grid$mean )
  filter( grid, abs( EY - mean ) / mean > 0.02 ) %>%
    print( n = 100 )
  
  summary( grid$cv - grid$cvY )
  qplot( grid$cv - grid$cvY )
  filter( grid, abs( cv - cvY ) > 0.05 ) %>%
    print( n = 100 )
  
}





# Generate random probabilities ----



generate_treatment_probabilities <- function(n, mean, sd) {
  var = sd^2
  
  alpha <- mean * ((mean * (1 - mean) / var) - 1)
  beta <- (1 - mean) * ((mean * (1 - mean) / var) - 1)
  
  if (alpha <= 0 || beta <= 0) {
    stop( glue::glue( "Invalid parameters of {mean}, {var}: Check that mean and variance are feasible within [0,1]") )
  }
  
  rbeta(n, alpha, beta)
}


if ( FALSE ) {
  # Example: Generate 1000 numbers with mean = 0.5 and variance = 0.02
  set.seed(123)
  random_numbers <- generate_treatment_probabilities(100000, mean = 0.5, var = 0.2)
  
  # Check mean and variance
  mean(random_numbers)
  var(random_numbers)
  
}

