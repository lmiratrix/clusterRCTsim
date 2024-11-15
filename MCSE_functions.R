
# Utility functions called by "table_Section5.4.R" to calculate monte
# carlo standard errors.

# Written by Xiang Meng


kurtosis <- function(x){
  S_T = sd(x)
  kurt = mean( (x - mean(x))^4 ) / S_T^4
  return(kurt)
}


MCSE_SE <- function(x){
  S_T = sd(x)
  R <- length(x)
  k_T <- kurtosis(x)
  return( sqrt( S_T^2 * sqrt( (k_T-1)/R )) )
}

MCSE_bias <- function(x){
  S_T2 = var(x)
  R <- length(x)
  return( sqrt( S_T2 / R ) )
}
