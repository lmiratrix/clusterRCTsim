
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

# Secure server data root; override via env var for a different mount point.
IES_CRT_ROOT <- Sys.getenv( "IES_CRT_DATA_ROOT", "/data/share/pse/IES CRT" )

# analyze cleaned real data ----
if ( FALSE ) {
  f <- file.path( IES_CRT_ROOT, "sd7/clean study datasets/beaconmentor.RDS" )
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
  if ( FALSE ) {
    # For looking at DB issues
    debug( clusterRCT:::schochet_variance_formula_block )
  }
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

