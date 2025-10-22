

# This file contains functions to generate data for simulations.

# In particular the code has a bunch of simulation utilities for
# making data with specific block structures

# It also has code for running canonical simulations of various types.




library( tidyverse )
library( furrr )


# PUMP based DGP code ----

# Default set of parameters for the DGP
model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 10                          # number of schools
    , K = 4                          # number of districts
    , nbar = 10                       # number of individuals per school
    , S.id = NULL                     # N-length vector of school assignments
    , D.id = NULL                     # N-length vector of district assignments
    , Xi0 = 0                         # scalar grand mean outcome under no treatment
    , MDES = 0.2            # minimum detectable effect size
    , R2.3 = 0.1              # percent of district variation
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.2           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.2          # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)


# Update parameters for specific treatment effect and desired
# heterogeneity
update_model_params <- function( ATE = 0.2, tx_het = FALSE ) {
    params = model.params.list
    params$MDES = ATE
    if ( !tx_het ) {
        params$omega.2 = 0
        params$omega.3 = 0
    }
    params
}

print_list <- function( lst ) {
    lst <- lst %>%
        map( ~ ifelse( is.null(.x), "NULL", .x ) )
    as_tibble( lst ) 
}

print_list(model.params.list)


#### Function to make data via PUMP methods ####

#'
#' Make blocked cluster RCT where each block has two sizes of cluster
#' defined by n1 and n2, with different treatment effects.
#'
#' @param n1 Size of small clusters.
#' @param n2 Size of large clusters
#' @param ATE1 Average treatment effect for small clusters.
#' @param ATE2 Average treatment effect for large clusters.
#' @param blocks TRUE means multiple districts (blocks).  FALSE means
#'   single block.
#' @param het_block TRUE means blocks are different shapes and styles,
#'   according to a specific and brittle recipe.  WARNING: this is
#'   brittle and may not work for all values of params passsed.
#'
#'   NOTE: THIS METHOD IS CURRENTLY LESS PREFERRED
make_two_tier_PUMP_data <- function( n1, ATE1, n2 = 5*n1, ATE2 = 3*ATE1,
                                     blocks = TRUE, params = model.params.list,
                                     het_block = FALSE ) {
    
    if ( het_block && !blocks ) {
        stop( "Invalid--can't have het blocks with no blocks" )
    }
    
    if ( !blocks ) {
        # Force the three level stats out
        params$ICC.3 = NULL
        params$R2.3 = NULL
        params$omega.3 = NULL
        design = "d2.2"
        params$J = params$J * params$K
        params$K = NULL
    } else {
        # Blocked, cluster-randomized design
        design = "d3.2_m3ff2rc"
    }
    
    # Make small clusters across the districts
    N = max(n1, n2)
    if ( n2 < n1 ) {
        stop( "n2 must be larger than n1" )
    }
    params$MDES = 0
    params$nbar = N
    
    sim.data <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.30,
                                    include_POs = TRUE )
    sim.data.sm = filter( sim.data, as.numeric( S.id ) %% 2 == 0 )
    sim.data.lg = filter( sim.data, as.numeric( S.id ) %% 2 == 1 )
    
    if ( !blocks ) {
        sim.data.sm$D.id = "sing"
        sim.data.lg$D.id = "sing"
    }    
    
    sim.data.sm <- sim.data.sm %>%
        group_by( D.id, S.id ) %>%
        slice_head( n = n1 ) %>%
        mutate( Y1 = Y0 + ATE1,
                Yobs = ifelse( T.x == 1, Y1, Y0 ) ) %>%
        ungroup()
    
    sim.data.lg <- sim.data.lg %>%
        mutate( Y1 = Y0 + ATE2,
                Yobs = ifelse( T.x == 1, Y1, Y0 ),
                S.id = paste0( "b-", S.id ) )
    
    
    # Bind so each block has small and large clusters
    sim.data = bind_rows( sim.data.sm, sim.data.lg )
    
    sim.data$Tx = sim.data$T.x
    sim.data$T.x = NULL
    
    if ( het_block ) {
        stop( "This code broken" )
        # Add a third size of blocks in District 1, and also change
        # the treatment proportion
        params$K = 1
        params$J = round( 1.5*params$J )
        params$MDES = ATE1 + ATE2
        params$nbar = n1 + 2
        sim.data3 <- PUMP::gen_sim_data( d_m = design, params,
                                         Tbar = 0.70, include_POs = TRUE )
        
        sim.data3$D.id = sim.data.small$D.id[[1]]
        sim.data3$S.id = paste0( "h-", sim.data3$S.id )
        sim.data = bind_rows( sim.data, sim.data3 )
        
        
        sim.data <- filter( sim.data, 
                            !( S.id %in% c( paste0( "b-", 11:16 ) ) ) )
        
        # browser()
        if ( FALSE ) {
            ss <- sim.data %>% group_by( D.id, S.id ) %>%
                summarise( n = n(), .groups = "drop" )
            table( ss$n, ss$D.id )
            sim.data %>% group_by( S.id, D.id, Tx ) %>%
                summarise( n = n(), .groups = "drop",
                           Tx = mean( Tx ) ) %>%
                group_by( D.id ) %>%
                summarise( p_tx = mean( Tx ),
                           sizes = paste0( unique( n ), collapse=", " ),
                           nT = sum( n ),
                           J = n(), .groups = "drop" )
        }
        
        # Rerandomize units within the blocks with different numbers randomized.
        sim.data$Tx = randomizr::block_and_cluster_ra( blocks = sim.data$D.id,
                                                       clusters = sim.data$S.id,
                                                       block_m = c( 10, 6, 10, 5 ) )
        sim.data <- mutate( sim.data,
                            Yobs = ifelse( Tx == 1, Y1, Y0 ) )
    } 
    
    # Standardize outcome to effect size units
    sdY0 = sd( sim.data$Y0 )
    sim.data <- mutate( sim.data,
                        Yobs = Yobs / sdY0,
                        Y0 = Y0 / sdY0,
                        Y1 = Y1 / sdY0 )
    
    
    as_tibble( sim.data )
}




##### Testing make_two_tier_PUMP_data #####
if ( FALSE ) {
    dat = make_two_tier_PUMP_data( 10, 0.2, het_block = FALSE)
    model.params.list
    head( dat )
    slice_sample( dat, n=10 )
    table( dat$D.id )
    table( table( dat$S.id ) )
    sd( dat$Y0 )
    
    describe_clusterRCT( Yobs ~ Tx | S.id | D.id, data=dat )
    compare_methods(  Yobs ~ Tx | S.id | D.id, data=dat ) %>%
        arrange( ATE_hat )
    
    # Person average ATE
    mean( dat$Y1 - dat$Y0 )
    table( round( dat$Y1 - dat$Y0, digits=5 ) )
    dat %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ) ) %>%
        pull( tauk ) %>% mean()
    
    
    
    # More varied data to tease apart methods
    dat = make_two_tier_PUMP_data( 10, 0.2, het_block = TRUE)
    dat
    table( dat$D.id )
    table( table( dat$S.id ) )
    
    ad <- dat %>% group_by( D.id, S.id ) %>%
        summarise( n = n(),
                   Z = mean( Tx ) )
    table( ad$D.id, ad$n )
    sd( dat$Y0 )
    
    describe_clusterRCT( Yobs ~ Tx | S.id | D.id, data=dat )
    compare_methods(  Yobs ~ Tx | S.id | D.id, data=dat ) %>%
        arrange( ATE_hat )
    
    
    
    
    # Person average ATE
    mean( dat$Y1 - dat$Y0 )
    table( round( dat$Y1 - dat$Y0, digits=5 ) )
    dat %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ) ) %>%
        pull( tauk ) %>% mean()
    
    
    
    
    dat = make_two_tier_PUMP_data( 10, 0.2, blocks=FALSE )
    dat
    dat$D.id
    table( table( dat$S.id ) )
    sd( dat$Y0 )
}


# Block structure based DGP code ----


#### Function to make data by specifying mixture of blocks and cluster sizes ####

#' Make data with very different blocks that have cluster mixtures
#' that correlate with treatment effects.
#'
#' Each block will be a mix of clusters of specified sizes, with
#' specified number of clusters of each size.  The clusters will each
#' have a mean outcome for that cluster, and a cluster-specific
#' treatment effect.  I.e., blocks are a mixture of a discrete number
#' of cluster types.
#'
#' The initial "Randomizations" are perfectly balanced, so estimates
#' will focus on bias, not estimation error.
#'
#' @param n List of lists of types of cluster size within each block,
#'   or a pre-computed dataframe of the block and cluster statistics.
#' @param J List of lists of number of clusters of each type.
#' @param Y0 List of lists of mean control-side outcome of each type
#'   of cluster.
#' @param ATE List of lists of average treatment effects for each
#'   cluster type.
#' @param ptx Proportion treated.  Can be a list of what proportion to
#'   treat within each block.
#'
#' 
make_block_data <- function( n, J, Y0 = 0, ATE, ptx = 0.5 ) {
    
    statblock = NA
    if ( is.data.frame( n ) ) {
        statblock = n
        
        if ( !( "ptx" %in% colnames( statblock ) ) ) {
            statblock$ptx = 0.5
        }
    } else {
        
        if ( length(ATE) == 1 ) {
            ATE = rep( ATE, length(J) )
        }
        if ( length(ptx) == 1 ) {
            ptx = rep( ptx, length(J) )
        }
        if ( length(Y0) == 1 ) {
            Y0 = rep( Y0, length(J) )
        }
        
        stopifnot( length(n) == length(J) )
        stopifnot( length(n) == length(Y0) )
        stopifnot( length(n) == length(ATE) )
        stopifnot( length(n) == length(ptx) )
        
        statblock = list( n=n, J=J, Y0=Y0, ATE=ATE, ptx=ptx )
        
    }
    
    # Make J clusters of size n with an ATE of given ATE.
    make_blocks <- function( n, J, Y0, ATE, ptx = 0.5, prefix ) {
        dat = tibble( S.id = rep( 1:J, each=n ),
                      Y0 = Y0 + rnorm( n*J, sd=1 ) )
        dat$Y1 = dat$Y0 + ATE
        dat$Tx = 0 + (dat$S.id %% round(1/ptx) == 0)
        dat$Yobs = ifelse( dat$Tx, dat$Y1, dat$Y0 )
        dat$S.id = paste0( prefix, "-", dat$S.id )
        dat
    }
    
    # Given vectors of n, J, Y0, and ATE, make batches of clusters of
    # given sizes via pmap
    make_blocks_v <- function( n, J, Y0, ATE, ptx=0.5 ) {
        stopifnot( length(n) == length(J) )
        if ( length(ATE) == 1 ) {
            ATE = rep( ATE, length(J) )
        }
        stopifnot( length(n) == length(ATE) )
        if ( length(Y0) == 1 ) {
            Y0 = rep( Y0, length(J) )
        }
        stopifnot( length(n) == length(Y0) )
        if ( length(ptx) == 1 ) {
            ptx = rep( ptx, length(J) )
        }
        stopifnot( length(n) == length(ptx) )
        if ( length( Y0 ) == 1 ) {
            Y0 = rep( Y0, length(n) )
        }
        
        blks = tibble( n=n, J=J, Y0=Y0, ATE=ATE, ptx=ptx,
                       prefix=letters[ 1:length(n) ] )
        res <- pmap_df( blks,
                        make_blocks )
        res
    }
    
    dats = pmap( statblock, make_blocks_v )
    names(dats) = LETTERS[ 1:length(dats) ]
    dat = bind_rows( dats, .id="D.id" )
    
    # Standardize
    sdY0 = sd( dat$Y0 )
    sim.data <- mutate( dat,
                        Yobs = Yobs / sdY0,
                        Y0 = Y0 / sdY0,
                        Y1 = Y1 / sdY0 )
    
    # Add some covariates
    sim.data$X.jk = sim.data$Y0 + rnorm( nrow(sim.data), sd=1)
    sim.data <- sim.data %>%
        group_by( S.id ) %>%
        mutate( X.jk = mean( X.jk ) ) %>%
        ungroup()
    sim.data$C.ijk = sim.data$Y0 + rnorm( nrow(sim.data), sd=1)
    
    as_tibble( sim.data )
}


make_block_data_homog <- function( n, J, K, Y0 = 0, ATE = 0, ptx = 0.5 ) {
    
    make_block_data( n = rep( list(n), K ),
                     J = rep( list(J), K ),
                     Y0 = Y0,
                     ATE = ATE,
                     ptx = ptx )
    
}


# Test make block data methods ----
if ( FALSE ) {
    
    tst <- make_block_data( c( 2, 3 ), J = c( 4, 2 ),
                            Y0 = c( 0, 100 ),
                            ATE = c( 50, 20 ),
                            ptx = c(0.25, 0.5 ) )
    knitr::kable( tst, digits = 1 )
    
    
    tst <- make_block_data( n = list( c( 4, 2 ), 3 ), 
                            J = list( c( 2, 2 ), 4 ),
                            Y0 = c( 0, 100 ),
                            ATE = c( 50, 20 ),
                            ptx = c(0.25, 0.5 ) )
    
    tst %>%
        group_by( D.id, S.id ) %>%
        summarise( Tx = mean( Tx ),
                   Yobs = mean( Yobs ),
                   n = n(), .groups = "drop" ) %>%
        knitr::kable( digits = 1 )
    
    
    
    dd <- tibble(
        n = list( c( 200, 400 ), 300 ),
        J = list( c( 4, 2 ), 2 ),
        Y0 = list( c( 0, 0 ), 0 ),
        ATE = list( c( 1, 2 ), 3 ),
        ptx = c( 0.5, 0.5 ) )
    dd
    
    dd %>%
        make_block_data( ) %>%
        group_by( D.id, S.id ) %>%
        summarise( Tx = mean( Tx ),
                   Yobs = mean( Yobs ),
                   n = n(), .groups = "drop" ) %>%
        print( n = 100 )
    
    
    make_block_data_homog( n=c(2,4), J = c(2,1), K = 2 )
}







#' Make data with very different block structure and cluster mixtures
#' that correlate with treatment effects.
#'
#' "Randomizations" will be perfectly balanced (assuming the
#' proportion of treatment evenly divides each block group size) to
#' focus on bias, not estimation error.
#'
#' @param var_tx Vary proportion units treated.
#' 
make_simple_block_example <- function( version = "A",
                                       tx_scale = 1, var_tx = FALSE ) {
    
    ptx = rep( 0.5, 3 )
    if ( var_tx ) {
        ptx = c( 1/2, 1/6, 1/2 )
    }
    
    if ( version == "A" ) {
        blks <- tribble( ~n, ~ J, ~ Y0, ~ ATE,
                         c( 20, 80 ), c( 2, 4 ), c( 1, 0.5 ), c( 1.4, 0.2 ),
                         c( 20, 80 ), c( 6, 12 ), c( 0, 1.5 ), c( 1.2, 2.4 ),
                         c( 5, 100 ), c( 20, 4 ), c( 4, 0 ), c(  4, -4 ) )
    } else {
        blks <- tribble( ~n, ~ J, ~ Y0, ~ ATE,
                         c( 20, 80 ), c( 2, 4 ), c( 1, 0 ),  c( 0, -2 ),
                         c( 10, 30 ), c( 6, 6 ), c( 0, 2 ), c( 3, 1 ),
                         c( 5, 10 ), c( 10, 10 ), c( 2, 0 ), c(  4, 2 ) )
    }
    blks$ptx = ptx
    blks$ATE = map( blks$ATE, ~ .x * tx_scale )
    
    make_block_data( blks )
}








##### Test make_simple_block_example() #####

if ( FALSE ) {
    dat = make_simple_block_example( version = "B", var_tx = TRUE )
    dat
    
    ds = dat %>%
        group_by( D.id, S.id, Tx ) %>%
        summarise( tauk = mean( Y1 - Y0 ),
                   n = n(), .groups = "drop" )
    ds
    ff <- filter( ds, D.id == "A" )
    ff
    mean( ff$tauk )
    weighted.mean( ff$tauk , ff$n )
    
    # Make table of block-average and person-average impacts by block
    ds <- ds %>%
        group_by( D.id ) %>%
        summarise( tau_c = mean( tauk ),
                   tau_p = weighted.mean( x=tauk, w=n ),
                   #       tau_pp = mean( tauk * n ) / mean(n),
                   p_tx = mean( Tx ),
                   J = n(),
                   Jwt = J * p_tx * (1-p_tx),
                   n = sum( n ),
                   nwt = n * p_tx * (1-p_tx ) )
    ds %>% knitr::kable( digits = 3 )
    
    ds %>%
        summarise( tau_cb = mean( tau_c ),
                   tau_pb = mean( tau_p ),
                   tau_p_FE = weighted.mean( tau_p, w = n*p_tx*(1-p_tx) ),
                   tau_c_FE = weighted.mean( tau_c, w = n*p_tx*(1-p_tx) ),
                   tau_c = weighted.mean( tau_c, w=J ),
                   tau_p = weighted.mean( tau_p, n ),
                   
                   J = sum( J ),
                   n = sum( n ) )
}


# Code for running simulations using the above DGPs ----


#'
#' @param simple If FALSE call PUMP simulation for two-tiered data.
#'   Otherwise is version, e.g., "A", or "B", and will call the
#'   make_simple_block_example() method for DGP.
#' @param blocks Have multiple blocks (only for dgp=PUMP).
#' @param het_block Have heterogeneous blocks (different sizes).  Only
#'   applicable to two-tiered data simulation (dgp=PUMP).
#' @param ICC.2 Intraclass correlation for the second level. Only
#'   applicable to two-tiered data simulation (dgp=PUMP).
#' @param var_tx Vary the proportion of units treated across blocks
#'   (for dgp=A or B only)
#' @param include_covariates Include covariates in the estimating
#'   models.
one_run <- function( dgp = "PUMP", 
                     blocks = TRUE,  
                     het_block = FALSE,
                     ICC.2 = NULL,
                     var_tx = FALSE,
                     include_covariates = FALSE,
                     model.params = model.params.list,
                     add_canonical_weights = TRUE,
                     return_data = FALSE ) {
    
    
    
    # Make the data
    if ( is.function(dgp) ) {
        stopifnot( is.null( ICC.2 ) )
        sim.data = dgp()
    } else if ( dgp != "PUMP" ) {
        stopifnot( is.null( ICC.2 ) )
        sim.data = make_simple_block_example( version = dgp, var_tx = var_tx )
    } else {
        # Make some data
        p = model.params
        if ( !is.null( ICC.2 ) ) {
            p$ICC.2 = ICC.2
        }
        
        sim.data = make_two_tier_PUMP_data( n1 = 10, ATE1 = 0,
                                            n2 = 40, ATE2 = 0.4,
                                            params = p,
                                            blocks = blocks, 
                                            het_block = het_block )
    }
    
    if ( return_data ) {
        return( sim.data )
    }
    
    # Analyze the data
    form = Yobs ~ Tx | S.id | D.id
    if ( !blocks ) {
        form = Yobs ~ Tx | S.id
    }
    control_form = NULL
    if ( include_covariates ) {
        if ( blocks ) {
            control_form = ~ X.jk + C.ijk
        } else {
            control_form = ~ C.ijk
        }
    }
    
    c1 <- clusterRCT::compare_methods( form, data=sim.data,
                                       control_formula = control_form,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = TRUE,
                                       include_dumb = FALSE,
                                       include_method_characteristics = FALSE )
    
    
    # Add the "canonical weight" estimators as a point of comparison
    if ( add_canonical_weights ) {
        c2 <- clusterRCT::canonical_weight_estimators( form, sim.data )
        
        c1 = bind_rows( c1, c2 )
    }
    
    if ( !blocks ) {
        sim.data$D.id = "sing"
    }
    
    
    # Add in the true ATEs
    ds = sim.data %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ),
                   n = n(), .groups = "drop" ) %>%
        group_by( D.id ) %>%
        summarise( tau_c = mean( tauk ),
                   tau_p = weighted.mean( tauk, w=n ),
                   J = n(),
                   n = sum( n ) ) %>%
        summarise( tau_cb = mean( tau_c ),
                   tau_c = weighted.mean( tau_c, w=J ),
                   tau_pb = mean( tau_p ),
                   tau_p = weighted.mean( tau_p, w=n ),
                   J = sum( J ),
                   n = sum( n ) )
    
    c1$tau_c = ds$tau_c
    c1$tau_p = ds$tau_p
    #c1$tau_cb = ds$tau_cb
    #c1$tau_pb = ds$tau_pb
    c1$J = ds$J
    c1$n = ds$n
    c1$sdY0 = sd( sim.data$Y0 )
    
    # Return the results
    c1 %>%
        dplyr::select( -any_of( c( "weight", "df" ) ) )
}



if ( FALSE ) {
    # testing
    one_run( dgp="A" )
    
    one_run( blocks = TRUE, het_block = TRUE )
    
    
    sim.data = make_two_tier_PUMP_data( n1 = 10, ATE1 = 0,
                                        n2 = 40, ATE2 = 0.4,
                                        params = model.params.list,
                                        blocks = TRUE, het_block = TRUE )
    
    c1 <- clusterRCT::compare_methods( form = Yobs ~ Tx | S.id | D.id,
                                       data=sim.data,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = TRUE,
                                       include_method_characteristics = FALSE,
                                       include_dumb = TRUE )
    
    c1
    
    
    # demo of custom dgp function
    my_dgp <- function( ) {
        make_block_data_homog( n=c(2,4), J = c(2,1), K = 2 )
    }
    my_dgp()
    
    one_run( dgp=my_dgp )
    
}


# Code to run full simulation ----

setup_parallel_if_needed <- function() {
    if ( !requireNamespace( "future", quietly = TRUE ) ) {
        stop( "Please install the future package to run simulations in parallel." )
    }
    
    if ( !requireNamespace( "furrr", quietly = TRUE ) ) {
        stop( "Please install the furrr package to run simulations in parallel." )
    }
    
    # Launch system if needed.
    library( furrr )
    if ( !inherits(future::plan(), "multisession") ) {
        ncores = parallel::detectCores() - 2
        cli::cli_alert("Starting {ncores} cores for parallel")
        plan(multisession, workers = ncores )
    } else {
        cli::cli_alert("Already running in parallel, not changing plan")
        
    }
}


calc_performance <- function( rps ) {
    res <- rps %>%
        group_by( method ) %>%
        summarise( EATE = mean( ATE_hat ),
                   tau_p = mean( tau_p ),
                   tau_c = mean( tau_c ),
                   #tau_cb = mean( tau_cb ),
                   #tau_pb = mean( tau_pb ),
                   SE = sd( ATE_hat ),
                   SE_E = sd( ATE_hat ) / sqrt(n()),
                   ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
        mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )
    res
}


     
#' Run a simple simulation with the above DGPs.
#' 
#' Wrapper for one_run()
run_simple_simulation <- function( dgp = "PUMP", 
                                   one_run_func = NULL,
                                   ICC.2 = 0.4, 
                                   omega.2 = 1,
                                   add_canonical_weights = TRUE,
                                   blocks = TRUE,
                                   het_block = TRUE,
                                   R = 1000,
                                   parallel = TRUE,
                                   calc_performance = TRUE ) {
    
    mpl = NULL
    if ( is.character(dgp) && dgp == "PUMP" ) {
        mpl = model.params.list
        mpl$ICC.2 = ICC.2
        mpl$omega.2 = omega.2
    }
    
    if ( parallel ) {
        setup_parallel_if_needed()
    }
    
    if ( is.null( one_run_func ) ) {
        one_run_func = \(r) {
            one_run_func( dgp = dgp,
                          blocks = blocks,
                          het_block = het_block,
                          add_canonical_weights = add_canonical_weights,
                          model.params = mpl )
        }
    }
    
    tictoc::tic()
    rps = future_map_dfr( 1:R, one_run_func,
                          .id = "runID", .progress=TRUE,
                          .options = furrr_options( seed=TRUE ) )
    tictoc::tic()
    
    #head( rps )
    #sd( rps$tau_c)
    
    # Calculate performance metrics and whatnot.
    if ( calc_performance ) {
        calc_performance( rps )
        
        res
    } else {
        rps 
    }
}


if ( FALSE) {
    
    
    my_dgp <- function( ) {
        make_block_data_homog( n=c(2,4), J = c(2,1), K = 2 )
    }
    my_dgp()
    one_run( dgp=my_dgp )
    
    run_simple_simulation( dgp=my_dgp, parallel = FALSE, R = 3 )
}



