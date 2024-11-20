

# This file contains functions to generate data for simulations.

library( tidyverse )


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



#### Function to make data via PUMP methods ####


#' @param blocks TRUE means multiple districts (blocks).  FALSE means
#'   single block.
#' @param het_block TRUE means blocks are different shapes and styles
#'
make_two_tier_data <- function( n1, ATE1, n2 = 5*n1, ATE2 = 3*ATE1,
                                blocks = TRUE, params = model.params.list,
                                het_block = FALSE ) {

    design = "d3.2_m3ff2rc"
    if ( het_block && !blocks ) {
        stop( "Invalid--can't have het blocks with no blocks" )
    }

    if ( !blocks ) {
        params$ICC.3 = 0
        params$R2.3 = 0
        params$omega.3 = 0
        design = "d2.2"
        params$J = params$J * params$K
        params$K = NULL
    }

    # Make small blocks across the districts
    params$nbar = n1
    params$MDES = ATE1
    sim.data.small <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.30, include_POs = TRUE )

    # Make large blocks across the districts
    params$nbar = n2
    params$MDES = ATE2
    sim.data <- PUMP::gen_sim_data( d_m = design, params, Tbar = 0.30,
                                    include_POs = TRUE )
    sim.data$S.id = paste0( "b-", sim.data$S.id )

    # Bind so each district has small and large blocks
    sim.data = bind_rows( sim.data.small, sim.data )

    if ( het_block ) {
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

        sim.data <- filter( sim.data, !( S.id %in% c( paste0( "b-", 11:16 ) ) ) )

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

        # Rerandomize units within the blocks
        sim.data$Tx = randomizr::block_and_cluster_ra( blocks = sim.data$D.id,
                                                       clusters = sim.data$S.id,
                                                       block_m = c( 10, 6, 10, 5 ) )
        sim.data <- mutate( sim.data,
                            Yobs = ifelse( Tx == 1, Y1, Y0 ) )
    } else {
        sim.data$Tx = sim.data$T.x
    }

    # Standardize outcome to effect size units
    sdY0 = sd( sim.data$Y0 )
    sim.data <- mutate( sim.data,
                        Yobs = Yobs / sdY0,
                        Y0 = Y0 / sdY0,
                        Y1 = Y1 / sdY0 )
    sim.data$T.x = NULL

    as_tibble( sim.data )
}




##### Testing make_two_tier_data #####
if ( FALSE ) {
    dat = make_two_tier_data( 10, 0.2, het_block = TRUE)
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
    dat = make_two_tier_data( 10, 0.2, het_block = TRUE)
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




    dat = make_two_tier_data( 10, 0.2, blocks=FALSE )
    dat
    dat$D.id
    table( table( dat$S.id ) )
    sd( dat$Y0 )
}



#### Function to make data by specifying mixture of blocks and cluster sizes ####

#' Make data with very different block structure and cluster mixtures
#' that correlate with treatment effects.
#'
#' Each block will be a mix of clusters of specified sizes, with
#' specified number of clusters of each size.  The clusters will each
#' have a mean outcome for that cluster, and a cluster-specific
#' treatment effect.  I.e., blocks are a mixture of a discrete number
#' of cluster types.
#'
#' "Randomizations" will be perfectly balanced to focus on bias, not
#' estimation error.
#'
#' @param n List lists of cluster sizes within each block
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

        statblock = list( n, J, Y0, ATE, ptx )

    }

    # Make J clusters of size n with an ATE of given ATE.
    make_blocks <- function( n, J, Y0, ATE, ptx = 0.5, prefix ) {
        dat = tibble( S.id = rep( 1:J, each=n ),
                      Y0 = Y0 + rnorm( n*J, sd=1 ) )
        dat$Y1 = dat$Y0 + ATE
        dat$Tx = 0 + (dat$S.id %% (1/ptx) == 0)
        dat$Yobs = ifelse( dat$Tx, dat$Y1, dat$Y0 )
        dat$S.id = paste0( prefix, "-", dat$S.id )
        dat
    }

    # Given vectors of n, J, Y0, and ATE, make batches of clusters of
    # given sizes via pmap
    make_blocks_v <- function( n, J, Y0, ATE, ptx=0.5 ) {
        stopifnot( length(n) == length(J) )
        stopifnot( length(n) == length(Y0) )
        stopifnot( length(n) == length(ATE) )
        if ( length(ptx) == 1 ) {
            ptx = rep( ptx, length(J) )
        }
        stopifnot( length(n) == length(ptx) )
        if ( length( Y0 ) == 1 ) {
            Y0 = rep( Y0, length(n) )
        }

        res <- pmap_df( list( n=n, J=J, Y0=Y0, ATE=ATE, ptx=ptx,
                              prefix=letters[ 1:length(n) ] ),
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



if ( FALSE ) {

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
}







#' Make data with very different block structure and cluster mixtures
#' that correlate with treatment effects.
#'
#' "Randomizations" will be perfectly balanced (assuming the
#' proportion of treatment evenly divides each block group size) to
#' focus on bias, not estimation error.
#'
#' @param blocks TRUE means multiple districts (blocks)
#' @param het_block TRUE means blocks are different shapes and styles
#'
make_simple_block_example <- function( version = "A", tx_scale = 1, var_tx = FALSE ) {

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



