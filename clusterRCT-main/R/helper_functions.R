##
## Utility functions for the rest of the package
##



scat = function( str, ... ) {
    dots = list( ... )
    nulls = sapply( dots, is.null )
    if ( any( nulls ) ) {

        # Convert list to string with names and values
        result <- paste(1:length(dots), dots, sep = ": ")

        # Concatenate all elements into a single string
        final_string <- paste(result, collapse = "\n")

        # Print the final string
        warning( paste0( "In scat, null values: ", final_string ) )
    }
    cat( sprintf( str, ... ) )
}



# To go with glue, e.g., catn( glue::glue( "blah" ) )
#
# (glue strips newlines)
catn <- function( str ) {
    str = paste0( str, "\n" )
    cat( str )
}



two_sided_p <- function( ATE_hat, SE_hat, df ) {
    2 * pt( -abs(ATE_hat)/SE_hat, df=df )
}




#### Expand control variables to dummies as needed ####



#' Take data with possibly categorical control variables and expand to
#' dummy variables.
#'
#' @return List of three objects: the revised data frame, revised
#'   control_formula, and final list of names of control variables.
expand_control_variables <- function( data, control_formula ) {
    require( formula.tools )

    n = nrow(data)
    # Get control variables and expand control matrix with dummy variables, etc.
    # This preserves rows with missing values.  Weird code, but is how it has to be.
    controls = model.matrix(control_formula, model.frame(control_formula, data, na.action=na.pass))

    stopifnot( nrow(controls) == n )
    stopifnot( colnames(controls)[1] == "(Intercept)" )
    controls = controls[,-1, drop=FALSE]

    c_names = make.names( colnames(controls), unique=TRUE )
    colnames(controls) = c_names

    # Update dataframe to expanded controls
    data = data %>%
        dplyr::select( -all_of( formula.tools::rhs.vars(control_formula) ) ) %>%
        cbind( as.data.frame( controls ) )

    control_formula = as.formula( paste0( "~ ", paste( c_names, collapse = " + " ) ) )

    return( list( data = data,
                  control_formula = control_formula,
                  c_names = c_names ) )
}


# Calculate number of variables in the control side.
number_controls <- function( control_formula ) {

    if ( is.null( control_formula ) ) {
        return( 0 )
    }

    if (length(formula.tools::lhs.vars(control_formula)) != 0 |
        length(formula.tools::rhs.vars(control_formula)) < 1) {
        stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
    }

    c_names <- formula.tools::rhs.vars(control_formula)
    return( length( c_names ) )
}




#### Aggregate data to clusters ####


#' Aggregate data to the cluster level.
#'
#' Also implements checks to ensure that all clusters are entirely
#' treated or not treated.
#'
#' If covariates specified, categorical covariates will be converted
#' to dummy variables and then averaged (even if they are only level-2
#' or level-3 covariates)
#'
#' @param data  Dataframe with 'clusterID' and 'Z' as columns.
#'   'blockID' optional column.  This is data in the "canonical form".
#'
#' @return tibble of cluster-aggregated data, including Ybar, n,
#'   blockID, clusterID, and Z
#'
#' @export
aggregate_data <- function( data, control_formula = NULL ) {
    # Aggregation to the cluster level

    datagg = NA
    stopifnot( all( c( "clusterID", "Z", "Yobs" ) %in% names(data) ) )

    if ( !is.null( control_formula ) ) {

        dd = expand_control_variables( data, control_formula )
        data = dd$data

        my.vars = dd$c_names
        datagg <-
            datagg <- data %>%
            group_by( across( any_of( c( "blockID", "clusterID", "Z" ) ) ) ) %>%
            summarise( Ybar = mean( Yobs ),
                       n = n(),
                       across( all_of( my.vars ), mean ),
                       .groups = "drop" )

        attr( datagg, "control_formula" ) <- dd$control_formula
        attr( datagg, "c_names" ) <- my.vars

    } else {
        datagg <-
            datagg <- data %>%
            group_by( across( any_of( c( "blockID", "clusterID", "Z" ) ) ) ) %>%
            summarise( Ybar = mean( Yobs ),
                       n = n(),
                       .groups = "drop" )

    }

    J = nrow( datagg )
    stopifnot( J == length(unique(data$clusterID) ) )

    return( datagg )
}






#' Make a canonical regression formula, possibly with control
#' variables.
#'
#' @param Yobs Name of outcome variable (assumed to exist in data)
#' @param Z Name of treatment variable (assumed to exist in data)
#' @param clusterID Name of cluster ID variable (assumed to exist in
#'   data)
#' @param blockID Name of block ID variable (assumed to exist in data,
#'   if this is not null).
#' @param control_formula What variables to control for, in the form
#'   of "~ X1 + X2".
#' @param data Dataframe holding all variables to be used in formula.
#' @param interacted  TRUE means include treatment by block
#'   interactions.  Will override FE flag and set FE to true.
#' @param FE TRUE means include block dummy variables.
#' @param cluster_RE Add a random effect term for clusterID
#'
#' @return Something like "Yobs ~ 1 + Z" or "Yobs ~ 1 + Z + X1 + X2"
#'
#' @export
make_regression_formula = function( Yobs = "Yobs", Z = "Z",
                                    clusterID = "clusterID",
                                    blockID = "blockID",
                                    control_formula = NULL,
                                    interacted = FALSE,
                                    FE = FALSE,
                                    cluster_RE = FALSE,
                                    data = NULL ) {

    require( formula.tools )

    if ( interacted ) {
        FE = TRUE
    }

    # Formula types
    # Yobs ~ 1 + Z + X
    # Yobs ~ 0 + blockID + Z + X
    # Yobs ~ 0 + blockID + Z:blockID + X

    stem = "1"
    if ( FE ) {
        stem = "0"
    }

    txS = "Z"
    if ( interacted ) {
        txS = paste0( Z, ":", blockID )
    }

    if ( FE || interacted ) {
        txS = paste0( txS, " + ", blockID )
    }

    if( !is.null( control_formula ) ) {
        if (length(formula.tools::lhs.vars(control_formula)) != 0 |
            length(formula.tools::rhs.vars(control_formula)) < 1) {
            stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
        }

        if (!is.null(data)) {
            control.vars <- formula.tools::get.vars(control_formula, data = data)
            if (any(!(control.vars %in% colnames(data)))) {
                stop("Some variables in control_formula are not present in your data.")
            }
        }

        c_names <- formula.tools::rhs.vars(control_formula)
        c_names <- paste( c_names, collapse =" + " )

        new.form <- sprintf( "%s ~ %s + %s + %s", Yobs, stem, txS, c_names)
    } else {
        new.form <- sprintf( "%s ~ %s + %s", Yobs, stem, txS )
    }

    if ( cluster_RE ) {
        new.form = paste0( new.form, " + (1|", clusterID, ")" )
    }
    return( as.formula( new.form ) )

}



#### Make canconical data and process variable name arguments ####


#' Process formula to get variable names
#'
#'
deconstruct_var_formula <- function( formula, data ) {
    formula <- as.formula(formula)
    rhs <- formula.tools::rhs.vars( formula )
    parts <- strsplit( rhs, split=" \\| " )[[1]]
    n_part <- length( parts )
    outcome <- formula.tools::lhs.vars(formula)
    if ( (  n_part != 2 && n_part != 3 ) ||
         (length(outcome) != 1 || length(rhs) != 1) ) {
        stop("The formula argument must be of the form outcome ~ treatment | clusterID or outcome ~ treatment | clusterID | blockID.", call. = FALSE )
    }

    clusterID = parts[[2]]
    blockID = NULL
    if ( n_part == 3 ) {
        blockID = parts[[3]]
    }
    main.vars = c( outcome, parts )
    if (!is.null( data ) && any(!(main.vars %in% colnames(data)))) {
        missing_vars = setdiff( main.vars, colnames(data) )
        missing_vars = paste0( missing_vars, collapse = ", " )
        stop( glue::glue( "The following outcome, treatment variables, or cluster variables in formula are not present in your data: {missing_vars}") )
    }



    list( outcome = outcome,
          Z = parts[[1]],
          clusterID = clusterID,
          blockID = blockID )
}

#' Get variables listed in control_formula formula.
#'
#' @return list of variable names
deconstruct_control_formula <- function( control_formula, data ) {

    if(length(formula.tools::lhs.vars(control_formula)) != 0 |
       length(formula.tools::rhs.vars(control_formula)) < 1) {
        stop("The control_formula argument must be of the form ~ X1 + X2 + ... + XN. (nothing on left hand side of ~)")
    }
    control.vars <- formula.tools::get.vars(control_formula, data = data)
    if (any(!(control.vars %in% colnames(data)))) {
        missing_vars = setdiff( control.vars, colnames(data) )
        missing_vars = paste0( missing_vars, collapse = ", " )
        stop( glue::glue( "Following variables in control_formula are not present in your data: {missing_vars}") )
    }

    if ( any( control.vars %in% c( "Y", "Z", "clusterID", "blockID" ) ) ) {
        stop( "Control variables with canonical names of Y, Z, clusterID, or blockID are not allowed." )
    }

    return( control.vars )
}



#' Make canonical data for analysis
#'
#' Given a formula, assemble dataframe of outcome, treatment, and any
#' clustering variables all with canonical names (Y, Z, clusterID,
#' blockID)
#'
#' Some notes on weird behavior to help implement default arguments in
#' rest of package.
#'
#' 1) If formula is NULL, then assume data is already in canonical
#' form, and return it.
#'
#' 2) If formula is not NULL, and data is NULL, then assume data is
#' stored in formula, and is in canonical form, and return it (after
#' verifying it is a data.frame)
#'
#' @param formula Notation for Y ~ Z | clusterID | blockID (| blockID
#'   is optional).
#' @param data Dataframe to pull data from.
#' @param patch_data If TRUE patch data using patch_data_set() after
#'   converting to canonical form.  Overrides drop_missing flag.
#'
#' @return canonical dataframe.
#'
#' @return Dataset with now-canonical variable names, and no
#'   extraneous variables.  Canonical datanames are (Y, Z, clusterID,
#'   blockID).
#'
#' @noRd
make_canonical_data <- function(formula, control_formula = NULL, data,
                                give_default_block = FALSE,
                                drop_missing = TRUE,
                                warn_missing = TRUE,
                                patch_data = FALSE ) {

    if ( patch_data ) {
        new_dat <- patch_data_set(formula, data, control_formula, warn_missing=warn_missing)
        if ( give_default_block && is.null( new_dat$blockID ) ) {
            new_dat$blockID = "Sall"
        }
        return( new_dat )
    }

    # If formula is NULL, or is a dataset, then assume we already are
    # canonical and just return.
    if ( !exists( "data" ) || is.null( data ) ) {
        if ( is.data.frame(formula) ) {
            return( formula )
        }
    } else if ( is.null( formula ) ) {
        return( data )
    }

    parts = deconstruct_var_formula(formula, data)

    Y = data[[ parts$outcome ]] # model.part( formula, data=data, lhs=1, drop = TRUE)
    Z = data[[ parts$Z ]] # model.part( formula, data=data, rhs=1, drop = TRUE)
    clusterID = data[[ parts$clusterID ]] # model.part( formula, data=data, rhs=2, drop=TRUE )
    clusterID = as.factor( clusterID )

    new_dat = data.frame( Yobs = Y, Z = Z, clusterID = clusterID )
    if ( !is.null( parts$blockID ) ) {
        blockID = data[[ parts$blockID ]]
        new_dat$clusterID = interaction(blockID, clusterID, drop=TRUE)
        new_dat$blockID = as.factor( blockID )
    }

    if ( !is.null( control_formula ) ) {
        control.vars = deconstruct_control_formula(control_formula, data)
        xd = data %>% dplyr::select( dplyr::all_of( control.vars ) )
        new_dat = bind_cols( new_dat, xd )
    }

    nr = nrow( new_dat )
    if ( sum( is.na( new_dat ) ) > 0 ) {
        if ( drop_missing ) {
            new_dat = na.omit(new_dat)
            new_dat <- new_dat %>%
                mutate( across( where( is.factor ), droplevels ) )
            drp = nr - nrow( new_dat )
            if ( warn_missing ) {
                warning( glue::glue( "{drp} rows with missing values dropped (of {nr} total rows)." ),
                         call. = FALSE )
            }
        } else {
            if ( warn_missing ) {
                warning( "Data has missing values.", call. = FALSE )
            }
        }
    }

    # Check for valid treatment variable
    n_levels = length( unique( new_dat$Z[ !is.na( new_dat$Z ) ] ) )
    if ( n_levels == 1 ) {
        warning( glue::glue( "Only single level in treatment variable {parts[[1]]}" ) )
    }
    if ( n_levels > 2 ) {
        stop( glue::glue( "Identified treatment variable '{parts$Z}' has more than two values. Did you swap treatment and block?" ))
    }

    if ( give_default_block && is.null( new_dat$blockID ) ) {
        new_dat$blockID = "Sall"
    }

    return( new_dat )
}


#' Any treatment arms with singleton cluster (tx or co)?
#'
#' Used to avoid calculating nonsensical standard errors.
#'
#' @param datagg Aggregated data with blockID and Z columns.
#'
#' @return TRUE or FALSE
has_singleton_clusters_agg <- function( datagg ) {
    tb = table( datagg$blockID, datagg$Z )
    return( any( tb <= 1 ) )
}



#' Identify all tx and all co blocks
#'
#' Identify blocks that are all tx or all co.  Rows with missing
#' treatment indicators are dropped.  Rows with missing block
#' identifiers are considered all unique blocks and thus put into the
#' table of results.
#'
#' @return Dataframe of block IDs corresponding to all tx or all co
#'   blocks along with number of units and tx status.  Returns NULL if
#'   nothing.
#'
#' @export
identify_singleton_blocks <- function( formula, data ) {

    parts = deconstruct_var_formula(formula, data)
    blockID = data[[parts$blockID]]

    blockID = as.character(blockID)
    Z = data[[parts$Z]]
    tb = table( Z, blockID, useNA = "always" )
    zros = apply( tb[-3,], 2, min )

    if ( sum( is.na(blockID) & !is.na(Z) ) == 0 ) {
        zros[ length(zros) ] = 1
    } else {
        zros[ length(zros) ] = 0
    }

    if ( sum( zros == 0 ) == 0 ) {
        return( NULL )
    }
    nms = colnames(tb)[ zros == 0 ]

    n = apply( tb[ 1:2, zros == 0, drop=FALSE ], 2, sum )
    n_tx = tb[ 2, zros == 0 ]
    ptx = n_tx / n
    tb <- tibble( blockID = nms,
            pTx = as.numeric( ptx ),
            n = as.numeric( n ) )
    tb
}



# # #
# if ( FALSE ) {

# make_base_formula()
# make_base_formula( control_formula = ~ X1 + X2 + X3 )
# make_FE_formula(  )
# make_FE_formula( control_formula = ~ X1 + X2 + X3 )

# make_FE_int_formula( )
# make_FE_int_formula( control_formula = ~ X1 + X2 + X3 )

# }




#### Functions to check data integrety ####


#' Check if one grouping variable nested in another
#'
#' Checks if clusterID nested in blockID, so each clusterID value
#' appears in only one blockID value.
#'
#' @param blockID List of blockIDs
#' @param clusterID List of clusterIDs
#'
#' @return TRUE if nested, FALSE otherwise.  If blockID NULL, always
#'   return TRUE.
#' @export
is_nested <- function( clusterID, blockID ) {
    if ( is.factor(clusterID) ) {
        clusterID = droplevels(clusterID)
    }
    if ( is.factor(blockID) ) {
        blockID = droplevels(blockID)
    }

    if ( is.null( blockID ) ) {
        return( TRUE )
    }

    blk = tapply( blockID, clusterID, function( x ) { length( unique( x ) ) } )
    return ( all( blk == 1 ) )
}




#' Check data integrity
#'
#' Ensure data has proper treatment variable, clusterIDs are defined,
#' and so forth.
check_data_integrity <- function( formula = NULL, data ) {

    if ( !is.null(formula) ) {
        if ( is.data.frame(formula) ) {
            data = formula
        } else {
            data = make_canonical_data( formula=formula,
                                        data=data )
        }
    }
    if ( is.null( data$blockID) ) {
        data$blockID = "single_block"
    }

    if ( is.null( data$Yobs ) ) {
        stop( "No outcome", call.=FALSE )
    }
    if ( is.null( data$Z ) ) {
        stop( "No treatment assignment vector", call.=FALSE )
    }
    if ( is.null( data$clusterID ) ) {
        stop( "No clustering ID", call.=FALSE )
    }

    if ( !is.numeric( data$Yobs ) ) {
        stop( "Outcome is not numeric", call.=FALSE )
    }

    if (is.numeric(data$Z) == FALSE ) {
        stop("Treatment indicator should be vector of ones and zeros", call.=FALSE)
    }

    n = nrow( data )
    if ((sum(data$Z == 1) + sum(data$Z == 0)) != n) {
        stop("Treatment indicator should be vector of ones and zeros", call.=FALSE)
    }

    ntx = sum(data$Z)
    if ( ntx == 0 || ntx ==n ) {
        stop( "Can't have all treated or all control across dataset", call.=FALSE )
    }


    sts <- data %>% group_by( blockID ) %>%
        summarise( ptx = mean( Z ) )
    if ( any( sts$ptx == 0 | sts$ptx == 1 ) ) {
        stop( "Some blocks have all treated or all control units", call.=FALSE )
    }


    if ( !is_nested( data$clusterID, data$blockID ) ) {
        stop( "SiteID not fully nested in ClusterID", call.=FALSE )
    }

    ptx = data %>% group_by( clusterID ) %>%
        summarize( ptx = mean( Z ) )
    if ( any( ptx$ptx != 0 & ptx$ptx != 1 ) ) {
        stop( "Treatment variation within clusters: not a cluster RCT", call.=FALSE )
    }

    return( TRUE )
}

