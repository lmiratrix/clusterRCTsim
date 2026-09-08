
#
# Load all the datasets and do some mild cleaning.
#


# load libraries
library(tidyverse)
library(glue)
library(clusterRCT)
library(readxl)
library(tictoc)

# clear memory
rm(list = ls())




# bring in metadata to get index on which csvs to load
metadata <- read_xlsx(here::here("prep_gilbert_data/study_metadata.xlsx"))
metadata

# Get all cluster RCTs listed in metadata file
# We will iterate over all of these.
ids <- metadata |>
    filter(cluster_rct == 1) |>
    pull(dataset_id)
ids

if ( FALSE ) {
    # Get IDs by bringing in all the csvs
    ids <- list.files(path = here::here( "prep_gilbert_data/data/raw/" ),
                      pattern = "crt_data_\\d+\\.csv", full.names = TRUE)
    id_id = parse_number(ids)
    id_id = id_id[ order( id_id ) ]
    ids = id_id
    rm( id_id )
}



source( here::here( "prep_gilbert_data/processing_functions.R" ) )



# Read in all the data we will analyze ----

poss_read_in = possibly( read_in, otherwise = NULL )

dfs <- map(ids, poss_read_in)
names( dfs ) <- ids
ids <- as.character(ids)

dfs = dfs[ !map_lgl( dfs, is.null ) ]




# Cleaning specific datasets ----

# WARNING: The rest of the code uses global variables given the list
# of datasets grabbed previously.  Not the best written code.


# some weird misclassification in kim - change treatment to modal
# treatment assignment for school
dfs[["2"]] <- dfs[["2"]] |>
    group_by(cluster_id) |>
    mutate(treat = round(mean(treat), 0)) |>
    ungroup()
describe_clusterRCT( y_std ~ treat | cluster_id | block_id, data = dfs[["2"]] )

# dataset 39 has what appears to be 1 misclassified
dfs[["49"]] <- dfs[["49"]] |>
    # id 78 appears to be in wrong cluster number
    mutate(cluster_id = if_else(s_id == 78, 63, cluster_id))
describe_clusterRCT( y_std ~ treat | cluster_id, data = dfs[["49"]] )

# one cluster with two people, one treat/one control, remove
dfs[["50"]]
get_structure( y_std ~ treat | cluster_id,  dfs[["50"]] ) |>
    arrange( n )
filter( dfs[["50"]], cluster_id == 211115 )
dfs[["50"]] <- dfs[["50"]] |>
    filter(cluster_id != 211115)
describe_clusterRCT( y_std ~ treat | cluster_id, data = dfs[["50"]] )

get_structure( y_std ~ treat | cluster_id,  dfs[["51"]] ) |>
    arrange( n )
dfs[["51"]] <- dfs[["51"]] |>
    # one cluster with two people, one treat/one control, remove
    filter(cluster_id != 211115)


# due to attrition, study 52 has block 1 with only treated units in a
# single cluster.  Filter it out.
blk_1 <- filter( dfs[["52"]], block_id == 1 )
blk_1
dfs[["52"]] <- dfs[["52"]] |>
    filter(block_id != 1)
rm( blk_1 )

# same modal treatment thing with 74 - same source so probably same issue
dfs[["74"]]
struct <- get_structure( y_std ~ treat | cluster_id, data = dfs[["74"]] )
if ( FALSE ) {
    # Most clusters are _almost_ all treated or _almost_ all control
    qplot( struct$Z )
}
dfs[["74"]] <- dfs[["74"]] |>
    group_by(cluster_id) |>
    mutate(treat = round(mean(treat), 0)) |>
    ungroup()
rm( struct )

describe_clusterRCT( y_std ~ treat | cluster_id, data = dfs[["74"]] )



# Drop some missing values
dfs[["16"]] <- dfs[["16"]] |>
    filter( !is.na( block_id ) )


ids = names( dfs )

cat( "Loading complete\n" )
