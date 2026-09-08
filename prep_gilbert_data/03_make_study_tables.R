
# Make descriptive tables of the datasets and save to excel sheets.


# load libraries
library(tidyverse)
library(glue)
library(clusterRCT)
library(readxl)
library(writexl)
library(tictoc)

# clear memory
rm(list = ls())

dir.create( here::here( "prep_gilbert_data/tabfig" ), showWarnings = FALSE, recursive = TRUE )

# bring in metadata to get index on which csvs to load
metadata <- read_xlsx( here::here( "prep_gilbert_data/study_metadata.xlsx") ) |>
    filter(cluster_rct == 1)

# we only care about the CRTs here
ids <- metadata |>
    pull(dataset_id)

# function to read in the csv and extract the relevant info
read_in <- function(i){

    # read in
    dat <- read_csv(here::here(glue("prep_gilbert_data/data/raw/crt_data_{i}.csv")))

    # get the number of people
    n <- nrow(dat)

    # number of clusters
    j <- dat |>
        distinct(cluster_id) |>
        nrow()

    # put these in a table
    out <- tibble(
        dataset_id = i,
        n = n,
        j = j
    )

    # add number of blocks, if applicable
    if(!is.null(dat$block_id)){

        k <- dat |>
            distinct(block_id) |>
            nrow()

        out <- out |>
            mutate(k = k)

    }

    # return
    return(out)

}

# map across the CRT ids
tab <- map(ids, read_in) |>
    bind_rows()

# combine with the other bits of the metadata that are relevant
tab <- tab |>
    left_join(metadata, by = "dataset_id") |>
    # group by study, not data set
    group_by(study_id) |>
    mutate(n_outcomes = n(),
           mean_n = mean(n) |> round(),
           mean_j = mean(j) |> round(),
           mean_k = mean(k) |> round(),
           min_n = min(n), max_n = max(n),
           min_j = min(j), max_j = max(j),
           outcome = paste(outcome, collapse = ", "),
           year = parse_number(year)
           ) |>
    distinct(study_id, .keep_all = TRUE) |>
    select(study_id, lead_author, year, treatment,
           intervention_name, n_outcomes, outcome, min_n, max_n, mean_n, sample, min_j, max_j, mean_j,
           cluster_unit, k, block_unit)

# export
write_xlsx(tab, here::here( "prep_gilbert_data/tabfig/study_summaries.xlsx" ))

# make a long version, like appendix table 2

# aggregated by study id
long <- tab |>
    mutate(report = glue("{lead_author} ({year})"),
           program = if_else(intervention_name == "NA", treatment, intervention_name),
           design = if_else(is.na(k), "CRT", "Blocked CRT")) |>
    distinct(study_id, .keep_all = TRUE) |>
    mutate(lab = glue("Study {study_id}: {program}"),
           k = replace_na(k, 0),
           size = glue("Individuals: {mean_n}, Clusters: {mean_j}, Blocks: {k}")) |>
    select(lab, sample, design, outcome, size, report) |>
    pivot_longer(sample:report, names_to = "characteristic", values_to = "value")

write_xlsx(long, here::here( "prep_gilbert_data/tabfig/study_long.xlsx" ))
