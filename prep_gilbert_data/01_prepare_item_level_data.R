

# Score raw item-level response data (Rasch model) into the
# person-by-outcome crt_data_{id}.csv format the rest of this pipeline
# expects.
#
# NOTE: The input to this script - a .Rdata file of item-level
# datasets exported from the IRT Warehouse - is not included in this
# repo (it lives only on the original processor's machine). This
# script is included for provenance/transparency, not because it can
# be independently rerun here. If you have your own copy of that
# export, set the JG_ITEM_DATA_PATH environment variable to its path.
# Otherwise, start from the already-scored files in data/raw/ (see
# README.md).


# load libraries
library(tidyverse)
library(lme4)
library(broom.mixed)
library(see)
library(haven)
library(glue)
library(mirt)
library(sjlabelled)
library(readxl)
library(janitor)
library(dtplyr)

# set ggplot themes
theme_set(theme_bw())
theme_update(legend.position = "bottom", legend.box = "vertical")

# clear memory
rm(list = ls())

# set seed for reproducibility
set.seed(2024)

# path to raw data
raw <- here::here( "prep_gilbert_data/data/raw" )
clean <- here::here( "prep_gilbert_data/data/clean" )
dir.create( raw, showWarnings = FALSE, recursive = TRUE )
dir.create( clean, showWarnings = FALSE, recursive = TRUE )

# load the item-level data
load( here::here( "prep_gilbert_data/datasets_list.Rdata" ) )

# take response time out of 2 datasets
datasets[[46]] <- datasets[[46]] |> select(-rt)
datasets[[47]] <- datasets[[47]] |> select(-rt)

crt_clean <- function(data, num){

  # if I forgot to remove any missing item responses, do so here
    data <- data |>
        drop_na(score)

  # get the minimum time in the data
  min_time <- min(data$time)

  # is the outcome polytomous?
  dichotomous <- is.null(data$polyscore)

  # if polytomous, use polyscore
  if(dichotomous == FALSE){
    data <- data |>
      mutate(resp = polyscore)
  } else {
    data <- data |>
      mutate(resp = score)
  }

  # if have baseline items,
  # construct baseline score
  # using a rasch model
  if(min_time == 0){

    # pivot to wide
    wide <- data |>
      filter(time == 0) |>
      select(s_id, item, resp) |>
      pivot_wider(names_from = item,
                  values_from = resp,
                  id_cols = s_id)

    # fit the model
    irt_mod <- wide |>
      select(-s_id) |>
      mirt(1, itemtype = "Rasch")

    # get the scores
    wide <- wide |>
      mutate(baseline = fscores(irt_mod)[,],
             std_baseline = scale(baseline)[,]) |>
      select(s_id, std_baseline)

    # add to df
    data <- data |>
      left_join(wide, by = "s_id")

  }

  # put everything together
  wide <- data |>
    # endline items only
    filter(time == 1) |>
    # remove time b/c everything is endline
    select(-time, -contains("score")) |>
    pivot_wider(names_from = item,
                values_from = resp,
                names_prefix = "post_item_") |>
    mutate(dataset_id = num)

  # get 1PL scores for the outcome
  rasch_mod <- wide |>
    select(contains("post_item")) |>
    mirt(1, itemtype = "Rasch", verbose = FALSE)

  wide <- wide |>
    mutate(y = fscores(rasch_mod)[,]) |>
    select(-contains("post_item"))

  # standardize according to control group mean and SD
  control_desc <- wide |>
    filter(treat == 0) |>
    summarise(sd = sd(y, na.rm=TRUE),
              mean = mean(y,na.rm=TRUE))

  wide <- wide |>
    mutate(y_std = (y - control_desc[[1,2]])/(control_desc[[1,1]]))

  # write csv
  write_csv(wide, glue("{raw}/crt_data_{num}.csv"))

}

# get number of datasets
num <- length(datasets)
cli::cli_inform( "Working with {num} datasets" )

# map across to export as csvs
# use walk b/c I just want the output
walk2(datasets, 1:num, crt_clean, .progress = TRUE)


