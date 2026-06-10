
# run_single_simulation.R
#
# Runs a mini simulation for a single scenario from the fullfactor grid.
# Useful for debugging or spot-checking a specific scenario.
#
# Output: results/single_ID.rds — raw replicate-level results from
#   run.scenario.fixed.block(), one row per (run x method).
#
# Usage: set SCENARIO_ID below, then source() this file from the
#   project root (or run via Rscript from the project root).


rm(list = ls())

library(here)
library(tidyverse)
library(purrr)

source(here::here("R/data_generation_code.R"))
source(here::here("R/simulation_support_code_fixed_block.R"))


# ---- Configuration ----

SCENARIO_ID <- 439       # raw ID in scenarios_fullfactor.rds (= "439_fullfactor" in sruns.clean)
NUM_REPS    <- 100        # number of synthetic datasets to generate
OUTPUT_FILE <- here::here( glue::glue( "results/single_{SCENARIO_ID}.rds") )


# ---- Load scenario grid ----

scenarios <- readRDS(here::here("results/scenarios.rds"))

table( scenarios$cv.cluster, scenarios$blocktype )

scenario  <- filter(scenarios, ID == SCENARIO_ID)
stopifnot("Scenario ID not found in scenarios_fullfactor.rds" = nrow(scenario) == 1)

cat(sprintf("\nScenario ID %d  (%d replicates)\n", SCENARIO_ID, NUM_REPS))

cat("DGP parameters:\n")
scenario[, !names(scenario) %in% "grid"] %>%
  t() %>%
  knitr::kable()

cat("\nBlock-level grid:\n")
grid <- scenario$grid[[1]]
grid %>%
  t() %>%
  knitr::kable( digits=2 )


describe_block_grid( grid )

# ---- Run simulation ----

set.seed(SCENARIO_ID^2 + SCENARIO_ID + 1)

results <- run.scenario.fixed.block(
  grid     = grid,
  R        = NUM_REPS,
  ID       = SCENARIO_ID,
  parallel = FALSE
)

cat(sprintf("\nResults: %d rows x %d columns\n", nrow(results), ncol(results)))
cat("Methods present:", paste(sort(unique(results$method)), collapse = ", "), "\n")


source( here::here( "R/group_estimators.R" ) )
results <- results %>%
  mutate( runID = paste0( run, subrun ) )
g <- group_estimators( results, var = "ATE.hat" )
g
core = sapply( g, function(x) sort(x)[[1]] )
core

cat("\nSE and SE.hat by method (across replicates):\n")
sts <- results %>%
  dplyr::filter( method %in% core ) %>%
  group_by(method) %>%
  summarise(trueSE = sd( ATE.hat ),
            mean_SE = mean(SE.hat, na.rm = TRUE),
            sd_SE   = sd(SE.hat, na.rm = TRUE),
            .groups = "drop") %>%
  mutate( relSE = trueSE / median(trueSE ) ) 

sts %>%
  arrange(mean_SE) %>%
  knitr::kable( digits=2 )

max( sts$trueSE ) / min( sts$trueSE )


# ---- Save ----

saveRDS(results, file = OUTPUT_FILE)
cat(sprintf("\nSaved to %s\n", OUTPUT_FILE))
