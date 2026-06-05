
# explore_person_bias.R
#
# Investigates bias in person-weighted estimators when cluster size is
# strongly correlated with treatment effect (one large outlier cluster).
#
# Output: kable table of bias by method and weight.
# Status: Runs with N_REPS = 5 as of 2026-06-04

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_block.R" ) )

N_REPS <- 5  # Set low for testing; increase for a full run

# Shared block-level args (cluster params differ between the two runs)
grid_base_args <- list(
  K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
  ATE = 0.2, sd.block = 0.3,        # sd.ATE -> sd.block
  ICC.block = 0.2,                    # ICC -> ICC.block
  cv.block = 0.6,                     # variable.J = TRUE -> cv.block
  cv.nbar = 0.33,                     # variable.nbar = TRUE -> cv.nbar
  variable.p = FALSE,
  num.cov.block = 3,                  # num.cov.district -> num.cov.block
  size.intercept.cor = -0.8,
  size.impact.cor = 0,
  nbar.J.cor = -0.9,
  empirical = TRUE,
  sd.cluster = 0.5,                   # sd.school -> sd.cluster
  ICC.cluster = 0.2,
  cv.cluster = 0.6                    # cv.school -> cv.cluster
)

# Grid with strong size-impact correlation at cluster level
grid_corr <- do.call( build_blocks,
                      c( grid_base_args, list( size.impact.cor.cluster = 0.75 ) ) ) %>%
  arrange( J )
grid_corr$J[1]        <- 4
grid_corr$n.bar[1]    <- 100
grid_corr$block.tau[1] <- -2          # district.tau -> block.tau
cor( grid_corr$J, grid_corr$n.bar )
describe_block_grid( grid_corr )      # describe_district_grid -> describe_block_grid

# Grid with no correlation
grid_nocorr <- do.call( build_blocks,
                        c( grid_base_args, list( size.impact.cor.cluster = 0 ) ) ) %>%
  arrange( J )
grid_nocorr$J[1]        <- 4
grid_nocorr$n.bar[1]    <- 100
grid_nocorr$block.tau[1] <- -2

# Describe data structure for correlated grid
ac <- single.CRT.fixed.block( grid = grid_corr, n.runs = N_REPS,
                               just.describe.data = TRUE )
ac

# Run simulation with no correlation
a <- single.CRT.fixed.block( grid = grid_nocorr, n.runs = N_REPS,
                              .progress = TRUE )
table( a$subrun )
a

a %>%
  filter( disfavored != 1 ) %>%
  group_by( method, weight ) %>%
  summarise( SE.b     = sd( ATE.hat ) / sqrt( n() ),
             bias.stu = mean( ATE.hat - ATE.finite.student ),
             bias.sch = mean( ATE.hat - ATE.finite.cluster ),  # ATE.finite.school -> ATE.finite.cluster
             .groups = "drop" ) %>%
  arrange( weight ) %>%
  knitr::kable( digits = 2 )
