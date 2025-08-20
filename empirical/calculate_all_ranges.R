

# Calculate all the ranges for the various configurations of estimators





# Function to calculate both ATE and SE summary statistics at once
calculate_summary_statistics <- function(data,
                                         ate_col = "ATE_hat_std", 
                                         se_col = "SE_hat_std", 
                                         ratio_cap = 5) {
  # Calculate ATE statistics
  ate_stats <- data %>%
    summarize(
      ATE_N     = sum(!is.na(.data[[ate_col]])),
      ATE_Nmiss = sum(is.na(.data[[ate_col]])),
      ATE_Min   = min(.data[[ate_col]], na.rm = TRUE),
      ATE_Max   = max(.data[[ate_col]], na.rm = TRUE),
      ATE_avg   = mean(.data[[ate_col]], na.rm = TRUE),
      SE_N      = sum(!is.na(.data[[se_col]])),
      SE_Nmiss  = sum(is.na(.data[[se_col]])),
      SE_Min    = min(.data[[se_col]], na.rm = TRUE),
      SE_Max    = max(.data[[se_col]], na.rm = TRUE),
      SE_avg    = mean(.data[[se_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Range = ATE_Max - ATE_Min,
      Ratio = pmin(ratio_cap, SE_Max / SE_Min)
    )
  
  return(ate_stats)
}


# Define common grouping variables
base_groups <- c("UniqueID", "class", "low_df", "has_single", "Collection", 
                 "Covariates")
estimand_groups <- c(base_groups, "Estimand")

# Get summary statistics grouped by base variables
summary_stats <- combined %>%
  group_by( across( all_of( base_groups ) ) ) %>%
  calculate_summary_statistics()
summary_stats



# Get summary statistics grouped by estimand
summary_stats_estimand <- combined %>%
  group_by( across( all_of( c( base_groups, "Estimand" ) ) ) ) %>%
  calculate_summary_statistics()
summary_stats_estimand


# Filter out fully interacted models
combined_exclFI <- 
  
  # Get summary statistics excluding fully interacted models
  summary_stats_exclFI <- combined %>%
  filter(!grepl('-FI', Estimator)) %>%
  group_by( across( all_of( base_groups ) ) ) %>%
  calculate_summary_statistics()


# Get summary statistics by estimand excluding fully interacted models
summary_stats_estimand_exclFI <- combined %>%
  filter(!grepl('-FI', Estimator)) %>%
  group_by( across( all_of( c( base_groups, "Estimand" ) ) ) ) %>%
  calculate_summary_statistics()


