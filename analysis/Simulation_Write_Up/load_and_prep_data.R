# Load and prep simulation results

 
require(tidyverse)
require(rlang)
require(stargazer)
require( ggthemes )
require(broom)
require(knitr)
require(kableExtra)
require(patchwork)
require(gridExtra)
require(pander)



load( here::here( "results_agg/processed_simulation_results_agg_full_factor_2025-05-21.RData" ) )



DGP_vars <- c("K", "J.bar", "p.bar", "sd.ATE", "ICC.district", "variable.p", 
              "n.bar", "sd.school", "ICC.school", "size.impact.cor", "size.impact.cor.school", 
              "proptx.impact.cor", "num.cov.district", "num.cov.school", "cv.school", 
              "num.cov.indiv", "empirical", "cv.district", "nbar.J.cor", 
              "size.intercept.cor", "cv.nbar", "R")


# Function to create labeled factors with zero as base level
create_labeled_factor <- function(x) {
  unique_values <- sort(unique(x))
  num_levels <- length(unique_values)
  
  # Check if 0 is among the values
  has_zero <- 0 %in% unique_values
  
  if (num_levels == 1) {
    # If only one level, just return as is
    return(factor(x))
  } else if (num_levels == 2) {
    if (has_zero) {
      # If one of the two values is 0, use "0" and "+"
      labels <- ifelse(unique_values == 0, "0", "+")
      return(factor(x, levels = unique_values, labels = labels))
    } else {
      # If no zero, use "S" and "L"
      labels <- c("S", "L")
      return(factor(x, levels = unique_values, labels = labels))
    }
  } else if (num_levels == 3) {
    if (has_zero) {
      # If one of the three values is 0, use "0", "S", and "L"
      zero_pos <- which(unique_values == 0)
      ordered_levels <- c(0, unique_values[unique_values != 0])
      labels <- c("0", "S", "L")
      return(factor(x, levels = ordered_levels, labels = labels))
    } else {
      # If no zero, use "S", "M", and "L"
      labels <- c("S", "M", "L")
      return(factor(x, levels = unique_values, labels = labels))
    }
  } else {
    # For more than three levels
    if (has_zero) {
      # Reorder so 0 is first
      ordered_levels <- c(0, unique_values[unique_values != 0])
      # First label is "0", rest are numbered L1, L2, etc.
      labels <- c("0", paste0("L", 1:(num_levels-1)))
      return(factor(x, levels = ordered_levels, labels = labels))
    } else {
      # If no zero value exists, just use S1, S2, etc.
      labels <- paste0("L", 1:num_levels)
      return(factor(x, levels = unique_values, labels = labels))
    }
  }
}


# Drop impossible combos ---

# Check combos of simulations
table( cv.school = sruns$cv.school, 
       sd.school = sruns$sd.school,
       size.impact.cor.school = sruns$size.impact.cor.school )


table( sruns$variable.p, sruns$proptx.impact.cor )

table( sruns$size.impact.cor, sruns$sd.school )


sruns <- sruns %>%
  filter( !(size.impact.cor.school != 0 & sd.school == 0) ) %>%
  filter( !(size.impact.cor.school != 0 & cv.school == 0) ) 







# Make clean set ----



sruns.clean <- sruns %>%  
  mutate(block_type = case_when(
    # Blocks all the same (btype1)
    cv.district == 0 &
      sd.ATE == 0 & 
      size.impact.cor == 0 & 
      nbar.J.cor == 0 &
      proptx.impact.cor == 0 & 
      size.intercept.cor == 0 &
      !variable.p &
      cv.nbar == 0 ~ 1,
    
    # Blocks all messy (btype2)
    cv.district == 0.85 &
      sd.ATE == 0.6 & 
      size.impact.cor == 0.8 & 
      nbar.J.cor == 0.8 &
      proptx.impact.cor == 0.8 & 
      size.intercept.cor == 0.8 &
      variable.p &
      cv.nbar == 0.33 ~ 3,
    
    # Blocks vaguely reasonable (btype3)
    cv.district == 0.33 &
      sd.ATE == 0.2 & 
      size.impact.cor == 0.2 & 
      nbar.J.cor == 0.3 &
      proptx.impact.cor == 0 & 
      size.intercept.cor == 0.3 &
      !variable.p &
      cv.nbar == 0.33 ~ 2,
    
    TRUE ~ NA_real_  # Catch any other combinations
  )) %>% 
  mutate(
    bias.s.estimand = case_when(
      weight == "Cluster" ~ bias.s.sch,
      weight == "Person" ~ bias.s.stu
    ),
    cover.s.estimand = case_when(
      weight == "Cluster" ~ cover.s.sch,
      weight == "Person" ~ cover.s.stu
    ),
    RMSE.s.estimand = case_when(
      weight == "Cluster" ~ RMSE.s.sch,
      weight == "Person" ~ RMSE.s.stu
    )
  ) %>%
  mutate(
    sd.ATE = create_labeled_factor(sd.ATE),
    size.impact.cor = create_labeled_factor(size.impact.cor),
    J.bar = create_labeled_factor(J.bar),
    p.bar = create_labeled_factor(p.bar),
    sd.school = create_labeled_factor(sd.school),
    ICC.school = create_labeled_factor(ICC.school),
    cv.school = create_labeled_factor(cv.school),
    num.cov.school = create_labeled_factor(num.cov.school),
    proptx.impact.cor = create_labeled_factor(proptx.impact.cor),
    size.intercept.cor = create_labeled_factor(size.intercept.cor),
    size.impact.cor.school = create_labeled_factor(size.impact.cor.school),
    # Fix for the factor problem - explicitly define levels based on your block types
    block_type = factor(block_type, 
                        levels = c(1, 2, 3),
                        labels = c("Same", "Reas", "Mess")),
    bias.s.estimand.mag = abs(bias.s.estimand),
    method = reorder(method, -biased),
    bias.ratio.sch = round(bias.f.sch/bias.s.sch, 4),
    bias.ratio.stu = round(bias.f.stu/bias.s.stu, 4),
    RMSE.ratio.sch = round(RMSE.f.sch/RMSE.s.sch, 4),
    RMSE.ratio.stu = round(RMSE.f.stu/RMSE.s.stu, 4),
    cover.ratio.stu = round(cover.f.stu/cover.s.stu, 4),
    cover.ratio.sch = round(cover.f.sch/cover.s.sch, 4),    
    true.se.ratio = round(SE.f/SE.s, 2)
  )  %>%
  mutate(
    is_FI = str_detect(method, "FI"),
    is_MLM = str_detect(method, "MLM"),
    is_cw = str_detect(method, "cw"),
    is_pw = str_detect(method, "pw"),
  ) 


# Check the levels of one of the variables
levels(sruns.clean$block_type)
table(sruns.clean$block_type)
levels(sruns.clean$J.bar)
levels(sruns.clean$ICC.school)

# Create the new variable with median RMSE.s.estimand for each ID
sruns.clean <- sruns.clean %>%
  group_by(ID) %>%
  mutate(RMSE.s.estimand.median = median(RMSE.s.estimand, na.rm = TRUE),
         SE.s.median = median(SE.s, na.rm = TRUE),
         CI.length.median = median(CI.length, na.rm = TRUE),
         sd.SE.hat.s.median=median(sd.SE.hat.s,na.rm=TRUE)) %>%
  ungroup()


# Create ordered factor levels
ordered_methods <- sruns.clean %>%
  # First get Cluster unbiased methods
  filter(weight == "Cluster" & biased == 0) %>% pull(method) %>% unique() %>%
  c(
    # Then Cluster biased methods
    sruns.clean %>% filter(weight == "Cluster" & biased == 1) %>% pull(method) %>% unique(),
    # Then Person unbiased methods
    sruns.clean %>% filter(weight == "Person" & biased == 0) %>% pull(method) %>% unique(),
    # Finally Person biased methods
    sruns.clean %>% filter(weight == "Person" & biased == 1) %>% pull(method) %>% unique()
  )

# Create the plot data with ordered methods
# 
# TODO: Just make this sruns.clean everywhere???
plot_data <- sruns.clean %>%
  mutate(
    weight_biased = paste(weight, biased, sep = "_"),
    method_ordered = factor(method, levels = ordered_methods)
  )


#figure out which estimators produce the same ATE results 

group_estimators <- function( rps, var ) {
  
  wid <- rps %>% 
    dplyr::select( ID, method, all_of(var) ) %>% 
    pivot_wider( names_from = method, values_from = all_of( var ) ) %>%
    dplyr::select( -ID )
  
  results <- combn(names(wid), 2, function(x) {
    mean_diff <- mean(wid[[x[1]]] - wid[[x[2]]])
    tibble(pair = paste(x, collapse = " - "), 
           mean_difference = round( mean_diff, digits=10 ) )
  }, simplify = FALSE) %>% bind_rows()
  
  
  # Create a matrix where each element is the mean difference between two columns
  mean_diff_matrix <- outer(1:ncol(wid), 1:ncol(wid), Vectorize(function(x, y) {
    abs( mean(wid[[x]] - wid[[y]]) )
  }))
  
  # Naming rows and columns for clarity
  colnames(mean_diff_matrix) <- names(wid)
  rownames(mean_diff_matrix) <- names(wid)
  
  mean_diff_matrix
  
  # Assuming mean_diff_matrix is your matrix of pairwise mean differences
  tolerance <- 10^-7
  
  # Create a binary matrix where 1 indicates differences below the tolerance
  binary_matrix <- 0 + (mean_diff_matrix < tolerance)
  binary_matrix
  
  # Convert this binary matrix into a list of groups
  groups <- list()
  already_included <- c()
  
  ests = sort( rownames(binary_matrix) )
  for (i in ests) {
    if (!( i %in% already_included)) {
      # Find indices where the estimator is within tolerance of others
      in_group <- which(binary_matrix[i, ] == 1)
      group_names <- colnames(binary_matrix)[in_group]
      groups[[length(groups) + 1]] <- group_names
      already_included <- c(already_included, group_names)
    }
  }
  groups
}

groups <- group_estimators(sruns, var = "mean.ATE.hat")


# List the groups
map_chr(groups, paste, collapse = ", ") %>%
  # Convert to a data frame with a proper column name
  data.frame(group_values = .) %>%
  # Now n() will work inside mutate()
  mutate(grp = 1:nrow(.)) %>%
  relocate(grp) %>%
  knitr::kable()

#this is the list of methods that don't produce the same ATE estimates
core.methods <- c("LRa-FE-db", "LRa-FIcw-db", "LRapw-FE-db", "LRapw-FIpw-db", "MLM-FE", "MLM-FIcw", "MLM-RE")


#this is the list of MLM methods
MLM.methods<-c( "MLM-FE", "MLM-FIcw", "MLM-RE" )

#this is the list of methods that don't produce the same se estimates
all.methods<-c("LRa-FE-db", "LRa-FE-het", "LRicw-FE-crve", "LRicw-FE-db",
               "LRa-FIcw-db","LRa-FIcw-het", "LRicw-FIcw-crve", 
               "LRicw-FIcw-db", "LRapw-FE-db", "LRapw-FE-het", 
               "LRi-FE-crve", "LRi-FE-db","LRapw-FIpw-db", 
               "LRapw-FIpw-het", "LRi-FIpw-crve", "LRi-FIpw-db", 
               "MLM-FE", "MLM-FIcw","MLM-RE")


groups <- group_estimators(sruns, var = "mean.SE.hat")


# List the groups
map_chr(groups, paste, collapse = ", ") %>%
  # Convert to a data frame with a proper column name
  data.frame(group_values = .) %>%
  # Now n() will work inside mutate()
  mutate(grp = 1:nrow(.)) %>%
  relocate(grp) %>%
  knitr::kable()

core.methods.se <- c("LRa-FE-db", "LRa-FE-het", "LRa-FIcw-db", "LRapw-FE-db", 
                     "LRapw-FE-het", "LRapw-FIpw-db", "LRapw-FIpw-het", 
                     "LRicw-FE-crve", "LRicw-FIcw-crve", "MLM-FE", "MLM-FIcw",
                     "MLM-RE")


# QUESTION: Can we delete these?
rm( groups )
rm( sruns )


