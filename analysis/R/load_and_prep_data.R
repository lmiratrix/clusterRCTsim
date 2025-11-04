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
require(blkvar)


load( here::here( "results_agg/combined_simulation_results.RData" ) )

# replace 'block' with 'block' for all names in sruns dataframe:
sruns <- sruns %>%
  rename_with(~ str_replace_all(., "district", "block")) %>%
  rename_with(~ str_replace_all(., ".dist", ".block")) %>%
  rename_with(~ str_replace_all(., "school", "cluster")) %>%
  rename_with(~ str_replace_all(., ".sch", ".clu"))
  
  
  


#figure out what didn't run in the low df 
lowdf<-sruns %>% filter(group=="lowdf")
with(lowdf,table(is.na(mean.SE.hat),method))
with(lowdf,table(is.na(mean.ATE.hat),method))
#there are estimates for all the methods for all the scenarios, but not SEs 

     
DGP_vars <- c("K", "J.bar", "p.bar", "sd.ATE", "ICC.block", "variable.p", 
              "n.bar", "sd.cluster", "ICC.cluster", "size.impact.cor", "size.impact.cor.cluster", 
              "proptx.impact.cor", "num.cov.block", "num.cov.cluster", "cv.cluster", 
              "num.cov.indiv", "empirical", "cv.block", "nbar.J.cor", 
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

# Function to create a single labeled factor variable based on group-specific labeling
create_labeled_factor_by_group <- function(data, variable, group_var) {
  # Convert to character to handle column names
  var_name <- deparse(substitute(variable))
  group_name <- deparse(substitute(group_var))
  
  # If passed as strings, use them directly
  if (is.character(substitute(variable))) var_name <- variable
  if (is.character(substitute(group_var))) group_name <- group_var
  
  # Create a helper function to apply labeling within each group
  label_within_group <- function(x) {
    unique_values <- sort(unique(x))
    num_levels <- length(unique_values)
    
    # Check if 0 is among the values in this group
    has_zero <- 0 %in% unique_values
    
    if (num_levels == 1) {
      # If only one level, just return as is
      factor(x)
    } else if (num_levels == 2) {
      if (has_zero) {
        # If one of the two values is 0, use "0" and "+"
        labels <- ifelse(unique_values == 0, "0", "+")
        factor(x, levels = unique_values, labels = labels)
      } else {
        # If no zero, use "S" and "L" (Small and Large)
        labels <- c("S", "L")
        factor(x, levels = unique_values, labels = labels)
      }
    } else if (num_levels == 3) {
      if (has_zero) {
        # If one of the three values is 0, use "0", "S", and "L"
        zero_pos <- which(unique_values == 0)
        ordered_levels <- c(0, unique_values[unique_values != 0])
        labels <- c("0", "S", "L")
        factor(x, levels = ordered_levels, labels = labels)
      } else {
        # If no zero, use "S", "M", and "L" (Small, Medium, Large)
        labels <- c("S", "M", "L")
        factor(x, levels = unique_values, labels = labels)
      }
    } else {
      # For more than three levels
      if (has_zero) {
        # Reorder so 0 is first
        ordered_levels <- c(0, unique_values[unique_values != 0])
        # First label is "0", rest are numbered L1, L2, etc.
        labels <- c("0", paste0("L", 1:(num_levels-1)))
        factor(x, levels = ordered_levels, labels = labels)
      } else {
        # If no zero value exists, just use L1, L2, etc.
        labels <- paste0("L", 1:num_levels)
        factor(x, levels = unique_values, labels = labels)
      }
    }
  }
  
  # Apply the labeling within each group and then combine
  data %>%
    group_by(!!sym(group_name)) %>%
    mutate(!!var_name := label_within_group(!!sym(var_name))) %>%
    ungroup()
}

# Example usage:
# sruns.clean <- create_labeled_factor_by_group(sruns.clean, "J.bar", "group")
# 
# This will create a single J.bar variable where:
# - Rows with group="full_factor": 10→S, 20→M, 30→L
# - Rows with group="lowdf": 3→S, 4→L
# - All in the same J.bar column


# Drop impossible combos ---

# Check combos of simulations
table( cv.cluster = sruns$cv.cluster, 
       sd.cluster = sruns$sd.cluster,
       size.impact.cor.cluster = sruns$size.impact.cor.cluster )


table( sruns$variable.p, sruns$proptx.impact.cor )

table( sruns$size.impact.cor, sruns$sd.cluster )


sruns <- sruns %>%
  filter( !(size.impact.cor.cluster != 0 & sd.cluster == 0) ) %>%
  filter( !(size.impact.cor.cluster != 0 & cv.cluster == 0) ) 







# Make clean set ----

#I'm going to keep the two different simulation sets with separate block type variables for now 
#I need to fix the levels by the groups 

sruns.clean <- sruns %>%  
  # Create new ID combining old ID and group
  mutate(ID = paste(ID, group, sep = "_")) %>%
  # Create a different block type variable for the full factor and the low df data sets 
  mutate(block_type = case_when(
    # Blocks all the same (btype1)
    cv.block == 0 &
      sd.ATE == 0 & 
      size.impact.cor == 0 & 
      nbar.J.cor == 0 &
      proptx.impact.cor == 0 & 
      size.intercept.cor == 0 &
      !variable.p &
      cv.nbar == 0 &
      group=="full_factor" ~ 1,
    
    # Blocks all messy (btype2)
    cv.block == 0.85 &
      sd.ATE == 0.6 & 
      size.impact.cor == 0.8 & 
      nbar.J.cor == 0.8 &
      proptx.impact.cor == 0.8 & 
      size.intercept.cor == 0.8 &
      variable.p &
      cv.nbar == 0.33 & group=="full_factor" ~ 3,
    
    # Blocks vaguely reasonable (btype3)
    cv.block == 0.33 &
      sd.ATE == 0.2 & 
      size.impact.cor == 0.2 & 
      nbar.J.cor == 0.3 &
      proptx.impact.cor == 0 & 
      size.intercept.cor == 0.3 &
      !variable.p &
      cv.nbar == 0.33 & group=="full_factor"~ 2,
    
    TRUE ~ NA_real_  # Catch any other combinations
  )) %>% 
  mutate(block_type_lowdf = case_when(
    # Blocks all the same (btype1)
    sd.ATE == 0 & 
      size.impact.cor == 0 & 
      size.intercept.cor == 0 & group=="lowdf" ~ 1,
    
    # Blocks all messy (btype2)
    sd.ATE == 0.6 & 
      size.impact.cor == 0.8 & 
      size.intercept.cor == 0.8 & group=="lowdf"~ 3,
    
    # Blocks vaguely reasonable (btype3)
    sd.ATE == 0.2 & 
      size.impact.cor == 0.2 & 
      size.intercept.cor == 0.3 & group=="lowdf"~ 2,
    
    TRUE ~ NA_real_  # Catch any other combinations
  )) %>% 
  mutate(
    bias.s.estimand = case_when(
      weight == "Cluster" ~ bias.s.clu,
      weight == "Person" ~ bias.s.stu
    ),
    cover.s.estimand = case_when(
      weight == "Cluster" ~ cover.s.clu,
      weight == "Person" ~ cover.s.stu
    ),
    RMSE.s.estimand = case_when(
      weight == "Cluster" ~ RMSE.s.clu,
      weight == "Person" ~ RMSE.s.stu
    )
  ) %>%
  # Apply group-based labeling to variables that differ by group
  create_labeled_factor_by_group("sd.ATE", "group") %>%
  create_labeled_factor_by_group("size.impact.cor", "group") %>%
  create_labeled_factor_by_group("J.bar", "group") %>%
  create_labeled_factor_by_group("sd.cluster", "group") %>%
  create_labeled_factor_by_group("ICC.cluster", "group") %>%
  create_labeled_factor_by_group("cv.cluster", "group") %>%
  create_labeled_factor_by_group("num.cov.cluster", "group") %>%
  create_labeled_factor_by_group("proptx.impact.cor", "group") %>%
  create_labeled_factor_by_group("size.intercept.cor", "group") %>%
  create_labeled_factor_by_group("size.impact.cor.cluster", "group") %>%
  mutate(
    # Fix for the factor problem - explicitly define levels based on your block types
    block_type = factor(block_type, 
                        levels = c(1, 2, 3),
                        labels = c("Same", "Reas", "Mess")),
    block_type_lowdf = factor(block_type_lowdf, 
                              levels = c(1, 2, 3),
                              labels = c("Same", "Reas", "Mess")),
    p.bar = factor(p.bar,
                   levels=c(.5,.7),
                   labels = c("Bal","UnBal")),
    bias.s.estimand.mag = abs(bias.s.estimand),
    method = reorder(method, -biased),
    bias.ratio.clu = round(bias.f.clu/bias.s.clu, 4),
    bias.ratio.stu = round(bias.f.stu/bias.s.stu, 4),
    RMSE.ratio.clu = round(RMSE.f.clu/RMSE.s.clu, 4),
    RMSE.ratio.stu = round(RMSE.f.stu/RMSE.s.stu, 4),
    cover.ratio.stu = round(cover.f.stu/cover.s.stu, 4),
    cover.ratio.clu = round(cover.f.clu/cover.s.clu, 4),    
    true.se.ratio = round(SE.f/SE.s, 2)
  ) %>%
  mutate(
    is_FI = str_detect(method, "FI"),
    is_MLM = str_detect(method, "MLM"),
    is_cw = str_detect(method, "cw"),
    is_pw = str_detect(method, "pw"),
  )

# Check the levels of one of the variables
levels(sruns.clean$block_type)
table(sruns.clean$block_type)
levels(sruns.clean$block_type_lowdf)
table(sruns.clean$block_type_lowdf)
levels(sruns.clean$J.bar)
levels(sruns.clean$ICC.cluster)
levels(sruns.clean$sd.ATE)
levels(sruns.clean$sd.cluster)
table(sruns.clean$ID)


# Create the new variable with median RMSE.s.estimand for each ID
sruns.clean <- sruns.clean %>%
  group_by(ID) %>%
  mutate(RMSE.s.estimand.median = median(RMSE.s.estimand, na.rm = TRUE),
         SE.s.median = median(SE.s, na.rm = TRUE),
         CI.length.median = median(CI.length, na.rm = TRUE),
         sd.SE.hat.s.median=median(sd.SE.hat.s,na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(ID,weight) %>%
  mutate(SE.s.median.weight = median(SE.s, na.rm = TRUE)) %>%
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

sruns.clean <- sruns.clean %>%
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


groups <- group_estimators(sruns.clean , var = "mean.ATE.hat")




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

####Pick Up Here

#the missing SE's for some methods in the low df group  is doing weird stuff, so run them separately 

groups.ff <- group_estimators(sruns.clean %>% filter(group=="full_factor"), var = "mean.SE.hat")
groups.lowdf <- group_estimators(sruns.clean %>% filter(group=="lowdf"), var = "mean.SE.hat")


# List the groups
map_chr(groups.ff, paste, collapse = ", ") %>%
  # Convert to a data frame with a proper column name
  data.frame(group_values = .) %>%
  # Now n() will work inside mutate()
  mutate(grp = 1:nrow(.)) %>%
  relocate(grp) %>%
  knitr::kable()

map_chr(groups.lowdf, paste, collapse = ", ") %>%
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


# QUESTION: Can we delete these?, SL: yes! 
rm( groups )
rm(groups.ff)
rm(groups.lowdf)
rm( sruns )


