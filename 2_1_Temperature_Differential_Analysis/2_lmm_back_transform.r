# Author: Worasit Sangjan
# Date: September 6, 2025
# Version 2

# ============================================================================
# BACK-TRANSFORMATION FOR LOG-TRANSFORMED LMM RESULTS
# Converts log-transformed EMM results back to original temperature scale (°C)
# 
# Purpose: Back-transform log-scale LMM results for interpretation and reporting
# Note: Statistical analysis remains on log scale - this is interpretation only
# ============================================================================

# Load required libraries
library(dplyr)

# ============================================================================
# SHIFT VALUE CALCULATION
# ============================================================================

cat("=== BACK-TRANSFORMATION: LOG TO ORIGINAL SCALE ===\n")

# Calculate shift values from original data
df <- read.csv("grow_tube_categorized.csv")
DELTA_TA_SHIFT <- abs(min(df$delta_Ta, na.rm = TRUE)) + 0.01
DELTA_TN_SHIFT <- abs(min(df$delta_Tn, na.rm = TRUE)) + 0.01

cat("Shift values: Delta_Ta =", round(DELTA_TA_SHIFT, 4), "| Delta_Tn =", round(DELTA_TN_SHIFT, 4), "\n\n")

# ============================================================================
# BACK-TRANSFORMATION FUNCTIONS
# ============================================================================

# Log back-transformation with delta method for confidence intervals
back_transform_log_ci <- function(log_emmean, log_se, shift_value, alpha = 0.05) {
  # Back-transform point estimate
  original_emmean <- exp(log_emmean) - shift_value
  
  # Delta method for standard error
  derivative <- exp(log_emmean)
  original_se <- log_se * derivative
  
  # Calculate confidence intervals
  z_crit <- qnorm(1 - alpha/2)
  lower_ci <- original_emmean - z_crit * original_se
  upper_ci <- original_emmean + z_crit * original_se
  
  return(list(
    emmean_original = original_emmean,
    se_original = original_se,
    lower_ci = lower_ci,
    upper_ci = upper_ci
  ))
}

# Process CSV results
back_transform_csv <- function(csv_file, shift_value, variable_name) {
  
  if (!file.exists(csv_file)) {
    cat("Warning: File", csv_file, "not found\n")
    return(NULL)
  }
  
  results <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Check if log transformation was used
  if (!any(grepl("log", results$Variable_Used, ignore.case = TRUE))) {
    cat("Warning:", variable_name, "was not log-transformed\n")
    return(NULL)
  }
  
  cat("Processing:", variable_name, "with", nrow(results), "rows\n")
  
  # Back-transform all rows
  back_transformed <- results
  
  for (i in 1:nrow(results)) {
    bt_result <- back_transform_log_ci(
      log_emmean = results$EMM[i],
      log_se = results$SE[i], 
      shift_value = shift_value
    )
    
    # Add back-transformed values
    back_transformed$EMM_Original_C[i] <- round(bt_result$emmean_original, 4)
    back_transformed$SE_Original[i] <- round(bt_result$se_original, 4)
    back_transformed$Lower_CI_Original_C[i] <- round(bt_result$lower_ci, 4)
    back_transformed$Upper_CI_Original_C[i] <- round(bt_result$upper_ci, 4)
    
    # Add temperature effect description
    temp_effect <- bt_result$emmean_original
    if (is.na(temp_effect)) {
      back_transformed$Temperature_Effect[i] <- "Missing_data"
    } else if (temp_effect > 0) {
      back_transformed$Temperature_Effect[i] <- paste0("+", round(temp_effect, 4), "°C")
    } else if (temp_effect < 0) {
      back_transformed$Temperature_Effect[i] <- paste0(round(temp_effect, 4), "°C")
    } else {
      back_transformed$Temperature_Effect[i] <- "0°C"
    }
  }
  
  return(back_transformed)
}

# ============================================================================
# PROCESS DATASETS
# ============================================================================

# Back-transform both datasets
delta_ta_back <- back_transform_csv("Combined_LMM_Delta_Ta_Complete_Results.csv", DELTA_TA_SHIFT, "delta_Ta")
delta_tn_back <- back_transform_csv("Combined_LMM_Delta_Tn_Complete_Results.csv", DELTA_TN_SHIFT, "delta_Tn")

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n=== SAVING RESULTS ===\n")

if (!is.null(delta_ta_back)) {
  write.csv(delta_ta_back, "Back_Transformed_Delta_Ta_Results.csv", row.names = FALSE)
  cat("Saved: Back_Transformed_Delta_Ta_Results.csv\n")
}

if (!is.null(delta_tn_back)) {
  write.csv(delta_tn_back, "Back_Transformed_Delta_Tn_Results.csv", row.names = FALSE)
  cat("Saved: Back_Transformed_Delta_Tn_Results.csv\n")
}

# Create comparison table if both datasets processed
if (!is.null(delta_ta_back) && !is.null(delta_tn_back)) {
  
  # Extract best model results
  ta_best <- delta_ta_back[delta_ta_back$Best_Model == "YES", ]
  tn_best <- delta_tn_back[delta_tn_back$Best_Model == "YES", ]
  
  # Create treatment comparison
  treatments <- intersect(unique(ta_best$Treatment), unique(tn_best$Treatment))
  
  comparison_table <- data.frame(
    Treatment = treatments,
    Delta_Ta_Effect_C = NA,
    Delta_Ta_CI_Lower = NA,
    Delta_Ta_CI_Upper = NA,
    Delta_Tn_Effect_C = NA,
    Delta_Tn_CI_Lower = NA,
    Delta_Tn_CI_Upper = NA,
    stringsAsFactors = FALSE
  )
  
  for (treatment in treatments) {
    # Delta Ta values
    ta_rows <- ta_best[ta_best$Treatment == treatment, ]
    if (nrow(ta_rows) > 0) {
      comparison_table$Delta_Ta_Effect_C[comparison_table$Treatment == treatment] <- 
        round(mean(ta_rows$EMM_Original_C, na.rm = TRUE), 3)
      comparison_table$Delta_Ta_CI_Lower[comparison_table$Treatment == treatment] <- 
        round(mean(ta_rows$Lower_CI_Original_C, na.rm = TRUE), 3)
      comparison_table$Delta_Ta_CI_Upper[comparison_table$Treatment == treatment] <- 
        round(mean(ta_rows$Upper_CI_Original_C, na.rm = TRUE), 3)
    }
    
    # Delta Tn values
    tn_rows <- tn_best[tn_best$Treatment == treatment, ]
    if (nrow(tn_rows) > 0) {
      comparison_table$Delta_Tn_Effect_C[comparison_table$Treatment == treatment] <- 
        round(mean(tn_rows$EMM_Original_C, na.rm = TRUE), 3)
      comparison_table$Delta_Tn_CI_Lower[comparison_table$Treatment == treatment] <- 
        round(mean(tn_rows$Lower_CI_Original_C, na.rm = TRUE), 3)
      comparison_table$Delta_Tn_CI_Upper[comparison_table$Treatment == treatment] <- 
        round(mean(tn_rows$Upper_CI_Original_C, na.rm = TRUE), 3)
    }
  }
  
  write.csv(comparison_table, "Temperature_Effects_Comparison.csv", row.names = FALSE)
  cat("Saved: Temperature_Effects_Comparison.csv\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== BACK-TRANSFORMATION COMPLETE ===\n")
cat("Results converted from log scale to original temperature scale (°C)\n")
cat("Statistical analysis remains on log scale - use these for interpretation only\n")