# Author: Worasit Sangjan
# Date: July 7, 2025

# ============================================================================
# BACK-TRANSFORMATION FOR CUBE ROOT TRANSFORMED LMM RESULTS
# Converts transformed results back to original temperature scale (°C)
# Works with both delta_Ta and delta_Tn analyses
# ============================================================================

# Load required libraries
library(emmeans)

# ============================================================================
# BACK-TRANSFORMATION FUNCTIONS
# ============================================================================

# Simple back-transformation for point estimates
back_transform_cube_root <- function(transformed_value) {
  return(sign(transformed_value) * abs(transformed_value)^3)
}

# Proper confidence intervals using delta method
get_corrected_ci_delta_method <- function(emmean_trans, se_trans, alpha = 0.05) {
  # Delta method for cube root: derivative = 3 * Y^2
  derivative <- 3 * emmean_trans^2
  se_original <- abs(se_trans * derivative)
  emmean_original <- back_transform_cube_root(emmean_trans)
  
  z_crit <- qnorm(1 - alpha/2)
  lower_ci <- emmean_original - z_crit * se_original
  upper_ci <- emmean_original + z_crit * se_original
  
  return(list(
    emmean_original = emmean_original,
    se_original = se_original,
    lower_ci = lower_ci,
    upper_ci = upper_ci
  ))
}

# Main back-transformation function
back_transform_emmeans <- function(emm_object, method_name = "unknown") {
  cat("Back-transforming", method_name, "...\n")
  
  # Try regrid() method first (preferred)
  tryCatch({
    emm_original <- regrid(emm_object, transform = "identity")
    cat("Yes: Used regrid() method for", method_name, "\n")
    return(list(method = "regrid", result = emm_original, success = TRUE))
  }, error = function(e) {
    cat("Using delta method for", method_name, "\n")
    
    # Apply delta method manually
    emm_summary <- summary(emm_object)
    result_df <- emm_summary
    
    for(i in 1:nrow(emm_summary)) {
      corrected <- get_corrected_ci_delta_method(emm_summary$emmean[i], emm_summary$SE[i])
      result_df$emmean_original[i] <- corrected$emmean_original
      result_df$SE_original[i] <- corrected$se_original
      result_df$lower_CI_original[i] <- corrected$lower_ci
      result_df$upper_CI_original[i] <- corrected$upper_ci
    }
    
    cat("Yes: Used delta method for", method_name, "\n")
    return(list(method = "delta_method", result = result_df, success = TRUE))
  })
}

# ============================================================================
# BACK-TRANSFORM ALL MODELS
# ============================================================================

cat("=== BACK-TRANSFORMING LMM RESULTS TO ORIGINAL SCALE ===\n")

# Initialize results list
back_transformed_results <- list()

# Model 1: Treatment × Thermal Phase
if(exists("emm_1")) {
  result_1 <- back_transform_emmeans(emm_1, "Model 1")
  if(result_1$success) {
    if(result_1$method == "regrid") {
      back_transformed_results$model_1 <- summary(result_1$result)
      # Standardize column names
      if("response" %in% colnames(back_transformed_results$model_1)) {
        back_transformed_results$model_1$emmean_original <- back_transformed_results$model_1$response
      }
      if("lower.CL" %in% colnames(back_transformed_results$model_1)) {
        back_transformed_results$model_1$lower_CI_original <- back_transformed_results$model_1$lower.CL
        back_transformed_results$model_1$upper_CI_original <- back_transformed_results$model_1$upper.CL
      }
    } else {
      back_transformed_results$model_1 <- result_1$result
    }
  }
}

# Model 2: Treatment × Time Block
if(exists("emm_2")) {
  result_2 <- back_transform_emmeans(emm_2, "Model 2")
  if(result_2$success) {
    if(result_2$method == "regrid") {
      back_transformed_results$model_2 <- summary(result_2$result)
      if("response" %in% colnames(back_transformed_results$model_2)) {
        back_transformed_results$model_2$emmean_original <- back_transformed_results$model_2$response
      }
      if("lower.CL" %in% colnames(back_transformed_results$model_2)) {
        back_transformed_results$model_2$lower_CI_original <- back_transformed_results$model_2$lower.CL
        back_transformed_results$model_2$upper_CI_original <- back_transformed_results$model_2$upper.CL
      }
    } else {
      back_transformed_results$model_2 <- result_2$result
    }
  }
}

# Model 3: Full Factorial
if(exists("emm_3")) {
  result_3 <- back_transform_emmeans(emm_3, "Model 3")
  if(result_3$success) {
    if(result_3$method == "regrid") {
      back_transformed_results$model_3 <- summary(result_3$result)
      if("response" %in% colnames(back_transformed_results$model_3)) {
        back_transformed_results$model_3$emmean_original <- back_transformed_results$model_3$response
      }
      if("lower.CL" %in% colnames(back_transformed_results$model_3)) {
        back_transformed_results$model_3$lower_CI_original <- back_transformed_results$model_3$lower.CL
        back_transformed_results$model_3$upper_CI_original <- back_transformed_results$model_3$upper.CL
      }
    } else {
      back_transformed_results$model_3 <- result_3$result
    }
  }
}

# ============================================================================
# DISPLAY RESULTS
# ============================================================================

# Function to display results nicely
display_back_transformed_results <- function(results, model_name) {
  if(is.null(results)) {
    cat("No results available for", model_name, "\n")
    return()
  }
  
  cat("\n=== ", model_name, " - BACK-TRANSFORMED RESULTS (°C) ===\n")
  
  # Check for missing values
  valid_rows <- !is.na(results$emmean_original) & 
                !is.na(results$lower_CI_original) & 
                !is.na(results$upper_CI_original)
  
  if(sum(valid_rows) == 0) {
    cat("No valid results to display\n")
    return()
  }
  
  results <- results[valid_rows, ]
  
  # For Model 3 with both thermal_phase and time_block
  if("thermal_phase" %in% colnames(results) && "time_block" %in% colnames(results)) {
    # Group by thermal_phase first, then time_block
    for(phase in unique(results$thermal_phase)) {
      cat("\n--- ", phase, " ---\n")
      phase_data <- results[results$thermal_phase == phase, ]
      
      for(time_block in unique(phase_data$time_block)) {
        if(!is.na(time_block) && time_block != "All_Combined") {
          cat("   ", time_block, ":\n")
        }
        
        block_data <- phase_data[phase_data$time_block == time_block, ]
        
        for(i in 1:nrow(block_data)) {
          treatment <- block_data$treatment[i]
          temp_effect <- block_data$emmean_original[i]
          lower_ci <- block_data$lower_CI_original[i]
          upper_ci <- block_data$upper_CI_original[i]
          
          # Handle missing values safely
          if(!is.na(temp_effect) && !is.na(lower_ci) && !is.na(upper_ci)) {
            effect_desc <- if(temp_effect > 0) {
              sprintf("%.4f°C WARMER", temp_effect)
            } else if(temp_effect < 0) {
              sprintf("%.4f°C COOLER", abs(temp_effect))
            } else {
              "No difference"
            }
            
            cat(sprintf("    %s: %s [95%% CI: %.4f to %.4f°C]\n", 
                        treatment, effect_desc, lower_ci, upper_ci))
          }
        }
      }
    }
  }
  # Group by thermal_phase only (Model 1)
  else if("thermal_phase" %in% colnames(results)) {
    for(phase in unique(results$thermal_phase)) {
      cat("\n--- ", phase, " ---\n")
      phase_data <- results[results$thermal_phase == phase, ]
      
      for(i in 1:nrow(phase_data)) {
        treatment <- phase_data$treatment[i]
        temp_effect <- phase_data$emmean_original[i]
        lower_ci <- phase_data$lower_CI_original[i]
        upper_ci <- phase_data$upper_CI_original[i]
        
        if(!is.na(temp_effect) && !is.na(lower_ci) && !is.na(upper_ci)) {
          effect_desc <- if(temp_effect > 0) {
            sprintf("%.4f°C WARMER", temp_effect)
          } else if(temp_effect < 0) {
            sprintf("%.4f°C COOLER", abs(temp_effect))
          } else {
            "No difference"
          }
          
          cat(sprintf("  %s: %s [95%% CI: %.4f to %.4f°C]\n", 
                      treatment, effect_desc, lower_ci, upper_ci))
        }
      }
    }
  }
  # Group by time_block only (Model 2)
  else if("time_block" %in% colnames(results)) {
    for(time_block in unique(results$time_block)) {
      cat("\n--- ", time_block, " ---\n")
      block_data <- results[results$time_block == time_block, ]
      
      for(i in 1:nrow(block_data)) {
        treatment <- block_data$treatment[i]
        temp_effect <- block_data$emmean_original[i]
        lower_ci <- block_data$lower_CI_original[i]
        upper_ci <- block_data$upper_CI_original[i]
        
        if(!is.na(temp_effect) && !is.na(lower_ci) && !is.na(upper_ci)) {
          effect_desc <- if(temp_effect > 0) {
            sprintf("%.4f°C WARMER", temp_effect)
          } else if(temp_effect < 0) {
            sprintf("%.4f°C COOLER", abs(temp_effect))
          } else {
            "No difference"
          }
          
          cat(sprintf("  %s: %s [95%% CI: %.4f to %.4f°C]\n", 
                      treatment, effect_desc, lower_ci, upper_ci))
        }
      }
    }
  }
}

# Display all results
display_back_transformed_results(back_transformed_results$model_1, "MODEL 1")
display_back_transformed_results(back_transformed_results$model_2, "MODEL 2") 
display_back_transformed_results(back_transformed_results$model_3, "MODEL 3")

# ============================================================================
# CREATE OUTPUT FILES
# ============================================================================

cat("\n=== CREATING OUTPUT FILES ===\n")

# Function to create standardized CSV
create_back_transformed_csv <- function(results, model_name, variable_name) {
  if(is.null(results)) return(NULL)
  
  output_table <- data.frame(
    Model = model_name,
    Treatment = results$treatment,
    Thermal_Phase = if("thermal_phase" %in% colnames(results)) results$thermal_phase else "All_Combined",
    Time_Block = if("time_block" %in% colnames(results)) results$time_block else "All_Combined",
    EMM_degrees_C = round(results$emmean_original, 4),
    SE_original_scale = if("SE_original" %in% colnames(results)) round(results$SE_original, 4) else NA,
    Lower_CI_degrees_C = round(results$lower_CI_original, 4),
    Upper_CI_degrees_C = round(results$upper_CI_original, 4),
    Variable_Analyzed = variable_name,
    stringsAsFactors = FALSE
  )
  
  # Add effect description
  output_table$Temperature_Effect <- sapply(output_table$EMM_degrees_C, function(x) {
    if(is.na(x)) return("Missing_Data")
    if(x > 0) return(sprintf("%.4f°C_Warmer", x))
    if(x < 0) return(sprintf("%.4f°C_Cooler", abs(x)))
    return("No_Effect")
  })
  
  return(output_table)
}

# Detect which variable was analyzed
if(exists("selected_description")) {
  variable_analyzed <- selected_description
} else if(exists("best_var_name")) {
  variable_analyzed <- best_var_name
} else {
  variable_analyzed <- "delta_temperature"
}

# Create CSV files
if(!is.null(back_transformed_results$model_1)) {
  table_1 <- create_back_transformed_csv(back_transformed_results$model_1, "Model_1_Thermal_Phase", variable_analyzed)
  write.csv(table_1, "Model_1_BackTransformed_Results.csv", row.names = FALSE)
  cat("Saved: Model_1_BackTransformed_Results.csv\n")
}

if(!is.null(back_transformed_results$model_2)) {
  table_2 <- create_back_transformed_csv(back_transformed_results$model_2, "Model_2_Time_Block", variable_analyzed)
  write.csv(table_2, "Model_2_BackTransformed_Results.csv", row.names = FALSE)
  cat("Saved: Model_2_BackTransformed_Results.csv\n")
}

if(!is.null(back_transformed_results$model_3)) {
  table_3 <- create_back_transformed_csv(back_transformed_results$model_3, "Model_3_Full_Factorial", variable_analyzed)
  write.csv(table_3, "Model_3_BackTransformed_Results.csv", row.names = FALSE)
  cat("Saved: Model_3_BackTransformed_Results.csv\n")
}

# Create treatment comparison
if(length(back_transformed_results) > 0) {
  treatments <- c("Paper_low", "Paper_raise", "Plastic_low", "Plastic_raise")
  
  treatment_comparison <- data.frame(
    Treatment = treatments,
    Model_1_Mean_Effect = NA,
    Model_2_Mean_Effect = NA,
    Model_3_Mean_Effect = NA,
    stringsAsFactors = FALSE
  )
  
  for(treatment in treatments) {
    if(!is.null(back_transformed_results$model_1)) {
      effects_1 <- back_transformed_results$model_1$emmean_original[back_transformed_results$model_1$treatment == treatment]
      if(length(effects_1) > 0) {
        treatment_comparison$Model_1_Mean_Effect[treatment_comparison$Treatment == treatment] <- 
          round(mean(effects_1, na.rm = TRUE), 4)
      }
    }
    
    if(!is.null(back_transformed_results$model_2)) {
      effects_2 <- back_transformed_results$model_2$emmean_original[back_transformed_results$model_2$treatment == treatment]
      if(length(effects_2) > 0) {
        treatment_comparison$Model_2_Mean_Effect[treatment_comparison$Treatment == treatment] <- 
          round(mean(effects_2, na.rm = TRUE), 4)
      }
    }
    
    if(!is.null(back_transformed_results$model_3)) {
      effects_3 <- back_transformed_results$model_3$emmean_original[back_transformed_results$model_3$treatment == treatment]
      if(length(effects_3) > 0) {
        treatment_comparison$Model_3_Mean_Effect[treatment_comparison$Treatment == treatment] <- 
          round(mean(effects_3, na.rm = TRUE), 4)
      }
    }
  }
  
  # Add overall ranking
  best_col <- if(!is.null(back_transformed_results$model_3)) "Model_3_Mean_Effect" else 
              if(!is.null(back_transformed_results$model_2)) "Model_2_Mean_Effect" else "Model_1_Mean_Effect"
  
  treatment_comparison$Overall_Rank <- rank(-treatment_comparison[[best_col]], ties.method = "min", na.last = "keep")
  
  write.csv(treatment_comparison, "Treatment_Comparison_BackTransformed.csv", row.names = FALSE)
  cat("Saved: Treatment_Comparison_BackTransformed.csv\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== BACK-TRANSFORMATION COMPLETE ===\n")
cat("Variable analyzed:", variable_analyzed, "\n")
cat("Models back-transformed:", length(back_transformed_results), "\n")
cat("Method: Statistically correct confidence intervals\n")

cat("\nFiles created:\n")
if(!is.null(back_transformed_results$model_1)) cat("- Model_1_BackTransformed_Results.csv\n")
if(!is.null(back_transformed_results$model_2)) cat("- Model_2_BackTransformed_Results.csv\n")
if(!is.null(back_transformed_results$model_3)) cat("- Model_3_BackTransformed_Results.csv\n")
if(length(back_transformed_results) > 0) cat("- Treatment_Comparison_BackTransformed.csv\n")

cat("\nThese results are ready for reporting in °C\n")
