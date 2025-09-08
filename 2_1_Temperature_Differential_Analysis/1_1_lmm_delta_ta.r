# Author: Worasit Sangjan
# Date: September 5, 2025
# Version 2

# ============================================================================
# LINEAR MIXED-EFFECTS MODEL ANALYSIS FOR GROW TUBE TEMPERATURE EFFECTS
# Analysis of delta_Ta (tube temperature - air temperature)
# 
# This script performs comprehensive LMM analysis with:
# 1. Unified transformation assessment for optimal model selection
# 2. Three modeling approaches with comprehensive diagnostics
# 3. Combined results output 
# ============================================================================

# Load required libraries
library(lme4)
library(emmeans)
library(car)
library(multcomp)
library(dplyr)

# ============================================================================
# DATA PREPARATION
# ============================================================================

cat("=== LOADING AND PREPARING DATA ===\n")

# Load data and remove control group
df <- read.csv("grow_tube_categorized.csv")
df <- df[df$Material != "Uncover", ]

# Create factor variables
df$treatment <- as.factor(paste(df$Material, df$Condition, sep="_"))
df$tube_id <- as.factor(paste(df$treatment, "Rep", df$Replicate, sep="_"))
df$date_factor <- as.factor(df$date)
df$thermal_phase <- as.factor(df$thermal_phase)
df$time_block <- as.factor(df$time_block)

cat("Dataset:", nrow(df), "observations\n")
cat("Treatments:", paste(levels(df$treatment), collapse = ", "), "\n")

# ============================================================================
# TRANSFORMATION ASSESSMENT (UNIFIED DATASET APPROACH)
# ============================================================================

cat("\n=== TRANSFORMATION ASSESSMENT ===\n")

# Model control parameters
control_opts <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5),
  calc.derivs = FALSE,
  check.nobs.vs.nlev = "ignore",
  check.nobs.vs.rankZ = "ignore"
)

# Calculate shift value for log/sqrt transformations
min_val <- min(df$delta_Ta, na.rm = TRUE)
shift_val <- abs(min_val) + 0.01

# Create all transformation variables
df$delta_Ta_log <- log(df$delta_Ta + shift_val)
df$delta_Ta_sqrt <- sqrt(df$delta_Ta + shift_val)
df$delta_Ta_cbrt <- sign(df$delta_Ta) * abs(df$delta_Ta)^(1/3)

# Test all transformations using unified model structure
model_original <- lmer(delta_Ta ~ treatment * thermal_phase + (1|tube_id) + (1|date_factor), 
                       data = df, REML = FALSE, control = control_opts)
model_log <- lmer(delta_Ta_log ~ treatment * thermal_phase + (1|tube_id) + (1|date_factor), 
                  data = df, REML = FALSE, control = control_opts)
model_sqrt <- lmer(delta_Ta_sqrt ~ treatment * thermal_phase + (1|tube_id) + (1|date_factor), 
                   data = df, REML = FALSE, control = control_opts)
model_cbrt <- lmer(delta_Ta_cbrt ~ treatment * thermal_phase + (1|tube_id) + (1|date_factor), 
                   data = df, REML = FALSE, control = control_opts)

# Test normality for all transformations
test_normality <- function(model) {
  residuals_sample <- sample(resid(model), min(5000, length(resid(model))))
  shapiro_test <- shapiro.test(residuals_sample)
  return(list(p_value = shapiro_test$p.value, statistic = shapiro_test$statistic))
}

norm_orig <- test_normality(model_original)
norm_log <- test_normality(model_log)
norm_sqrt <- test_normality(model_sqrt)
norm_cbrt <- test_normality(model_cbrt)

# Create transformation comparison table
transformation_summary <- data.frame(
  Transformation = c("Original", "Log", "SQRT", "CBRT"),
  AIC = c(AIC(model_original), AIC(model_log), AIC(model_sqrt), AIC(model_cbrt)),
  Shapiro_P = c(norm_orig$p_value, norm_log$p_value, norm_sqrt$p_value, norm_cbrt$p_value),
  Shapiro_W = c(norm_orig$statistic, norm_log$statistic, norm_sqrt$statistic, norm_cbrt$statistic)
)

cat("Transformation comparison:\n")
print(transformation_summary)

# Select transformation with lowest AIC
best_aic_idx <- which.min(transformation_summary$AIC)
best_transformation <- transformation_summary$Transformation[best_aic_idx]

if(best_transformation == "Original") {
  selected_variable <- "delta_Ta"
  selected_description <- "original delta_Ta"
  selected_model <- model_original
} else if(best_transformation == "Log") {
  selected_variable <- "delta_Ta_log"
  selected_description <- "log transformed delta_Ta"
  selected_model <- model_log
} else if(best_transformation == "SQRT") {
  selected_variable <- "delta_Ta_sqrt"
  selected_description <- "square root transformed delta_Ta"
  selected_model <- model_sqrt
} else if(best_transformation == "CBRT") {
  selected_variable <- "delta_Ta_cbrt"
  selected_description <- "cube root transformed delta_Ta"
  selected_model <- model_cbrt
}

cat("*** SELECTED:", best_transformation, "transformation (lowest AIC =", round(min(transformation_summary$AIC), 2), ") ***\n")

# ============================================================================
# DIAGNOSTIC PLOTS
# ============================================================================

# Original data diagnostics
png("Diagnostic_Plots_Original_vs_Selected.png", width = 3600, height = 2400, res = 250)
par(mfrow = c(2, 3), 
    mar = c(5, 5, 4, 2),
    cex.main = 1.8,
    cex.lab = 1.5,
    cex.axis = 1.3,
    font.main = 2,
    font.lab = 2)

# Original data plots
plot(fitted(model_original), resid(model_original), 
     main = "Residuals vs Fitted (Original)", 
     xlab = "Fitted Values", ylab = "Residuals",
     pch = 16, cex = 0.5)
abline(h = 0, col = "red", lwd = 2)

qqnorm(resid(model_original), main = "Normal Q-Q Plot (Original)",
       pch = 16, cex = 0.5)
qqline(resid(model_original), col = "red", lwd = 2)

hist(resid(model_original), breaks = 50, 
     main = "Residuals Distribution (Original)", 
     xlab = "Residuals", col = "lightblue")

# Selected transformation plots
plot(fitted(selected_model), resid(selected_model), 
     main = paste("Residuals vs Fitted (", best_transformation, ")", sep=""), 
     xlab = "Fitted Values", ylab = "Residuals",
     pch = 16, cex = 0.5)
abline(h = 0, col = "red", lwd = 2)

qqnorm(resid(selected_model), main = paste("Normal Q-Q Plot (", best_transformation, ")", sep=""),
       pch = 16, cex = 0.5)
qqline(resid(selected_model), col = "red", lwd = 2)

hist(resid(selected_model), breaks = 50, 
     main = paste("Residuals Distribution (", best_transformation, ")", sep=""), 
     xlab = "Residuals", col = "lightgreen")

dev.off()

# ============================================================================
# COMPREHENSIVE MODEL ANALYSIS FUNCTION
# ============================================================================

analyze_model_comprehensive <- function(formula, model_name, data, variable_used, transformation_summary) {
  
  cat("\n=== ANALYZING", model_name, "===\n")
  
  # Fit model
  model <- lmer(formula, data = data, REML = FALSE, control = control_opts)
  
  # Model diagnostics
  residuals_sample <- sample(resid(model), min(5000, length(resid(model))))
  residual_normality <- shapiro.test(residuals_sample)
  
  # Random effects normality
  random_effects <- ranef(model)
  re_vals <- unlist(random_effects)
  if(length(re_vals) >= 3) {
    re_normality <- shapiro.test(re_vals)
    re_norm_p <- re_normality$p.value
  } else {
    re_norm_p <- NA
  }
  
  # Calculate ICC
  var_comp <- VarCorr(model)
  var_df <- as.data.frame(var_comp)
  if(nrow(var_df) > 1) {
    random_var <- var_df$vcov[1]
    residual_var <- var_df$vcov[nrow(var_df)]
    icc <- random_var / (random_var + residual_var)
  } else {
    icc <- 0
  }
  
  # ANOVA Type III tests
  anova_result <- Anova(model, type = "III")
  
  # Extract significance levels
  get_significance <- function(p_val) {
    if(is.na(p_val)) return("--")
    if(p_val < 0.001) return("***")
    if(p_val < 0.01) return("**")
    if(p_val < 0.05) return("*")
    return("NS")
  }
  
  treatment_p <- if("treatment" %in% rownames(anova_result)) anova_result["treatment", "Pr(>Chisq)"] else NA
  thermal_phase_p <- if("thermal_phase" %in% rownames(anova_result)) anova_result["thermal_phase", "Pr(>Chisq)"] else NA
  time_block_p <- if("time_block" %in% rownames(anova_result)) anova_result["time_block", "Pr(>Chisq)"] else NA
  
  # Find interaction effects
  interaction_effects <- grep(":", rownames(anova_result), value = TRUE)
  if(length(interaction_effects) > 0) {
    interaction_p <- anova_result[interaction_effects[1], "Pr(>Chisq)"]
  } else {
    interaction_p <- NA
  }
  
  # Post-hoc analysis
  formula_str <- paste(deparse(formula), collapse = " ")
  
  if(grepl("thermal_phase.*time_block|time_block.*thermal_phase", formula_str)) {
    emm <- emmeans(model, ~ treatment | thermal_phase * time_block, mode = "asymptotic")
    grouping_structure <- "thermal_phase * time_block"
  } else if(grepl("thermal_phase", formula_str)) {
    emm <- emmeans(model, ~ treatment | thermal_phase, mode = "asymptotic")
    grouping_structure <- "thermal_phase"
  } else if(grepl("time_block", formula_str)) {
    emm <- emmeans(model, ~ treatment | time_block, mode = "asymptotic")
    grouping_structure <- "time_block"
  } else {
    emm <- emmeans(model, ~ treatment, mode = "asymptotic")
    grouping_structure <- "none"
  }
  
  cld_result <- cld(emm, Letters = letters, decreasing = TRUE)
  emm_summary <- summary(emm)
  
  # Handle confidence intervals
  if("lower.CL" %in% colnames(emm_summary)) {
    lower_col <- "lower.CL"
    upper_col <- "upper.CL"
  } else if("asymp.LCL" %in% colnames(emm_summary)) {
    lower_col <- "asymp.LCL"
    upper_col <- "asymp.UCL"
  } else {
    emm_summary$lower.CL <- emm_summary$emmean - 1.96 * emm_summary$SE
    emm_summary$upper.CL <- emm_summary$emmean + 1.96 * emm_summary$SE
    lower_col <- "lower.CL"
    upper_col <- "upper.CL"
  }
  
  # Create comprehensive results table
  results_table <- data.frame(
    Model = model_name,
    Treatment = emm_summary$treatment,
    Thermal_Phase = if("thermal_phase" %in% colnames(emm_summary)) emm_summary$thermal_phase else "All_Combined",
    Time_Block = if("time_block" %in% colnames(emm_summary)) emm_summary$time_block else "All_Combined",
    EMM = round(emm_summary$emmean, 4),
    SE = round(emm_summary$SE, 4),
    Lower_CI = round(emm_summary[[lower_col]], 4),
    Upper_CI = round(emm_summary[[upper_col]], 4),
    Statistical_Group = as.character(cld_result$.group),
    
    # Transformation comparison (unified assessment)
    Original_AIC = transformation_summary$AIC[transformation_summary$Transformation == "Original"],
    Log_AIC = transformation_summary$AIC[transformation_summary$Transformation == "Log"],
    SQRT_AIC = transformation_summary$AIC[transformation_summary$Transformation == "SQRT"],
    CBRT_AIC = transformation_summary$AIC[transformation_summary$Transformation == "CBRT"],
    
    # Model diagnostics
    Residual_Normality_P = format(residual_normality$p.value, scientific = TRUE, digits = 2),
    Random_Effect_Normality_P = if(is.na(re_norm_p)) "NA" else format(re_norm_p, scientific = TRUE, digits = 2),
    ICC = format(icc, scientific = TRUE, digits = 2),
    
    # Type III test significance
    Treatment_Sig = get_significance(treatment_p),
    Thermal_Phase_Sig = get_significance(thermal_phase_p),
    Time_Block_Sig = get_significance(time_block_p),
    Interaction_Sig = get_significance(interaction_p),
    
    # Model info
    Variable_Used = variable_used,
    Grouping_Structure = grouping_structure,
    Model_AIC = round(AIC(model), 2),
    
    stringsAsFactors = FALSE
  )
  
  return(results_table)
}

# ============================================================================
# THREE MODELING APPROACHES
# ============================================================================

# Build model formulas using selected variable
model_formula_1 <- as.formula(paste(selected_variable, "~ treatment * thermal_phase + (1|tube_id) + (1|date_factor)"))
model_formula_2 <- as.formula(paste(selected_variable, "~ treatment * time_block + (1|tube_id) + (1|date_factor)"))
model_formula_3 <- as.formula(paste(selected_variable, "~ treatment * thermal_phase * time_block + (1|tube_id) + (1|date_factor)"))

# Run all three models with comprehensive diagnostics
results_1 <- analyze_model_comprehensive(model_formula_1, "Model_1_Thermal_Phase", df, selected_description, transformation_summary)
results_2 <- analyze_model_comprehensive(model_formula_2, "Model_2_Time_Block", df, selected_description, transformation_summary)
results_3 <- analyze_model_comprehensive(model_formula_3, "Model_3_Full_Factorial", df, selected_description, transformation_summary)

# ============================================================================
# COMBINE RESULTS AND SAVE
# ============================================================================

cat("\n=== CREATING COMBINED RESULTS ===\n")

# Combine all model results
combined_results <- rbind(results_1, results_2, results_3)

# Add model comparison information
combined_results$Best_Model <- ifelse(combined_results$Model_AIC == min(combined_results$Model_AIC), "YES", "NO")
combined_results$AIC_Rank <- match(combined_results$Model_AIC, sort(unique(combined_results$Model_AIC)))
combined_results$Delta_AIC <- combined_results$Model_AIC - min(combined_results$Model_AIC)

# Save the comprehensive combined results
write.csv(combined_results, "Combined_LMM_Delta_Ta_Complete_Results.csv", row.names = FALSE)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Data transformation used:", selected_description, "\n")
cat("Total observations:", nrow(df), "\n")
cat("Best model (lowest AIC):", combined_results$Model[combined_results$AIC_Rank == 1][1], "\n")

cat("\nFiles generated:\n")
cat("1. Combined_LMM_Delta_Ta_Complete_Results.csv - Complete results with diagnostics\n")
cat("2. Diagnostic_Plots_Original_vs_Selected.png - Model diagnostic comparison\n")

cat("\nCombined results include:\n")
cat("- EMM, SE, confidence intervals, and statistical groups for all models\n")
cat("- AIC values for all transformations (unified assessment)\n")
cat("- Residual and random effects normality tests\n")
cat("- ICC (Intraclass Correlation Coefficient)\n")
cat("- Type III test significance levels\n")
cat("- Model ranking based on AIC\n")