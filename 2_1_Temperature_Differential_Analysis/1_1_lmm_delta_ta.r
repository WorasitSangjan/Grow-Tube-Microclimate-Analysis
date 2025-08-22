# ============================================================================
# LINEAR MIXED-EFFECTS MODEL ANALYSIS FOR GROW TUBE TEMPERATURE EFFECTS
# Analysis of delta_Ta (tube temperature - air temperature)
# 
# WORKFLOW:
# 1. Data preparation
# 2. Normality assessment and transformation selection
# 3. Three modeling approaches using selected data
# 4. Model comparison and results
# ============================================================================

# Load required libraries
library(lme4)
library(emmeans)
library(car)
library(multcomp)

# ============================================================================
# STEP 1: DATA LOADING AND PREPARATION
# ============================================================================

# Load data
df <- read.csv("grow_tube_categorized.csv")
cat("Original dataset:", nrow(df), "observations\n")

# Remove control group
df <- df[df$Material != "Uncover", ]
cat("After removing control:", nrow(df), "observations\n")

# Create factor variables
df$treatment <- as.factor(paste(df$Material, df$Condition, sep="_"))
df$tube_id <- as.factor(paste(df$treatment, "Rep", df$Replicate, sep="_"))
df$date_factor <- as.factor(df$date)
df$thermal_phase <- as.factor(df$thermal_phase)
df$time_block <- as.factor(df$time_block)

# Check data
cat("delta_Ta range:", range(df$delta_Ta, na.rm = TRUE), "\n")
cat("Missing values in delta_Ta:", sum(is.na(df$delta_Ta)), "\n")
cat("Treatments:", levels(df$treatment), "\n")

# ============================================================================
# STEP 2: NORMALITY ASSESSMENT AND DATA TRANSFORMATION SELECTION
# ============================================================================

# Model control parameters
control_opts <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5),
  calc.derivs = FALSE,
  check.nobs.vs.nlev = "ignore",
  check.nobs.vs.rankZ = "ignore"
)

# Test Model A: Original scale
model_original <- lmer(delta_Ta ~ treatment * thermal_phase + (1|tube_id) + (1|date_factor), 
                       data = df, REML = FALSE, control = control_opts)

# Create residual diagnostic plots for original model
png("residuals_original_delta_Ta.png", width = 3600, height = 1200, res = 250)
par(mfrow = c(1, 3), 
    mar = c(5, 5, 4, 2),
    cex.main = 2.5,
    cex.lab = 2.2,
    cex.axis = 1.9,
    font.main = 2,
    font.lab = 2,
    lwd = 2)

plot(fitted(model_original), resid(model_original), 
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", ylab = "Residuals",
     pch = 16, cex = 0.8)
abline(h = 0, col = "red", lwd = 3)

qqnorm(resid(model_original), main = "Normal Q-Q Plot",
       pch = 16, cex = 0.8)
qqline(resid(model_original), col = "red", lwd = 3)

hist(resid(model_original), breaks = 50, 
     main = "Residuals Distribution", 
     xlab = "Residuals", col = "lightblue")
dev.off()

# Test normality of original model residuals
residuals_orig <- sample(resid(model_original), min(5000, length(resid(model_original))))
shapiro_orig <- shapiro.test(residuals_orig)
cat("Original model AIC:", round(AIC(model_original), 2), "\n")
cat("Original model Shapiro p-value:", format.pval(shapiro_orig$p.value), "\n")

# Test Model B: Cube root transformed scale
df$delta_Ta_trans <- sign(df$delta_Ta) * abs(df$delta_Ta)^(1/3)

model_transformed <- lmer(delta_Ta_trans ~ treatment * thermal_phase + (1|tube_id) + (1|date_factor), 
                          data = df, REML = FALSE, control = control_opts)

# Create residual diagnostic plots for transformed model
png("residuals_transformed_delta_Ta.png", width = 3600, height = 1200, res = 250)
par(mfrow = c(1, 3), 
    mar = c(5, 5, 4, 2),
    cex.main = 2.5,
    cex.lab = 2.2,
    cex.axis = 1.9,
    font.main = 2,
    font.lab = 2,
    lwd = 2)

plot(fitted(model_transformed), resid(model_transformed), 
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", ylab = "Residuals",
     pch = 16, cex = 0.8)
abline(h = 0, col = "red", lwd = 3)

qqnorm(resid(model_transformed), main = "Normal Q-Q Plot",
       pch = 16, cex = 0.8)
qqline(resid(model_transformed), col = "red", lwd = 3)

hist(resid(model_transformed), breaks = 50, 
     main = "Residuals Distribution", 
     xlab = "Residuals", col = "lightgreen")
dev.off()

# Test normality of transformed model residuals
residuals_trans <- sample(resid(model_transformed), min(5000, length(resid(model_transformed))))
shapiro_trans <- shapiro.test(residuals_trans)
cat("Transformed model AIC:", round(AIC(model_transformed), 2), "\n")
cat("Transformed model Shapiro p-value:", format.pval(shapiro_trans$p.value), "\n")

# DECISION: Select the best data transformation
cat("\n=== DATA TRANSFORMATION SELECTION ===\n")
cat("Original model     - AIC:", round(AIC(model_original), 2), "| Shapiro p-value:", format.pval(shapiro_orig$p.value), "\n")
cat("Transformed model  - AIC:", round(AIC(model_transformed), 2), "| Shapiro p-value:", format.pval(shapiro_trans$p.value), "\n")

# Selection criteria: prefer normality, then AIC
if(shapiro_trans$p.value > shapiro_orig$p.value && shapiro_trans$p.value > 0.001) {
  selected_variable <- "delta_Ta_trans"
  selected_description <- "cube root transformed delta_Ta"
  cat("*** SELECTED: Transformed data (better normality) ***\n")
} else if(shapiro_orig$p.value > 0.05) {
  selected_variable <- "delta_Ta"
  selected_description <- "original delta_Ta"
  cat("*** SELECTED: Original data (adequate normality) ***\n")
} else {
  if(AIC(model_transformed) < AIC(model_original)) {
    selected_variable <- "delta_Ta_trans"
    selected_description <- "cube root transformed delta_Ta"
    cat("*** SELECTED: Transformed data (lower AIC) ***\n")
  } else {
    selected_variable <- "delta_Ta"
    selected_description <- "original delta_Ta"
    cat("*** SELECTED: Original data (lower AIC) ***\n")
  }
}

# ============================================================================
# STEP 3: THREE MODELING APPROACHES USING SELECTED DATA
# ============================================================================

# Build the model formula using selected variable
model_formula_1 <- as.formula(paste(selected_variable, "~ treatment * thermal_phase + (1|tube_id) + (1|date_factor)"))
model_formula_2 <- as.formula(paste(selected_variable, "~ treatment * time_block + (1|tube_id) + (1|date_factor)"))
model_formula_3 <- as.formula(paste(selected_variable, "~ treatment * thermal_phase * time_block + (1|tube_id) + (1|date_factor)"))

# ============================================================================
# MODEL 1: TREATMENT × THERMAL_PHASE
# ============================================================================

cat("\n=== MODEL 1: TREATMENT x THERMAL_PHASE ===\n")
model_1 <- lmer(model_formula_1, data = df, REML = FALSE, control = control_opts)
cat("Model 1 AIC:", round(AIC(model_1), 2), "\n")

anova_1 <- Anova(model_1, type = "III")
print(anova_1)

# Post-hoc analysis
emm_1 <- emmeans(model_1, ~ treatment | thermal_phase, mode = "asymptotic")
print(emm_1)

pairs_1 <- pairs(emm_1, adjust = "tukey")
print(pairs_1)

cld_1 <- cld(emm_1, Letters = letters, decreasing = TRUE)
print(cld_1)

# ============================================================================
# MODEL 2: TREATMENT × TIME_BLOCK
# ============================================================================

cat("\n=== MODEL 2: TREATMENT x TIME_BLOCK ===\n")
model_2 <- lmer(model_formula_2, data = df, REML = FALSE, control = control_opts)
cat("Model 2 AIC:", round(AIC(model_2), 2), "\n")

anova_2 <- Anova(model_2, type = "III")
print(anova_2)

# Post-hoc analysis
emm_2 <- emmeans(model_2, ~ treatment | time_block, mode = "asymptotic")
print(emm_2)

pairs_2 <- pairs(emm_2, adjust = "tukey")
print(pairs_2)

cld_2 <- cld(emm_2, Letters = letters, decreasing = TRUE)
print(cld_2)

# ============================================================================
# MODEL 3: FULL FACTORIAL
# ============================================================================

cat("\n=== MODEL 3: FULL FACTORIAL ===\n")
model_3 <- lmer(model_formula_3, data = df, REML = FALSE, control = control_opts)
cat("Model 3 AIC:", round(AIC(model_3), 2), "\n")

anova_3 <- Anova(model_3, type = "III")
print(anova_3)

# Check for significant three-way interaction
three_way_p <- anova_3["treatment:thermal_phase:time_block", "Pr(>Chisq)"]
cat("Three-way interaction p-value:", format.pval(three_way_p), "\n")

if(three_way_p < 0.05) {
  cat("*** THREE-WAY INTERACTION IS SIGNIFICANT ***\n")
  emm_3 <- emmeans(model_3, ~ treatment | thermal_phase * time_block, mode = "asymptotic")
} else {
  emm_3 <- emmeans(model_3, ~ treatment | thermal_phase, mode = "asymptotic")
}

print(emm_3)
pairs_3 <- pairs(emm_3, adjust = "tukey")
print(pairs_3)

# ============================================================================
# MODEL COMPARISON
# ============================================================================

cat("\n=== MODEL COMPARISON ===\n")

model_comparison <- data.frame(
  Model = c("Model_1_Thermal_Phase", "Model_2_Time_Block", "Model_3_Full_Factorial"),
  Formula = c("treatment × thermal_phase", "treatment x time_block", "treatment x thermal_phase x time_block"),
  AIC = c(AIC(model_1), AIC(model_2), AIC(model_3)),
  Purpose = c("Weather-based effects", "Daily timing effects", "Comprehensive interactions")
)
model_comparison$AIC_Rank <- rank(model_comparison$AIC, ties.method = "min")
model_comparison$Delta_AIC <- model_comparison$AIC - min(model_comparison$AIC)

print(model_comparison)

best_model_name <- model_comparison$Model[which.min(model_comparison$AIC)]
cat("Best model (lowest AIC):", best_model_name, "\n")

# ============================================================================
# GENERATE OUTPUT FILES
# ============================================================================

# Function to create summary table
create_summary_table <- function(emm_result, model_name) {
  emm_summary <- summary(emm_result)
  
  # Handle confidence interval column names
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
  
  # Create base table
  result_table <- data.frame(
    Model = model_name,
    Treatment = emm_summary$treatment,
    Thermal_Phase = if("thermal_phase" %in% colnames(emm_summary)) emm_summary$thermal_phase else "All_Combined",
    Time_Block = if("time_block" %in% colnames(emm_summary)) emm_summary$time_block else "All_Combined",
    EMM = round(emm_summary$emmean, 4),
    SE = round(emm_summary$SE, 4),
    Lower_CI = round(emm_summary[[lower_col]], 4),
    Upper_CI = round(emm_summary[[upper_col]], 4),
    Variable_Used = selected_description,
    stringsAsFactors = FALSE
  )
  
  return(result_table)
}

# Generate summary tables
table_1 <- create_summary_table(emm_1, "Model_1_Thermal_Phase")
table_2 <- create_summary_table(emm_2, "Model_2_Time_Block") 
table_3 <- create_summary_table(emm_3, "Model_3_Full_Factorial")

# Save individual model results
write.csv(table_1, "Model_1_Treatment_x_Thermal_Phase_deltaTA.csv", row.names = FALSE)
write.csv(table_2, "Model_2_Treatment_x_Time_Block_deltaTA.csv", row.names = FALSE)
write.csv(table_3, "Model_3_Full_Factorial_deltaTA.csv", row.names = FALSE)

# Save model comparison
write.csv(model_comparison, "Model_Comparison_deltaTA.csv", row.names = FALSE)

# Create overall treatment ranking
treatments <- levels(df$treatment)
treatment_summary <- data.frame(
  Treatment = treatments,
  Model_1_Mean = sapply(treatments, function(t) mean(table_1$EMM[table_1$Treatment == t], na.rm = TRUE)),
  Model_2_Mean = sapply(treatments, function(t) mean(table_2$EMM[table_2$Treatment == t], na.rm = TRUE)),
  Model_3_Mean = sapply(treatments, function(t) mean(table_3$EMM[table_3$Treatment == t], na.rm = TRUE))
)
treatment_summary$Best_Model_Rank <- rank(-treatment_summary[[paste0("Model_", which.min(model_comparison$AIC), "_Mean")]], ties.method = "min")

write.csv(treatment_summary, "Treatment_Summary_deltaTA.csv", row.names = FALSE)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Data transformation used:", selected_description, "\n")
cat("Total observations:", nrow(df), "\n")
cat("Best model:", best_model_name, "\n")
cat("Three-way interaction significant:", ifelse(three_way_p < 0.05, "YES", "NO"), "\n")

cat("\nFiles generated:\n")
cat("- residuals_original_delta_Ta.png\n")
cat("- residuals_transformed_delta_Ta.png\n")
cat("- Model_1_Treatment_x_Thermal_Phase_deltaTA.csv\n")
cat("- Model_2_Treatment_x_Time_Block_deltaTA.csv\n")
cat("- Model_3_Full_Factorial_deltaTA.csv\n")
cat("- Model_Comparison_deltaTA.csv\n")
cat("- Treatment_Summary_deltaTA.csv\n")