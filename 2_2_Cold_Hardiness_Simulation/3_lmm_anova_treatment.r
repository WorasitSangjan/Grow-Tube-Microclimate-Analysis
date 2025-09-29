# Author: Worasit Sangjan
# Date: September 29, 2025
# Version 2.2

# ============================================================================
# COLD HARDINESS ANALYSIS - BY TREATMENT
# Compare varieties within each treatment
# ============================================================================

# ============================================================================
# USER CONFIGURATION
# ============================================================================
SEASON_TO_ANALYZE <- "2024-2025"
ANALYSIS_TYPE <- "by_treatment"

cat("===============================================================\n")
cat("COLD HARDINESS ANALYSIS - BY TREATMENT\n")
cat("Season:", SEASON_TO_ANALYZE, "\n")
cat("===============================================================\n")

# ============================================================================
# LOAD LIBRARIES
# ============================================================================
required_packages <- c("lme4", "emmeans", "car", "ggplot2", "dplyr", "gridExtra", 
                      "RColorBrewer", "grid", "multcomp", "agricolae")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================
cat("\n=== LOADING DATA ===\n")

dynamic_file <- paste0("dynamic_slopes_", ANALYSIS_TYPE, ".csv")
winter_file <- paste0("winter_means_", ANALYSIS_TYPE, ".csv")
df_dynamic <- read.csv(dynamic_file)
df_winter <- read.csv(winter_file)

df_dynamic_season <- df_dynamic[df_dynamic$Season == SEASON_TO_ANALYZE, ]
df_winter_season <- df_winter[df_winter$Season == SEASON_TO_ANALYZE, ]

df_dynamic_season$Treatment <- gsub("Uncover_uncover", "No-Tube", df_dynamic_season$Treatment)
df_dynamic_season$Treatment <- gsub("Paper_raise", "Paper_high", df_dynamic_season$Treatment)
df_dynamic_season$Treatment <- gsub("Plastic_raise", "Plastic_high", df_dynamic_season$Treatment)

df_winter_season$Treatment <- gsub("Uncover_uncover", "No-Tube", df_winter_season$Treatment)
df_winter_season$Treatment <- gsub("Paper_raise", "Paper_high", df_winter_season$Treatment)
df_winter_season$Treatment <- gsub("Plastic_raise", "Plastic_high", df_winter_season$Treatment)

treatment_order <- c("Paper_low", "Paper_high", "No-Tube", "Plastic_low", "Plastic_high")
variety_order <- c("Cabernet Sauvignon", "Chardonnay", "Concord", "Mourvedre")

factor_cols <- c("Season", "Variety", "Treatment", "Phase", "Replicate")
for(col in factor_cols) {
  if(col %in% colnames(df_dynamic_season)) {
    if(col == "Treatment") {
      df_dynamic_season[[col]] <- factor(df_dynamic_season[[col]], levels = treatment_order)
    } else if(col == "Variety") {
      df_dynamic_season[[col]] <- factor(df_dynamic_season[[col]], levels = variety_order)
    } else {
      df_dynamic_season[[col]] <- as.factor(df_dynamic_season[[col]])
    }
  }
  if(col %in% colnames(df_winter_season)) {
    if(col == "Treatment") {
      df_winter_season[[col]] <- factor(df_winter_season[[col]], levels = treatment_order)
    } else if(col == "Variety") {
      df_winter_season[[col]] <- factor(df_winter_season[[col]], levels = variety_order)
    } else {
      df_winter_season[[col]] <- as.factor(df_winter_season[[col]])
    }
  }
}

df_dynamic_season$Vine_ID <- with(df_dynamic_season, paste(Variety, Treatment, Replicate, sep = "_"))
results_dir <- paste0("results_", SEASON_TO_ANALYZE, "_by_treatment")
plots_dir <- paste0("plots_", SEASON_TO_ANALYZE, "_by_treatment")

if(!dir.exists(results_dir)) dir.create(results_dir)
if(!dir.exists(plots_dir)) dir.create(plots_dir)

cat("Dynamic data loaded:", nrow(df_dynamic_season), "slope observations\n")
cat("Winter data loaded:", nrow(df_winter_season), "mean hardiness observations\n")

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

test_transformations <- function(treatment_data, control_opts) {
  transformation_results <- list()
  
  model_original <- lmer(Slope_C_per_day ~ Variety * Phase + (1|Replicate), 
                        data = treatment_data, REML = FALSE, control = control_opts)
  
  min_slope <- min(treatment_data$Slope_C_per_day, na.rm = TRUE)
  shift_value <- abs(min_slope) + 0.001
  treatment_data$Slope_log <- log(treatment_data$Slope_C_per_day + shift_value)
  
  model_log <- tryCatch({
    lmer(Slope_log ~ Variety * Phase + (1|Replicate), 
         data = treatment_data, REML = FALSE, control = control_opts)
  }, error = function(e) NULL)
  
  treatment_data$Slope_sqrt <- sqrt(treatment_data$Slope_C_per_day + shift_value)
  
  model_sqrt <- tryCatch({
    lmer(Slope_sqrt ~ Variety * Phase + (1|Replicate), 
         data = treatment_data, REML = FALSE, control = control_opts)
  }, error = function(e) NULL)
  
  treatment_data$Slope_cbrt <- sign(treatment_data$Slope_C_per_day) * abs(treatment_data$Slope_C_per_day)^(1/3)
  
  model_cbrt <- tryCatch({
    lmer(Slope_cbrt ~ Variety * Phase + (1|Replicate), 
         data = treatment_data, REML = FALSE, control = control_opts)
  }, error = function(e) NULL)
  
  transformation_results$original <- list(
    model = model_original,
    aic = AIC(model_original),
    shapiro_p = shapiro.test(residuals(model_original))$p.value
  )
  
  if(!is.null(model_log)) {
    transformation_results$log <- list(
      model = model_log,
      aic = AIC(model_log),
      shapiro_p = shapiro.test(residuals(model_log))$p.value,
      shift_value = shift_value
    )
  }
  
  if(!is.null(model_sqrt)) {
    transformation_results$sqrt <- list(
      model = model_sqrt,
      aic = AIC(model_sqrt),
      shapiro_p = shapiro.test(residuals(model_sqrt))$p.value,
      shift_value = shift_value
    )
  }
  
  if(!is.null(model_cbrt)) {
    transformation_results$cbrt <- list(
      model = model_cbrt,
      aic = AIC(model_cbrt),
      shapiro_p = shapiro.test(residuals(model_cbrt))$p.value
    )
  }
  
  return(transformation_results)
}

save_diagnostics <- function(treatment, results, transformation_results, file_path) {
  sink(file_path)
  
  cat("=== LMM DIAGNOSTICS FOR", toupper(treatment), "DYNAMIC ANALYSIS ===\n\n")
  
  cat("TRANSFORMATION COMPARISON (AIC and Normality):\n")
  cat("Model\t\tAIC\t\tShapiro p-value\tInterpretation\n")
  cat("----------------------------------------------------------------\n")
  
  for(trans_name in names(transformation_results)) {
    trans_result <- transformation_results[[trans_name]]
    interpretation <- ifelse(trans_result$shapiro_p > 0.05, "NORMAL", "NON-NORMAL")
    cat(sprintf("%-12s\t%.2f\t\t%.2e\t%s\n", 
                toupper(trans_name), trans_result$aic, trans_result$shapiro_p, interpretation))
  }
  
  best_aic <- min(sapply(transformation_results, function(x) x$aic))
  best_trans <- names(transformation_results)[which.min(sapply(transformation_results, function(x) x$aic))]
  
  cat("\nBEST MODEL BY AIC:", toupper(best_trans), "(AIC =", best_aic, ")\n")
  cat("USING ORIGINAL DATA FOR BIOLOGICAL INTERPRETATION\n\n")
  
  cat("MODEL SUMMARY (ORIGINAL DATA):\n")
  print(summary(results$model))
  
  cat("\n\nANOVA TYPE III TESTS:\n")
  print(results$anova)
  
  cat("\n\nLMM ASSUMPTION CHECKS:\n")
  cat("1. NORMALITY OF RESIDUALS:\n")
  cat("   Shapiro-Wilk: W =", results$shapiro_residuals$statistic, ", p =", results$shapiro_residuals$p.value, "\n")
  cat("   Interpretation:", ifelse(results$shapiro_residuals$p.value > 0.05, "NORMAL", "NON-NORMAL"), "\n\n")
  
  cat("2. NORMALITY OF RANDOM EFFECTS:\n")
  cat("   Shapiro-Wilk: W =", results$shapiro_random$statistic, ", p =", results$shapiro_random$p.value, "\n")
  cat("   Interpretation:", ifelse(results$shapiro_random$p.value > 0.05, "NORMAL", "NON-NORMAL"), "\n\n")
  
  sink()
}

# ============================================================================
# RUN ANALYSIS
# ============================================================================
cat("\n=== RUNNING ANALYSIS ===\n")

control_opts <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5),
  calc.derivs = FALSE
)

dynamic_results <- list()
winter_results <- list()
treatments <- levels(df_dynamic_season$Treatment)

# Dynamic analysis
for(treatment in treatments) {
  cat("Analyzing dynamic phases for", treatment, "\n")
  treatment_data <- df_dynamic_season[df_dynamic_season$Treatment == treatment, ]
  
  if(nrow(treatment_data) == 0) {
    cat("  No data for treatment:", treatment, "\n")
    next
  }
  
  transformation_results <- test_transformations(treatment_data, control_opts)
  
  model <- transformation_results$original$model
  emm <- emmeans(model, ~ Variety | Phase)
  pairs_result <- pairs(emm, adjust = "tukey")
  cld_result <- cld(emm, Letters = letters, decreasing = TRUE)
  
  model_residuals <- residuals(model)
  fitted_values <- fitted(model)
  random_effects <- ranef(model)$Replicate[,1]
  shapiro_residuals <- shapiro.test(model_residuals)
  shapiro_random <- shapiro.test(random_effects)
  
  dynamic_results[[treatment]] <- list(
    model = model,
    emm = emm,
    pairs = pairs_result,
    cld = cld_result,
    anova = Anova(model, type = "III"),
    residuals = model_residuals,
    fitted_values = fitted_values,
    random_effects = random_effects,
    shapiro_residuals = shapiro_residuals,
    shapiro_random = shapiro_random
  )
  
  lmm_diagnostic_file <- paste0(results_dir, "/LMM_Diagnostics_", gsub("_", "-", treatment), ".txt")
  save_diagnostics(treatment, dynamic_results[[treatment]], transformation_results, lmm_diagnostic_file)
}

# Winter analysis with HSD groups
for(treatment in treatments) {
  cat("Analyzing winter maintenance for", treatment, "\n")
  treatment_winter <- df_winter_season[df_winter_season$Treatment == treatment, ]
  
  if(nrow(treatment_winter) == 0) {
    cat("  No winter data for treatment:", treatment, "\n")
    next
  }
  
  anova_model <- aov(Mean_Hardiness ~ Variety, data = treatment_winter)
  
  emm_winter <- emmeans(anova_model, ~ Variety)
  pairs_winter <- pairs(emm_winter, adjust = "tukey")
  cld_winter <- cld(emm_winter, Letters = letters, decreasing = TRUE)
  
  tukey_result <- TukeyHSD(anova_model)
  hsd_groups <- HSD.test(anova_model, "Variety")
  
  descriptive_stats <- treatment_winter %>%
    group_by(Variety) %>%
    summarise(
      n = n(),
      mean = mean(Mean_Hardiness),
      sd = sd(Mean_Hardiness),
      se = sd(Mean_Hardiness) / sqrt(n()),
      .groups = 'drop'
    )
  
  residuals_anova <- residuals(anova_model)
  shapiro_test <- shapiro.test(residuals_anova)
  levene_test <- car::leveneTest(Mean_Hardiness ~ Variety, data = treatment_winter)
  bartlett_test <- bartlett.test(Mean_Hardiness ~ Variety, data = treatment_winter)
  
  winter_results[[treatment]] <- list(
    model = anova_model,
    emm = emm_winter,
    pairs = pairs_winter,
    cld = cld_winter,
    summary = summary(anova_model),
    tukey = tukey_result,
    hsd_groups = hsd_groups,
    means = aggregate(Mean_Hardiness ~ Variety, treatment_winter, mean),
    descriptive_stats = descriptive_stats,
    raw_data = treatment_winter,
    residuals = residuals_anova,
    shapiro_test = shapiro_test,
    levene_test = levene_test,
    bartlett_test = bartlett_test
  )
  
  diagnostic_file <- paste0(results_dir, "/Winter_Diagnostics_", gsub("_", "-", treatment), ".txt")
  sink(diagnostic_file)
  
  cat("=== WINTER ANALYSIS DIAGNOSTICS FOR", toupper(treatment), "===\n\n")
  
  cat("DESCRIPTIVE STATISTICS:\n")
  print(descriptive_stats)
  
  cat("\n\nANOVA ASSUMPTIONS:\n")
  cat("1. NORMALITY OF RESIDUALS (Shapiro-Wilk test):\n")
  cat("   W =", shapiro_test$statistic, ", p-value =", shapiro_test$p.value, "\n")
  cat("   Interpretation:", ifelse(shapiro_test$p.value > 0.05, "NORMAL (p > 0.05)", "NON-NORMAL (p <= 0.05)"), "\n\n")
  
  cat("2. HOMOGENEITY OF VARIANCE:\n")
  cat("   Levene Test: F =", levene_test[1,2], ", p-value =", levene_test[1,3], "\n")
  cat("   Interpretation:", ifelse(levene_test[1,3] > 0.05, "EQUAL VARIANCES (p > 0.05)", "UNEQUAL VARIANCES (p <= 0.05)"), "\n")
  cat("   Bartlett Test: Chi-sq =", bartlett_test$statistic, ", p-value =", bartlett_test$p.value, "\n")
  cat("   Interpretation:", ifelse(bartlett_test$p.value > 0.05, "EQUAL VARIANCES (p > 0.05)", "UNEQUAL VARIANCES (p <= 0.05)"), "\n\n")
  
  cat("3. RESIDUAL ANALYSIS:\n")
  cat("   Range of residuals:", range(residuals_anova)[1], "to", range(residuals_anova)[2], "\n")
  cat("   Mean of residuals:", mean(residuals_anova), "\n")
  cat("   SD of residuals:", sd(residuals_anova), "\n\n")
  
  cat("ANOVA RESULTS:\n")
  print(summary(anova_model))
  
  cat("\n\nTUKEY HSD P-VALUES:\n")
  print(tukey_result)
  
  cat("\n\nEMMEANS CLD GROUPINGS:\n")
  print(cld_winter)
  
  cat("\n\nHSD GROUPINGS (agricolae):\n")
  print(hsd_groups$groups)
  
  cat("\n\nHSD STATISTICS:\n")
  cat("Alpha level:", hsd_groups$alpha, "\n")
  cat("HSD value:", hsd_groups$HSD, "\n")
  cat("F-statistic:", hsd_groups$statistics$F, "\n")
  cat("P-value:", hsd_groups$statistics$P, "\n")
  
  cat("\n\nRAW DATA CHECK:\n")
  print(treatment_winter[, c("Variety", "Mean_Hardiness")])
  
  sink()
  cat("Saved diagnostics for", treatment, "to", diagnostic_file, "\n")
}

# ============================================================================
# CREATE VISUALIZATION
# ============================================================================
cat("\n=== CREATING VISUALIZATION ===\n")

variety_colors <- c(
  "Chardonnay" = "#1f77b4",
  "Cabernet Sauvignon" = "#2ca02c",
  "Concord" = "#ff7f0e",
  "Mourvedre" = "#d62728"
)

prepare_raw_data <- function() {
  all_raw_data <- data.frame()
  
  fall_data <- df_dynamic_season[df_dynamic_season$Phase %in% unique(df_dynamic_season$Phase)[grepl("fall", unique(df_dynamic_season$Phase), ignore.case = TRUE)], ]
  if(nrow(fall_data) > 0) {
    fall_processed <- fall_data %>%
      dplyr::select(Variety, Treatment, Slope_C_per_day, Replicate) %>%
      dplyr::mutate(
        Phase_Name = "Fall Acclimation",
        Y_Value = Slope_C_per_day,
        Y_Label = "Slope (°C/day)"
      )
    all_raw_data <- rbind(all_raw_data, fall_processed[, c("Variety", "Treatment", "Phase_Name", "Y_Value", "Y_Label", "Replicate")])
  }
  
  winter_processed <- df_winter_season %>%
    dplyr::select(Variety, Treatment, Mean_Hardiness, Replicate) %>%
    dplyr::mutate(
      Phase_Name = "Winter Maintenance",
      Y_Value = Mean_Hardiness,
      Y_Label = "Mean Hardiness (°C)"
    )
  all_raw_data <- rbind(all_raw_data, winter_processed[, c("Variety", "Treatment", "Phase_Name", "Y_Value", "Y_Label", "Replicate")])
  
  spring_data <- df_dynamic_season[df_dynamic_season$Phase %in% unique(df_dynamic_season$Phase)[grepl("spring", unique(df_dynamic_season$Phase), ignore.case = TRUE)], ]
  if(nrow(spring_data) > 0) {
    spring_processed <- spring_data %>%
      dplyr::select(Variety, Treatment, Slope_C_per_day, Replicate) %>%
      dplyr::mutate(
        Phase_Name = "Spring Deacclimation",
        Y_Value = Slope_C_per_day,
        Y_Label = "Slope (°C/day)"
      )
    all_raw_data <- rbind(all_raw_data, spring_processed[, c("Variety", "Treatment", "Phase_Name", "Y_Value", "Y_Label", "Replicate")])
  }
  
  all_raw_data$Phase_Name <- factor(all_raw_data$Phase_Name, 
                                   levels = c("Fall Acclimation", "Winter Maintenance", "Spring Deacclimation"))
  all_raw_data$Variety <- factor(all_raw_data$Variety, levels = variety_order)
  all_raw_data$Treatment <- factor(all_raw_data$Treatment, levels = treatment_order)
  
  return(all_raw_data)
}

get_statistical_letters <- function(phase_name) {
  letters_data <- data.frame()
  
  for(treatment in treatments) {
    if(phase_name == "Winter Maintenance") {
      if(treatment %in% names(winter_results)) {
        cld_df <- as.data.frame(winter_results[[treatment]]$cld)
        treatment_letters <- data.frame(
          Treatment = treatment,
          Variety = as.character(cld_df$Variety),
          .group = as.character(cld_df$.group),
          stringsAsFactors = FALSE
        )
        letters_data <- rbind(letters_data, treatment_letters)
      }
    } else {
      if(treatment %in% names(dynamic_results)) {
        cld_df <- as.data.frame(dynamic_results[[treatment]]$cld)
        
        if("Phase" %in% colnames(cld_df)) {
          if(phase_name == "Fall Acclimation") {
            phase_filter <- unique(cld_df$Phase)[grepl("fall", unique(cld_df$Phase), ignore.case = TRUE)][1]
          } else {
            phase_filter <- unique(cld_df$Phase)[grepl("spring", unique(cld_df$Phase), ignore.case = TRUE)][1]
          }
          
          cld_filtered <- cld_df[cld_df$Phase == phase_filter, ]
          
          if(nrow(cld_filtered) > 0) {
            treatment_letters <- data.frame(
              Treatment = treatment,
              Variety = as.character(cld_filtered$Variety),
              .group = as.character(cld_filtered$.group),
              stringsAsFactors = FALSE
            )
            letters_data <- rbind(letters_data, treatment_letters)
          }
        }
      }
    }
  }
  
  if(nrow(letters_data) == 0) {
    return(data.frame())
  }
  
  letters_data$Variety <- factor(letters_data$Variety, levels = variety_order)
  letters_data$Treatment <- factor(letters_data$Treatment, levels = treatment_order)
  
  n_varieties <- length(variety_order)
  variety_width <- 0.8
  
  letters_data$x_position <- NA
  
  for(i in 1:nrow(letters_data)) {
    treatment_idx <- which(treatment_order == as.character(letters_data$Treatment[i]))
    variety_idx <- which(variety_order == as.character(letters_data$Variety[i]))
    
    base_pos <- treatment_idx
    variety_offset <- (variety_idx - (n_varieties + 1)/2) * (variety_width / n_varieties)
    
    letters_data$x_position[i] <- base_pos + variety_offset
  }
  
  return(letters_data)
}

create_combined_plot <- function() {
  raw_data <- prepare_raw_data()
  phases <- c("Fall Acclimation", "Winter Maintenance", "Spring Deacclimation")
  y_labels <- c("Slope (°C/day)", "Mean hardiness (°C)", "Slope (°C/day)")
  
  plot_list <- list()
  
  for(i in 1:length(phases)) {
    phase_data <- raw_data[raw_data$Phase_Name == phases[i], ]
    letters_data <- get_statistical_letters(phases[i])
    
    max_values <- phase_data %>%
      group_by(Treatment, Variety) %>%
      summarise(max_val = max(Y_Value, na.rm = TRUE), .groups = "drop")
    
    letters_positioned <- data.frame()
    if(nrow(letters_data) > 0) {
      letters_positioned <- merge(letters_data, max_values, by = c("Treatment", "Variety"))
      letters_positioned$y_position <- letters_positioned$max_val * 1.08
    }
    
    p <- ggplot(phase_data, aes(x = Treatment, y = Y_Value, fill = Variety)) +
      geom_boxplot(position = position_dodge(width = 0.8), 
                   alpha = 0.7, outlier.shape = NA, width = 0.6) +
      geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
                 aes(color = Variety), size = 2, alpha = 0.9) +
      scale_fill_manual(values = variety_colors, name = "Variety",
                       labels = c("Cabernet Sauvignon", "Chardonnay", "Concord", "Mourvedre")) +
      scale_color_manual(values = variety_colors, name = "Variety",
                        labels = c("Cabernet Sauvignon", "Chardonnay", "Concord", "Mourvedre")) +
      scale_x_discrete(labels = c("Paper\nLow", "Paper\nHigh", "No-Tube", "Plastic\nLow", "Plastic\nHigh")) +
      labs(title = c("Fall acclimation", "Winter maintenance", "Spring deacclimation")[i],
           x = if(i == 2) "Treatment" else "",
           y = y_labels[i]) +
      theme_minimal(base_size = 16, base_family = "Arial") +
      theme(
        plot.title = element_text(size = 22, family = "Arial", hjust = 0.5, vjust = -1),
        axis.title = element_text(size = 21, family = "Arial"),
        axis.text.x = element_text(size = 17, family = "Arial"),
        axis.text.y = element_text(size = 17, family = "Arial"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
    
    if(nrow(letters_positioned) > 0) {
      p <- p + geom_text(data = letters_positioned, 
                aes(x = x_position, y = y_position, label = .group),
                size = 8, family = "Arial", vjust = 0, color = "black",
                inherit.aes = FALSE)
    }
    
    plot_list[[i]] <- p
  }
  
  legend_plot <- ggplot(raw_data, aes(x = Treatment, y = Y_Value, fill = Variety)) +
    geom_boxplot() +
    scale_fill_manual(values = variety_colors, name = "",
                     labels = c("Cabernet Sauvignon", "Chardonnay", "Concord", "Mourvedre")) +
    theme_minimal(base_family = "Arial") +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 20, family = "Arial"),
          legend.key.size = unit(1.2, "cm"),
          legend.margin = margin(t = 15),
          legend.box.spacing = unit(0.8, "cm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  legend <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]
  
  combined_plot <- grid.arrange(
    plot_list[[1]], plot_list[[2]], plot_list[[3]],
    legend,
    ncol = 3, nrow = 2,
    heights = c(4, 0.8),
    layout_matrix = rbind(c(1, 2, 3),
                         c(4, 4, 4))
  )
  
  return(combined_plot)
}

combined_plot <- create_combined_plot()

plot_filename <- paste0(plots_dir, "/Cold_Hardiness_", SEASON_TO_ANALYZE, ".png")
png(filename = plot_filename, width = 18, height = 6, units = "in", res = 600, bg = "white")
grid.draw(combined_plot)
dev.off()

cat("Plot saved as:", plot_filename, "\n")
cat("ANALYSIS COMPLETE\n")