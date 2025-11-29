# Author: Worasit Sangjan
# Date: 13 July 2025

# ============================================================================
# CHILLING UNITS AND GDD ANALYSIS
# Single combined graph with 3 panels: Chilling Units | GDD-like Heating | Traditional GDD
# ============================================================================

# ============================================================================
# USER CONFIGURATION
# ============================================================================
SEASON_TO_ANALYZE <- "2024_2025"  # Change to "2023_2024" or "2024_2025"

cat("===============================================================\n")
cat("CHILLING UNITS AND GDD ANALYSIS\n")
cat("Season:", SEASON_TO_ANALYZE, "\n")
cat("===============================================================\n")

# ============================================================================
# LOAD LIBRARIES
# ============================================================================
required_packages <- c("readr", "dplyr", "ggplot2", "car", "emmeans", 
                      "multcomp", "gridExtra", "grid")

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

# Set working directory
setwd(file.path("/Users/ws/Desktop/Paper_GrowthTube/Data_Analysis/1_data_analysis/6_chill_gdd", 
                paste0(SEASON_TO_ANALYZE, "_coldhardiness")))

# Define data structure
treatment_order <- c("Paper_low", "Paper_raise", "No-Tube", "Plastic_low", "Plastic_raise")
variety_order <- c("Cabernet Sauvignon", "Chardonnay", "Concord", "Mourvedre")

treatment_colors <- c(
  "Paper_low" = "#E74C3C",      # Red
  "Paper_raise" = "#F39C12",    # Orange  
  "No-Tube" = "#27AE60",        # Green
  "Plastic_low" = "#3498DB",    # Blue
  "Plastic_raise" = "#9B59B6"   # Purple
)

# Load data files
file_names <- c(
  "combined_detailed_results_Cabernet Sauvignon.csv",
  "combined_detailed_results_Chardonnay.csv", 
  "combined_detailed_results_Concord.csv",
  "combined_detailed_results_Mourvedre.csv"
)

all_data <- data.frame()
for(file in file_names) {
  if(file.exists(file)) {
    variety_data <- read_csv(file, show_col_types = FALSE)
    all_data <- rbind(all_data, variety_data)
    cat("Loaded:", file, "- Rows:", nrow(variety_data), "\n")
  } else {
    cat("WARNING: File not found:", file, "\n")
  }
}

# Calculate traditional GDD
all_data$T_daily_avg <- (all_data$T_max + all_data$T_min) / 2
all_data$Traditional_GDD_daily <- ifelse(all_data$T_daily_avg > 10, 
                                        all_data$T_daily_avg - 10, 0)

all_data <- all_data %>%
  group_by(Treatment, Replicate, Variety) %>%
  arrange(datetime) %>%
  mutate(Traditional_GDD_sum = cumsum(Traditional_GDD_daily)) %>%
  ungroup()

# Extract final values
final_metrics <- all_data %>%
  group_by(Treatment, Replicate, Variety) %>%
  summarise(
    Final_DD_Chilling = abs(min(DD_chilling_sum, na.rm = TRUE)),
    Final_DD_Heating = max(DD_heating_sum, na.rm = TRUE),
    Final_Traditional_GDD = max(Traditional_GDD_sum, na.rm = TRUE),
    .groups = "drop"
  )

# Clean treatment names and set factors
final_metrics$Treatment <- gsub("Uncover_uncover", "No-Tube", final_metrics$Treatment)
final_metrics$Treatment <- factor(final_metrics$Treatment, levels = treatment_order)
final_metrics$Variety <- factor(final_metrics$Variety, levels = variety_order)

# Create results directory
results_dir <- paste0("results_", SEASON_TO_ANALYZE, "_chill_gdd")
plots_dir <- paste0("plots_", SEASON_TO_ANALYZE, "_chill_gdd")

if(!dir.exists(results_dir)) dir.create(results_dir)
if(!dir.exists(plots_dir)) dir.create(plots_dir)

cat("Total combined data:", nrow(all_data), "observations\n")
cat("Final metrics dataset prepared with", nrow(final_metrics), "observations\n")

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

analyze_metric <- function(data, metric_column, metric_name, varieties) {
  cat("Analyzing", metric_name, "\n")
  all_results <- data.frame()
  
  for(variety in varieties) {
    cat("  Processing", variety, "\n")
    variety_data <- data[data$Variety == variety, ]
    if(nrow(variety_data) == 0) next
    
    # ANOVA model and assumption tests
    anova_model <- lm(as.formula(paste(metric_column, "~ Treatment")), data = variety_data)
    residuals_anova <- residuals(anova_model)
    shapiro_test <- shapiro.test(residuals_anova)
    levene_test <- car::leveneTest(variety_data[[metric_column]], variety_data$Treatment)
    bartlett_test <- bartlett.test(variety_data[[metric_column]], variety_data$Treatment)
    anova_result <- Anova(anova_model, type = "III")
    treatment_p <- anova_result[["Pr(>F)"]][1]
    
    # Post-hoc comparisons
    emm <- emmeans(anova_model, ~ Treatment)
    if(treatment_p < 0.05) {
      cld_result <- cld(emm, Letters = letters, decreasing = TRUE)
    } else {
      emm_summary <- summary(emm)
      cld_result <- data.frame(
        Treatment = emm_summary$Treatment,
        emmean = emm_summary$emmean,
        SE = emm_summary$SE,
        lower.CL = emm_summary$lower.CL,
        upper.CL = emm_summary$upper.CL,
        .group = rep("a", length(emm_summary$Treatment))
      )
    }
    
    # Descriptive statistics
    descriptive_stats <- variety_data %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        Mean = round(mean(.data[[metric_column]], na.rm = TRUE), 2),
        SD = round(sd(.data[[metric_column]], na.rm = TRUE), 2),
        SE = round(sd(.data[[metric_column]], na.rm = TRUE) / sqrt(n()), 2),
        .groups = "drop"
      )
    
    # Ensure alignment between descriptive stats and CLD results
    cld_df <- as.data.frame(cld_result)
    cld_df$Treatment <- as.character(cld_df$Treatment)
    descriptive_stats$Treatment <- as.character(descriptive_stats$Treatment)
    
    merged_results <- merge(descriptive_stats, cld_df[, c("Treatment", ".group")], 
                           by = "Treatment", all.x = TRUE)
    merged_results$.group[is.na(merged_results$.group)] <- "a"
    
    # Compile results
    variety_result <- data.frame(
      Metric = metric_name,
      Variety = variety,
      Treatment = merged_results$Treatment,
      Mean = merged_results$Mean,
      SE = merged_results$SE,
      F_value = round(anova_result[["F value"]][1], 3),
      P_value = round(treatment_p, 4),
      Significance = ifelse(treatment_p < 0.001, "***",
                           ifelse(treatment_p < 0.01, "**",
                                 ifelse(treatment_p < 0.05, "*", "ns"))),
      Statistical_Group = as.character(merged_results$.group),
      Shapiro_P = round(shapiro_test$p.value, 4),
      Shapiro_OK = ifelse(shapiro_test$p.value > 0.05, "Yes", "No"),
      Levene_P = round(levene_test[1,3], 4),
      Levene_OK = ifelse(levene_test[1,3] > 0.05, "Yes", "No"),
      Bartlett_P = round(bartlett_test$p.value, 4),
      Bartlett_OK = ifelse(bartlett_test$p.value > 0.05, "Yes", "No")
    )
    
    all_results <- rbind(all_results, variety_result)
  }
  
  return(all_results)
}

# ============================================================================
# RUN ANALYSIS
# ============================================================================
cat("\n=== RUNNING ANALYSIS ===\n")

varieties <- levels(final_metrics$Variety)

# Run analyses for each metric
dd_chilling_analysis <- analyze_metric(final_metrics, "Final_DD_Chilling", 
                                       "Ferguson_DD_Chilling", varieties)
dd_heating_analysis <- analyze_metric(final_metrics, "Final_DD_Heating", 
                                     "Ferguson_DD_Heating", varieties)
traditional_gdd_analysis <- analyze_metric(final_metrics, "Final_Traditional_GDD", 
                                          "Traditional_GDD", varieties)

# Combine results
combined_results <- rbind(dd_chilling_analysis, dd_heating_analysis, traditional_gdd_analysis)

# Save results
csv_filename <- paste0(results_dir, "/Combined_CU_GDD_Analysis_Results_", SEASON_TO_ANALYZE, ".csv")
write.csv(combined_results, csv_filename, row.names = FALSE)

cat("Results saved:", csv_filename, "\n")

# ============================================================================
# CREATE VISUALIZATION
# ============================================================================
cat("\n=== CREATING VISUALIZATION ===\n")

# Function to get statistical letters
get_statistical_letters <- function(metric_name) {
  stat_letters <- combined_results %>%
    dplyr::filter(Metric == metric_name) %>%
    dplyr::select(Variety, Treatment, Statistical_Group) %>%
    dplyr::rename(.group = Statistical_Group)
  
  stat_letters$Variety <- factor(stat_letters$Variety, levels = variety_order)
  stat_letters$Treatment <- factor(stat_letters$Treatment, levels = treatment_order)
  
  return(stat_letters)
}

# Create individual metric plots
create_metric_plot <- function(data, metric_column, metric_name, y_label, title_text) {
  stat_letters <- get_statistical_letters(metric_name)
  
  # Calculate max values for letter positioning
  max_values <- data %>%
    group_by(Variety, Treatment) %>%
    summarise(max_val = max(.data[[metric_column]], na.rm = TRUE), .groups = "drop")
  
  letters_positioned <- merge(stat_letters, max_values, by = c("Variety", "Treatment"))
  
  ggplot(data, aes(x = Variety, y = .data[[metric_column]], fill = Treatment)) +
    geom_boxplot(position = position_dodge(width = 0.8), 
                 alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
               aes(color = Treatment), size = 2, alpha = 0.9) +
    # geom_text(data = letters_positioned, 
    #           aes(x = Variety, y = max_val * 1.08, label = .group, group = Treatment),
    #           position = position_dodge(width = 0.8), 
    #           size = 8, family = "Arial", vjust = 0, color = "black") +
    scale_fill_manual(values = treatment_colors, name = "Treatment",
                     labels = c("Paper Buried", "Paper Raised", "No-Tube", "Plastic Buried", "Plastic Raised")) +
    scale_color_manual(values = treatment_colors, name = "Treatment",
                      labels = c("Paper Buried", "Paper Raised", "No-Tube", "Plastic Buried", "Plastic Raised")) +
    scale_x_discrete(labels = c("Cabernet\nSauvignon", "Chardonnay", "Concord", "Mourvedre")) +
    labs(title = title_text, 
         x = "", 
         y = y_label) +
    theme_minimal(base_size = 16, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 23, family = "Arial", hjust = 0.5, vjust = -1),
      axis.title = element_text(size = 22, family = "Arial"),
      axis.title.y = element_text(size = 22, family = "Arial", margin = margin(r = 15)),
      axis.text.x = element_text(size = 20, family = "Arial"),
      axis.text.y = element_text(size = 20, family = "Arial"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
}

# Create traditional GDD plot (collapsed across varieties)
create_traditional_gdd_plot <- function(data) {
  # Get statistical letters for Traditional GDD
  stat_letters_trad <- get_statistical_letters("Traditional_GDD")
  
  # Since we're collapsing across varieties, we need to get one set of letters per treatment
  # Take the first occurrence of each treatment (they should be the same across varieties for this analysis)
  stat_letters_collapsed <- stat_letters_trad %>%
    group_by(Treatment) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(All = "All")
  
  # Calculate max values for letter positioning
  data_collapsed <- data %>% mutate(All = "All")
  max_values_trad <- data_collapsed %>%
    group_by(All, Treatment) %>%
    summarise(max_val = max(Final_Traditional_GDD, na.rm = TRUE), .groups = "drop")
  
  letters_positioned_trad <- merge(stat_letters_collapsed, max_values_trad, by = c("All" = "All", "Treatment"))
  
  ggplot(data_collapsed, aes(x = All, y = Final_Traditional_GDD, fill = Treatment)) +
    geom_boxplot(position = position_dodge(width = 0.8), 
                 alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
               aes(color = Treatment), size = 2, alpha = 0.9) +
    # geom_text(data = letters_positioned_trad, 
    #           aes(x = All, y = max_val * 1.08, label = .group, group = Treatment),
    #           position = position_dodge(width = 0.8), 
    #           size = 8, family = "Arial", vjust = 0, color = "black") +
    scale_fill_manual(values = treatment_colors, name = "Treatment",
                     labels = c("Paper Buried", "Paper Raised", "No-Tube", "Plastic Buried", "Plastic Raised")) +
    scale_color_manual(values = treatment_colors, name = "Treatment",
                      labels = c("Paper Buried", "Paper Raised", "No-Tube", "Plastic Buried", "Plastic Raised")) +
    labs(title = "Traditional GDD", x = "", y = "GDD (°C·days)") +
    theme_minimal(base_size = 16, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 23, family = "Arial", hjust = 0.5, vjust = -1),
      axis.title = element_text(size = 22, family = "Arial"),
      axis.title.y = element_text(size = 22, family = "Arial", margin = margin(r = 15)),
      axis.text.x = element_text(size = 20, family = "Arial"),
      axis.text.y = element_text(size = 20, family = "Arial"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
}

# Create combined plot
create_combined_plot <- function() {
  # Create individual plots
  p_chilling <- create_metric_plot(final_metrics, "Final_DD_Chilling", 
                                  "Ferguson_DD_Chilling", "Chilling units (°C·days)", "Chilling units")
  p_heating <- create_metric_plot(final_metrics, "Final_DD_Heating", 
                                 "Ferguson_DD_Heating", "GDD (°C·days)", "GDD-like heating")
  p_traditional <- create_traditional_gdd_plot(final_metrics)
  
  # Create shared legend with proper spacing
  legend_plot <- ggplot(final_metrics, aes(x = Variety, y = Final_DD_Chilling, fill = Treatment)) +
    geom_boxplot() +
    scale_fill_manual(values = treatment_colors, name = "",
                     labels = c("Paper Buried", "Paper Raised", "No-Tube", "Plastic Buried", "Plastic Raised")) +
    theme_minimal(base_family = "Arial") +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 20, family = "Arial"),
          legend.key.size = unit(1.2, "cm"),
          legend.margin = margin(t = 15),
          legend.box.spacing = unit(0.8, "cm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  # Extract legend
  legend <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]
  
  # Combine plots
  combined_plot <- grid.arrange(
    p_chilling, p_heating, p_traditional,
    textGrob("Cultivar", gp = gpar(fontsize = 21, fontfamily = "Arial")),
    legend,
    ncol = 3, nrow = 3,
    heights = c(4, 0.1, 0.5),
    widths = c(1.65, 1.65, 0.8),  # First two panels wider, third panel narrower
    layout_matrix = rbind(c(1, 2, 3),
                         c(4, 4, 4),
                         c(5, 5, 5))
  )
  
  return(combined_plot)
}

# Create and save plot
combined_plot <- create_combined_plot()

plot_filename <- paste0(plots_dir, "/Chilling_GDD_Analysis_", SEASON_TO_ANALYZE, ".png")
png(filename = plot_filename, width = 18, height = 6, units = "in", res = 600, bg = "white")
grid.draw(combined_plot)
dev.off()

cat("Plot saved as:", plot_filename, "\n")
cat("ANALYSIS COMPLETE\n")