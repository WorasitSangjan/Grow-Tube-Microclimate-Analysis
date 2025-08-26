# Author: Worasit Sangjan
# Date: July 13, 2025

# ============================================================================
# FIGURE GENERATION: Temperature Effects by Time Block
# Creates publication-ready figure from LMM analysis results
# Run this after completing your LMM analysis (need emm_2 and back-transformed data)
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(emmeans)

cat("Creating time block effects figure...\n")

# ============================================================================
# CHECK PREREQUISITES
# ============================================================================

# Check for required objects from LMM analysis
if(!exists("emm_2")) {
  stop("Error: 'emm_2' object not found. Please run your LMM analysis first.")
}

# Check for back-transformed data file
celsius_file <- "Model_2_BackTransformed_Results.csv"
if(!file.exists(celsius_file)) {
  stop("Error: Back-transformed results file not found. Please run back-transformation first.")
}

cat("Prerequisites found\n")

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

# Load back-transformed results (°C values)
celsius_data <- read.csv(celsius_file)
cat("Loaded", nrow(celsius_data), "rows of °C data\n")

# Check and standardize column names
cat("Available columns:", paste(colnames(celsius_data), collapse = ", "), "\n")

# Find the temperature column (handle different naming conventions)
temp_col <- NULL
possible_temp_cols <- c("EMM_degrees_C", "EMM_delta_Tn_degrees_C", "EMM_delta_Ta_degrees_C", 
                       "emmean_original", "Temperature_Effect_C")

for (col in possible_temp_cols) {
  if (col %in% colnames(celsius_data)) {
    temp_col <- col
    break
  }
}

if (is.null(temp_col)) {
  stop("Error: Could not find temperature column in back-transformed data. Available columns: ", 
       paste(colnames(celsius_data), collapse = ", "))
}

cat("Using temperature column:", temp_col, "\n")

# Standardize column name
celsius_data$EMM_degrees_C <- celsius_data[[temp_col]]

# Get statistical results from emmeans
plot_data <- summary(emm_2)

# Clean treatment names
clean_treatment_names <- function(x) {
  case_when(
    x == "Paper_low" ~ "Paper\nLow",
    x == "Paper_raise" ~ "Paper\nRaise", 
    x == "Plastic_low" ~ "Plastic\nLow",
    x == "Plastic_raise" ~ "Plastic\nRaise",
    TRUE ~ as.character(x)
  )
}

# Clean time block names and set proper temporal order
clean_time_names <- function(x) {
  clean_names <- case_when(
    x == "Night Period" ~ "Night Period",
    x == "Morning Warming" ~ "Morning Warming",
    x == "Peak Heat" ~ "Peak Heat", 
    x == "Evening Cooling" ~ "Evening Cooling",
    TRUE ~ as.character(x)
  )
  factor(clean_names, levels = c("Morning Warming", "Peak Heat", "Evening Cooling", "Night Period"))
}

# Apply cleaning to plot data
plot_data$treatment_clean <- clean_treatment_names(plot_data$treatment)
plot_data$time_clean <- clean_time_names(plot_data$time_block)

# Apply cleaning to celsius data
celsius_data$treatment_clean <- clean_treatment_names(celsius_data$Treatment)
celsius_data$time_clean <- clean_time_names(celsius_data$Time_Block)

# Debug: Check data before merge
cat("Plot data treatments:", unique(plot_data$treatment_clean), "\n")
cat("Celsius data treatments:", unique(celsius_data$treatment_clean), "\n")
cat("Plot data time blocks:", levels(plot_data$time_clean), "\n")
cat("Celsius data time blocks:", levels(celsius_data$time_clean), "\n")

# Merge datasets
plot_data_final <- merge(
  plot_data, 
  celsius_data[, c("treatment_clean", "time_clean", "EMM_degrees_C")],
  by = c("treatment_clean", "time_clean"),
  all.x = TRUE
)

# Check if merge was successful
if (any(is.na(plot_data_final$EMM_degrees_C))) {
  cat("Warning: Some temperature values are missing after merge\n")
  print(plot_data_final[is.na(plot_data_final$EMM_degrees_C), c("treatment_clean", "time_clean")])
}

cat("Final dataset has", nrow(plot_data_final), "rows\n")

cat("Data preparation complete\n")

# ============================================================================
# GET SIGNIFICANCE LETTERS
# ============================================================================

# Generate compact letter display
cld_data <- summary(cld(emm_2, Letters = letters, decreasing = TRUE))

# Clean CLD data
cld_data$treatment_clean <- clean_treatment_names(cld_data$treatment)
cld_data$time_clean <- clean_time_names(cld_data$time_block)

# Clean significance letters and set position
cld_data$significance_letters <- gsub(" ", "", cld_data$.group)
cld_data$y_pos <- cld_data$emmean + cld_data$SE + 0.1

# Merge CLD data with temperature data for consistency
cld_data <- merge(
  cld_data,
  celsius_data[, c("treatment_clean", "time_clean", "EMM_degrees_C")],
  by = c("treatment_clean", "time_clean"),
  all.x = TRUE
)

cat("Significance letters prepared with", nrow(cld_data), "rows\n")

# ============================================================================
# CREATE PUBLICATION-READY FIGURE
# ============================================================================

p <- ggplot(plot_data_final, aes(x = treatment_clean, y = emmean, fill = EMM_degrees_C)) +
  
  # Main bars with error bars
  geom_col(width = 0.7, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, linewidth = 0.8, color = "black") +
  
  # Significance letters
  geom_text(data = cld_data, 
            aes(x = treatment_clean, y = y_pos, label = significance_letters),
            size = 6, fontface = "bold", color = "black") +
  
  # Reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", 
             alpha = 0.7, linewidth = 0.8) +
  
  # Panel layout
  facet_wrap(~ time_clean, scales = "free_y", ncol = 4) +
  
  # Temperature-based color scheme
  scale_fill_gradient2(
    low = "#2166AC",      # Blue for cooling
    mid = "white",        # White for neutral  
    high = "#B2182B",     # Red for warming
    midpoint = 0,
    name = "Temperature\nDifference (°C)",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  
  # Labels and titles
  labs(
    title = "Grow Tube Temperature Effects by Time Block",
    subtitle = "Letters indicate significant differences (Tukey HSD, Alpha = 0.05) | Colors show temperature impact",
    x = "Treatment",
    y = expression(paste(Delta, "T (transformed scale)")),
    caption = "Different letters within each panel indicate statistical significance\nBlue = cooling effect, Red = warming effect"
  ) +
  
  # Theme and styling
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 12, hjust = 0.5, face = "italic"),
    axis.text.x = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 13),
    strip.background = element_rect(fill = "lightblue", color = "black"),
    strip.text = element_text(face = "bold", size = 18),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

# Save high-quality figure
ggsave("Time_Block_Effects_Figure.png", p, 
       width = 16, height = 4.5, dpi = 300, bg = "white")

cat("Figure saved: Time_Block_Effects_Figure.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== FIGURE SUMMARY ===\n")
cat("Temperature-based colors: Blue (cooling) → White (neutral) → Red (warming)\n")
cat("Significance letters showing statistical differences\n") 
cat("Four time block panels in single row\n")
cat("Error bars showing ± 1 SE\n")
cat("High-resolution output (300 DPI)\n")

# Display temperature range
temp_range <- range(plot_data_final$EMM_degrees_C, na.rm = TRUE)
cat("Temperature range:", round(temp_range[1], 2), "to", round(temp_range[2], 2), "°C\n")

cat("\n=== TIME BLOCKS (TEMPORAL ORDER) ===\n")
cat("- Morning Warming\n")
cat("- Peak Heat\n")
cat("- Evening Cooling\n")
cat("- Night Period\n")

cat("\n=== READY FOR PUBLICATION! ===\n")
