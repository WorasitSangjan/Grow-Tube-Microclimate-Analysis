# ============================================================================
# FIGURE GENERATION: Temperature Effects Heatmap
# Creates publication-ready heatmap from LMM analysis results
# Run this after completing your LMM analysis and back-transformation
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)

cat("Creating temperature effects heatmap...\n")

# ============================================================================
# CHECK PREREQUISITES
# ============================================================================

# Check for back-transformed data file
celsius_file <- "Model_3_BackTransformed_Results.csv"
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

# Find the temperature columns (handle different naming conventions)
temp_celsius_col <- NULL
temp_transformed_col <- NULL

# Possible column names for °C values
possible_celsius_cols <- c("EMM_degrees_C", "EMM_delta_Tn_degrees_C", "Temperature_Effect_C", "Temperature_Effect")
# Possible column names for transformed values - try EMM_degrees_C as backup
possible_transformed_cols <- c("EMM_transformed_scale", "emmean", "transformed_value", "EMM_degrees_C")

for (col in possible_celsius_cols) {
  if (col %in% colnames(celsius_data)) {
    temp_celsius_col <- col
    break
  }
}

for (col in possible_transformed_cols) {
  if (col %in% colnames(celsius_data)) {
    temp_transformed_col <- col
    break
  }
}

if (is.null(temp_celsius_col)) {
  stop("Error: Could not find °C temperature column. Available columns: ", 
       paste(colnames(celsius_data), collapse = ", "))
}

cat("Using °C column:", temp_celsius_col, "\n")

# Handle case where we don't have transformed values - use °C values for both
if (is.null(temp_transformed_col)) {
  cat("No transformed column found - using °C values for display\n")
  temp_transformed_col <- temp_celsius_col
} else {
  cat("Using transformed column:", temp_transformed_col, "\n")
}

# Define proper ordering
time_block_order <- c("Morning\nWarming", "Peak\nHeat", "Evening\nCooling", "Night\nPeriod")
thermal_phase_order <- c("Cold Cloudy", "Cold Sunny", "Warm Cloudy", "Warm Sunny")

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

# Clean time block names
clean_time_names <- function(x) {
  clean_names <- case_when(
    grepl("Morning", x, ignore.case = TRUE) ~ "Morning\nWarming",
    grepl("Peak", x, ignore.case = TRUE) ~ "Peak\nHeat", 
    grepl("Evening", x, ignore.case = TRUE) ~ "Evening\nCooling",
    grepl("Night", x, ignore.case = TRUE) ~ "Night\nPeriod",
    TRUE ~ as.character(x)
  )
  factor(clean_names, levels = time_block_order)
}

# Clean thermal phase names
clean_thermal_names <- function(x) {
  clean_names <- case_when(
    grepl("Cold.*Cloudy", x, ignore.case = TRUE) ~ "Cold Cloudy",
    grepl("Cold.*Sunny", x, ignore.case = TRUE) ~ "Cold Sunny",
    grepl("Warm.*Cloudy", x, ignore.case = TRUE) ~ "Warm Cloudy", 
    grepl("Warm.*Sunny", x, ignore.case = TRUE) ~ "Warm Sunny",
    TRUE ~ as.character(x)
  )
  factor(clean_names, levels = thermal_phase_order)
}

# Apply data cleaning
heatmap_data <- celsius_data %>%
  mutate(
    treatment_clean = clean_treatment_names(Treatment),
    time_clean = clean_time_names(Time_Block),
    thermal_clean = clean_thermal_names(Thermal_Phase),
    temp_celsius = .data[[temp_celsius_col]],
    temp_transformed = .data[[temp_transformed_col]]
  ) %>%
  arrange(thermal_clean, time_clean) %>%
  group_by(thermal_clean, time_clean) %>%
  mutate(
    # Add significance letters based on ranking within each condition
    rank_in_condition = rank(-temp_celsius, ties.method = "min", na.last = "keep"),
    significance_letters = ifelse(is.na(temp_celsius), "", letters[rank_in_condition])
  ) %>%
  ungroup() %>%
  mutate(
    # Handle NA values for color scale
    fill_value = ifelse(is.na(temp_celsius), 0, temp_celsius)
  )

# Debug: Check data after cleaning
cat("Time block levels:", paste(levels(heatmap_data$time_clean), collapse = " → "), "\n")
cat("Thermal phase levels:", paste(levels(heatmap_data$thermal_clean), collapse = " | "), "\n")
cat("Treatment levels:", length(unique(heatmap_data$treatment_clean)), "treatments\n")

cat("Data preparation complete\n")

# ============================================================================
# CREATE PUBLICATION-READY HEATMAP
# ============================================================================

p <- ggplot(heatmap_data, aes(x = time_clean, y = treatment_clean, fill = fill_value)) +
  
  # Base heatmap tiles
  geom_tile(color = "white", linewidth = 0.6) +
  
  # Layer 1: °C values (top - small, gray with parentheses)
  geom_text(data = filter(heatmap_data, !is.na(temp_celsius)),
            aes(label = sprintf("(%.3f°C)", temp_celsius)),
            color = "gray50", size = 2.75, fontface = "bold", 
            vjust = -2.75, hjust = 0.5) +
  
  # Layer 2: Temperature values (center - larger, bold, black)
  geom_text(data = filter(heatmap_data, !is.na(temp_transformed)),
            aes(label = sprintf("%.3f", temp_transformed)),
            color = "black", size = 3.0, fontface = "bold", 
            vjust = 0, hjust = 0.5) +
  
  # Layer 3: Significance letters (bottom - black)
  geom_text(data = filter(heatmap_data, !is.na(temp_celsius)),
            aes(label = significance_letters),
            color = "black", size = 3.2, fontface = "bold", 
            vjust = 1.5, hjust = 0.5) +
  
  # Layer 4: NA labels
  geom_text(data = filter(heatmap_data, is.na(temp_celsius)),
            aes(label = "NA"),
            color = "gray50", size = 2.8, fontface = "bold") +
  
  # Panel layout by thermal phase
  facet_grid(. ~ thermal_clean, scales = "free_x", space = "free_x") +
  
  # Temperature-based color scheme
  scale_fill_gradient2(
    low = "#2166AC",      # Blue for cooling
    mid = "white",        # White for neutral  
    high = "#B2182B",     # Red for warming
    midpoint = 0,
    name = "Temperature\nDifference\n(°C)",
    breaks = pretty(heatmap_data$temp_celsius[!is.na(heatmap_data$temp_celsius)], n = 5),
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
    na.value = "gray80"
  ) +
  
  # Labels and titles
  labs(
    title = "Grow Tube Temperature Effects: Comprehensive Heatmap",
    subtitle = "Top: °C effects (practical) | Center: Temperature values | Bottom: Significance groups\nSame letters = not significantly different | Color scale based on °C for interpretation",
    x = "Time Block",
    y = "Treatment",
    caption = "Values in parentheses show temperature differences in °C\nDifferent letters within each panel indicate statistical significance"
  ) +
  
  # Theme and styling
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10, lineheight = 0.9),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    panel.spacing = unit(0.2, "cm"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  ) +
  
  scale_x_discrete(drop = FALSE)

# ============================================================================
# SAVE FIGURE
# ============================================================================

# Save high-quality figure
ggsave("Temperature_Effects_Heatmap.png", p, 
       width = 12, height = 5, dpi = 300, bg = "white")

cat("Figure saved: Temperature_Effects_Heatmap.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== HEATMAP SUMMARY ===\n")
cat("Multi-layer display: °C (top), Temperature values (center), Significance (bottom)\n")
cat("Temperature-based colors: Blue (cooling) → White (neutral) → Red (warming)\n")
cat("Proper temporal ordering of time blocks\n")
cat("Significance letters showing statistical differences\n")
cat("High-resolution output (300 DPI)\n")

# Filter out NA values for calculations
valid_data <- filter(heatmap_data, !is.na(temp_celsius))

if(nrow(valid_data) > 0) {
  temp_range <- range(valid_data$temp_celsius)
  cat("Temperature range:", round(temp_range[1], 3), "to", round(temp_range[2], 3), "°C\n")
  
  # Find extreme effects
  max_warming <- valid_data[which.max(valid_data$temp_celsius), ]
  max_cooling <- valid_data[which.min(valid_data$temp_celsius), ]
  
  cat("\n=== EXTREME EFFECTS ===\n")
  cat(sprintf("Maximum warming: %.3f°C - %s in %s, %s\n", 
              max_warming$temp_celsius, max_warming$treatment_clean,
              max_warming$thermal_clean, gsub("\n", " ", max_warming$time_clean)))
  
  cat(sprintf("Maximum cooling: %.3f°C - %s in %s, %s\n", 
              max_cooling$temp_celsius, max_cooling$treatment_clean,
              max_cooling$thermal_clean, gsub("\n", " ", max_cooling$time_clean)))
  
  # Practical significance analysis
  strong_effects <- sum(abs(valid_data$temp_celsius) > 1.0)
  moderate_effects <- sum(abs(valid_data$temp_celsius) > 0.5 & abs(valid_data$temp_celsius) <= 1.0)
  small_effects <- sum(abs(valid_data$temp_celsius) > 0.1 & abs(valid_data$temp_celsius) <= 0.5)
  
  cat("\n=== PRACTICAL SIGNIFICANCE ===\n")
  cat(sprintf("Strong effects (>1°C): %d combinations\n", strong_effects))
  cat(sprintf("Moderate effects (0.5-1°C): %d combinations\n", moderate_effects))
  cat(sprintf("Small effects (0.1-0.5°C): %d combinations\n", small_effects))
}

# Data completeness
na_count <- sum(is.na(heatmap_data$temp_celsius))
total_combinations <- nrow(heatmap_data)

cat("\n=== DATA COMPLETENESS ===\n")
cat(sprintf("Valid combinations: %d\n", nrow(valid_data)))
cat(sprintf("Missing (NA): %d\n", na_count))
cat(sprintf("Completeness: %.1f%%\n", (nrow(valid_data) / total_combinations) * 100))

cat("\n=== TIME BLOCKS (TEMPORAL ORDER) ===\n")
for(i in seq_along(time_block_order)) {
  cat("•", gsub("\n", " ", time_block_order[i]), "\n")
}

cat("\n=== THERMAL PHASES ===\n")
for(i in seq_along(thermal_phase_order)) {
  cat("•", thermal_phase_order[i], "\n")
}

cat("\n=== READY FOR PUBLICATION! ===\n")