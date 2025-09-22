# Author: Worasit Sangjan
# Date: 21 September 2025
# Version: 2.1

# ============================================================================
# FIGURE GENERATION: Full factorial
# Generates separate heatmaps for Delta_Ta and Delta_Tn
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)

cat("Creating temperature effects HEATMAPS for Delta_Ta and Delta_Tn...\n")

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

# File names for back-transformed results
delta_ta_file <- "Back_Transformed_Delta_Ta_Results.csv"
delta_tn_file <- "Back_Transformed_Delta_Tn_Results.csv"

# Check files exist
if (!file.exists(delta_ta_file) || !file.exists(delta_tn_file)) {
  stop("Error: Back-transformed files not found. Run back-transformation first.")
}

# Load Model 3 data only
delta_ta_data <- read.csv(delta_ta_file, stringsAsFactors = FALSE) %>%
  filter(Model == "Model_3_Full_Factorial")

delta_tn_data <- read.csv(delta_tn_file, stringsAsFactors = FALSE) %>%
  filter(Model == "Model_3_Full_Factorial")

cat("Delta_Ta Model 3 data:", nrow(delta_ta_data), "rows\n")
cat("Delta_Tn Model 3 data:", nrow(delta_tn_data), "rows\n")

# ============================================================================
# DATA CLEANING FUNCTION
# ============================================================================

clean_heatmap_data <- function(data, variable_name) {
  
  cleaned <- data %>%
    # Remove rows with missing temperature data
    filter(!is.na(EMM_Original_C)) %>%
    
    mutate(
      # Clean treatment names - sentence case
      Treatment_Clean = case_when(
        Treatment == "Paper_low" ~ "Paper\nLow",
        Treatment == "Paper_raise" ~ "Paper\nHigh", 
        Treatment == "Plastic_low" ~ "Plastic\nLow",
        Treatment == "Plastic_raise" ~ "Plastic\nHigh",
        TRUE ~ as.character(Treatment)
      ),
      
      # Clean time block names - sentence case
      Time_Block_Clean = case_when(
        grepl("Morning", Time_Block, ignore.case = TRUE) ~ "Morning\nWarming",
        grepl("Peak", Time_Block, ignore.case = TRUE) ~ "Peak\nHeat", 
        grepl("Evening", Time_Block, ignore.case = TRUE) ~ "Afternoon\nCooling",
        grepl("Night", Time_Block, ignore.case = TRUE) ~ "Night\nPeriod",
        TRUE ~ as.character(Time_Block)
      ),
      
      # Clean thermal phase names - sentence case
      Thermal_Phase_Clean = case_when(
        grepl("cold.*cloudy", Thermal_Phase, ignore.case = TRUE) ~ "Cold Cloudy",
        grepl("cold.*sunny", Thermal_Phase, ignore.case = TRUE) ~ "Cold Sunny",
        grepl("warm.*cloudy", Thermal_Phase, ignore.case = TRUE) ~ "Warm Cloudy", 
        grepl("warm.*sunny", Thermal_Phase, ignore.case = TRUE) ~ "Warm Sunny",
        TRUE ~ as.character(Thermal_Phase)
      ),
      
      # Use back-transformed temperature values
      Temperature_C = round(EMM_Original_C, 3),
      Log_EMM = round(EMM, 3),  # Original log-scale EMM for display
      
      Variable_Name = variable_name
    ) %>%
    
    # Set proper factor ordering with sentence case
    mutate(
      Treatment_Clean = factor(Treatment_Clean, levels = c("Paper\nLow", "Paper\nHigh", "Plastic\nLow", "Plastic\nHigh")),
      Time_Block_Clean = factor(Time_Block_Clean, levels = c("Morning\nWarming", "Peak\nHeat", "Afternoon\nCooling", "Night\nPeriod")),
      Thermal_Phase_Clean = factor(Thermal_Phase_Clean, levels = c("Cold Cloudy", "Cold Sunny", "Warm Cloudy", "Warm Sunny"))
    )
  
  return(cleaned)
}

# Clean both datasets
ta_clean <- clean_heatmap_data(delta_ta_data, "Delta_Ta")
tn_clean <- clean_heatmap_data(delta_tn_data, "Delta_Tn")

cat("Cleaned Delta_Ta:", nrow(ta_clean), "valid combinations\n")
cat("Cleaned Delta_Tn:", nrow(tn_clean), "valid combinations\n")

# ============================================================================
# HEATMAP CREATION FUNCTION
# ============================================================================

create_heatmap <- function(data, title_suffix) {
  
  # Create temperature range for consistent color scaling
  temp_range <- range(data$Temperature_C, na.rm = TRUE)
  
  p <- ggplot(data, aes(x = Time_Block_Clean, y = Treatment_Clean, fill = Temperature_C)) +
    
    # Base heatmap tiles
    geom_tile(color = "white", linewidth = 0.8) +
    
    # Layer 1: °C values (top - smaller, gray)
    geom_text(aes(label = sprintf("(%.3f°C)", Temperature_C)),
              color = "gray40", size = 3.1,
              vjust = -2.5, family = "Arial") +
    
    # Layer 2: Log EMM values (center - larger, black)
    geom_text(aes(label = sprintf("%.3f", Log_EMM)),
              color = "black", size = 3.4,
              vjust = 0.25, family = "Arial") +
    
    # Layer 3: Statistical groups (bottom)
    geom_text(aes(label = Statistical_Group),
              color = "black", size = 3.6,
              vjust = 1.75, family = "Arial") +
    
    # Panel layout by thermal phase
    facet_grid(. ~ Thermal_Phase_Clean, scales = "free_x", space = "free_x") +
    
    # Temperature-based color scheme with accessibility
    scale_fill_gradient2(
      low = "#2166AC",      # Blue for cooling
      mid = "white",        # White for neutral  
      high = "#B2182B",     # Red for warming
      midpoint = 0,
      name = "Temperature\ndifference (°C)",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
      limits = c(min(-abs(max(abs(temp_range)))), max(abs(temp_range)))  # Symmetric scale
    ) +
    
    # Labels and styling
    labs(
      x = "Time block",
      y = "Treatment"
    ) +
    
    # Clean theme with Arial font
    theme_minimal(base_size = 11, base_family = "Arial") +
    theme(
      # No titles for cleaner look
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 13),
      axis.title = element_text(size = 14),
      
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      
      strip.text = element_text(size = 13),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      
      legend.position = "left",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      
      panel.spacing = unit(0.3, "cm"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
  
  return(p)
}

# ============================================================================
# CREATE AND SAVE SEPARATE FIGURES
# ============================================================================

# Create Delta_Ta heatmap
cat("\nCreating Delta_Ta heatmap...\n")
p_ta <- create_heatmap(ta_clean, "ΔTa (Tube - Air)")
ggsave("Temperature_Effects_Delta_Ta_Model3.png", p_ta, 
       width = 12, height = 4.75, dpi = 600, bg = "white")
cat("Saved: Temperature_Effects_Delta_Ta_Model3.png\n")

# Create Delta_Tn heatmap
cat("Creating Delta_Tn heatmap...\n")
p_tn <- create_heatmap(tn_clean, "ΔTn (Tube - No-tube)")
ggsave("Temperature_Effects_Delta_Tn_Model3.png", p_tn, 
       width = 12, height = 4.75, dpi = 600, bg = "white")
cat("Saved: Temperature_Effects_Delta_Tn_Model3.png\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")

# Function to create summary for each dataset
create_summary <- function(data, var_name) {
  cat("\n", var_name, "SUMMARY:\n")
  
  temp_range <- range(data$Temperature_C, na.rm = TRUE)
  cat("Temperature range:", round(temp_range[1], 3), "to", round(temp_range[2], 3), "°C\n")
  
  # Find extremes
  max_warming <- data[which.max(data$Temperature_C), ]
  max_cooling <- data[which.min(data$Temperature_C), ]
  
  cat(sprintf("Max warming: %.3f°C - %s in %s, %s\n", 
              max_warming$Temperature_C, 
              gsub("\n", " ", max_warming$Treatment_Clean),
              max_warming$Thermal_Phase_Clean, 
              max_warming$Time_Block_Clean))
  
  cat(sprintf("Max cooling: %.3f°C - %s in %s, %s\n", 
              max_cooling$Temperature_C, 
              gsub("\n", " ", max_cooling$Treatment_Clean),
              max_cooling$Thermal_Phase_Clean, 
              max_cooling$Time_Block_Clean))
  
  # Effect distribution
  large_effects <- sum(abs(data$Temperature_C) > 2.0)
  medium_effects <- sum(abs(data$Temperature_C) > 1.0 & abs(data$Temperature_C) <= 2.0)
  small_effects <- sum(abs(data$Temperature_C) > 0.5 & abs(data$Temperature_C) <= 1.0)
  minimal_effects <- sum(abs(data$Temperature_C) <= 0.5)
  
  cat("Effect magnitudes:\n")
  cat("  Large (>2°C):", large_effects, "\n")
  cat("  Medium (1-2°C):", medium_effects, "\n") 
  cat("  Small (0.5-1°C):", small_effects, "\n")
  cat("  Minimal (<0.5°C):", minimal_effects, "\n")
}

# Generate summaries
create_summary(ta_clean, "DELTA_TA")
create_summary(tn_clean, "DELTA_TN")

# Treatment ranking comparison
cat("\nTREATMENT RANKING COMPARISON:\n")
ta_ranking <- ta_clean %>%
  group_by(Treatment_Clean) %>%
  summarise(Mean_Effect_Ta = round(mean(Temperature_C), 3), .groups = "drop") %>%
  arrange(desc(Mean_Effect_Ta))

tn_ranking <- tn_clean %>%
  group_by(Treatment_Clean) %>%
  summarise(Mean_Effect_Tn = round(mean(Temperature_C), 3), .groups = "drop") %>%
  arrange(desc(Mean_Effect_Tn))

comparison <- merge(ta_ranking, tn_ranking, by = "Treatment_Clean")
print(comparison)

cat("\n=== ANALYSIS COMPLETE ===\n")