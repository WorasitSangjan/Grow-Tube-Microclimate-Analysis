# Author: Worasit Sangjan
# Date: 21 September 2025
# Version: 2.1

# ============================================================================
# FIGURE GENERATION: Temperature Effects by Time Block
# Generates separate figures for Delta_Ta and Delta_Tn
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)

cat("Creating time block effects figures for Delta_Ta and Delta_Tn...\n")

# ============================================================================
# CHECK AND LOAD DATA FILES
# ============================================================================

# Define file paths
delta_ta_file <- "Back_Transformed_Delta_Ta_Results.csv"
delta_tn_file <- "Back_Transformed_Delta_Tn_Results.csv"

# Check files exist
if (!file.exists(delta_ta_file) || !file.exists(delta_tn_file)) {
  stop("Error: Back-transformed files not found. Please ensure both files are in working directory.")
}

# Load data files
delta_ta_data <- read.csv(delta_ta_file, stringsAsFactors = FALSE)
delta_tn_data <- read.csv(delta_tn_file, stringsAsFactors = FALSE)

cat("Loaded Delta_Ta data:", nrow(delta_ta_data), "rows\n")
cat("Loaded Delta_Tn data:", nrow(delta_tn_data), "rows\n")

# ============================================================================
# DATA PREPARATION FUNCTION
# ============================================================================

prepare_timeblock_data <- function(data, variable_name, model_filter = "Model_2_Time_Block") {
  
  # Filter for Model 2 (time block only) and remove missing values
  plot_data <- data %>%
    filter(Model == model_filter) %>%
    filter(!is.na(EMM_Original_C)) %>%
    filter(!is.na(Statistical_Group)) %>%
    mutate(
      # Clean treatment names - sentence case
      treatment_clean = case_when(
        Treatment == "Paper_low" ~ "Paper\nLow",
        Treatment == "Paper_raise" ~ "Paper\nHigh", 
        Treatment == "Plastic_low" ~ "Plastic\nLow",
        Treatment == "Plastic_raise" ~ "Plastic\nHigh",
        TRUE ~ as.character(Treatment)
      ),
      
      # Clean time block names - sentence case
      time_block_clean = case_when(
        grepl("Morning", Time_Block, ignore.case = TRUE) ~ "Morning Warming",
        grepl("Peak", Time_Block, ignore.case = TRUE) ~ "Peak Heat",
        grepl("Evening", Time_Block, ignore.case = TRUE) ~ "Afternoon Cooling", 
        grepl("Night", Time_Block, ignore.case = TRUE) ~ "Night Period",
        TRUE ~ as.character(Time_Block)
      ),
      
      # Set variable name for identification
      Variable_Name = variable_name,
      
      # Use back-transformed values
      emmean = EMM,  # Log-transformed EMM for y-axis
      temperature_effect = EMM_Original_C,  # Celsius values for coloring
      se_value = SE,
      significance_letters = trimws(Statistical_Group)
    ) %>%
    
    # Set factor levels for proper ordering
    mutate(
      treatment_clean = factor(treatment_clean, 
                              levels = c("Paper\nLow", "Paper\nHigh", "Plastic\nLow", "Plastic\nHigh")),
      time_block_clean = factor(time_block_clean, 
                               levels = c("Morning Warming", "Peak Heat", "Afternoon Cooling", "Night Period"))
    ) %>%
    
    # Calculate error bar positions
    mutate(
      y_min = emmean - se_value,
      y_max = emmean + se_value,
      letter_y_pos = emmean + se_value + 0.1
    )
  
  return(plot_data)
}

# Prepare data for both variables
ta_plot_data <- prepare_timeblock_data(delta_ta_data, "Delta_Ta")
tn_plot_data <- prepare_timeblock_data(delta_tn_data, "Delta_Tn")

cat("Prepared Delta_Ta plot data:", nrow(ta_plot_data), "rows\n")
cat("Prepared Delta_Tn plot data:", nrow(tn_plot_data), "rows\n")

# ============================================================================
# FIGURE CREATION FUNCTION
# ============================================================================

create_timeblock_figure <- function(plot_data, title_suffix, y_label, filename) {
  
  # Calculate temperature range for consistent color scaling
  temp_range <- range(plot_data$temperature_effect, na.rm = TRUE)
  max_abs_temp <- max(abs(temp_range))
  
  p <- ggplot(plot_data, aes(x = treatment_clean, y = emmean, fill = temperature_effect)) +
    
    # Main bars with black outline
    geom_col(width = 0.7, alpha = 0.8, color = "black", linewidth = 0.3) +
    
    # Error bars
    geom_errorbar(aes(ymin = y_min, ymax = y_max), 
                  width = 0.2, linewidth = 0.8, color = "black") +
    
    # Significance letters
    geom_text(aes(x = treatment_clean, y = letter_y_pos, label = significance_letters),
              size = 6, color = "black", vjust = 2.5, family = "Arial") +
    
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", 
               alpha = 0.7, linewidth = 0.8) +
    
    # Panel layout by time block
    facet_wrap(~ time_block_clean, scales = "free_y", ncol = 4) +
    
    # Temperature-based color scheme (symmetric around zero)
    scale_fill_gradient2(
      low = "#2166AC",      # Blue for cooling
      mid = "white",        # White for neutral  
      high = "#B2182B",     # Red for warming
      midpoint = 0,
      name = "Temperature\ndifference (°C)",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
      limits = c(-max_abs_temp, max_abs_temp)  # Symmetric scale
    ) +
    
    # Labels
    labs(
      x = "Treatment",
      y = y_label
    ) +
    
    # Theme and styling with Arial font
    theme_bw(base_size = 16, base_family = "Arial") +
    theme(
      # No titles, subtitles, or caption for cleanest look
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.text.x = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 18),
      legend.position = "left",
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "lightblue", color = "black"),  # Light blue to distinguish from thermal phase
      strip.text = element_text(size = 18),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    )
  
  # Save figure
  ggsave(filename, p, width = 16, height = 3.3, dpi = 600, bg = "white")
  
  return(p)
}

# ============================================================================
# CREATE AND SAVE FIGURES
# ============================================================================

# Create Delta_Ta figure
cat("\nCreating Delta_Ta time block figure...\n")
p_ta <- create_timeblock_figure(
  ta_plot_data, 
  "ΔTa (Tube - Air)",
  expression(paste(Delta, "Ta (transformed scale)")),
  "Temperature_Effects_Delta_Ta_TimeBlock.png"
)
cat("Saved: Temperature_Effects_Delta_Ta_TimeBlock.png\n")

# Create Delta_Tn figure  
cat("Creating Delta_Tn time block figure...\n")
p_tn <- create_timeblock_figure(
  tn_plot_data,
  "ΔTn (Tube - No-tube)", 
  expression(paste(Delta, "Tn (transformed scale)")),
  "Temperature_Effects_Delta_Tn_TimeBlock.png"
)
cat("Saved: Temperature_Effects_Delta_Tn_TimeBlock.png\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== FIGURE SUMMARY ===\n")

# Function to summarize each dataset
summarize_timeblock_effects <- function(data, var_name) {
  cat("\n", var_name, "SUMMARY:\n")
  
  temp_range <- range(data$temperature_effect, na.rm = TRUE)
  cat("Temperature range:", round(temp_range[1], 3), "to", round(temp_range[2], 3), "°C\n")
  
  # Find extremes
  max_effect <- data[which.max(abs(data$temperature_effect)), ]
  cat("Largest absolute effect:", round(max_effect$temperature_effect, 3), "°C -", 
      gsub("\n", " ", max_effect$treatment_clean), "in", gsub("\n", " ", max_effect$time_block_clean), "\n")
  
  # Count effects by time block
  timeblock_summary <- data %>%
    group_by(time_block_clean) %>%
    summarise(
      n_treatments = n(),
      avg_effect = round(mean(temperature_effect), 3),
      max_effect = round(max(abs(temperature_effect)), 3),
      .groups = "drop"
    )
  
  cat("Effects by time block:\n")
  print(timeblock_summary)
}

# Generate summaries
summarize_timeblock_effects(ta_plot_data, "DELTA_TA")
summarize_timeblock_effects(tn_plot_data, "DELTA_TN")

# Overall comparison
cat("\n=== TREATMENT RANKING COMPARISON ===\n")
ta_ranking <- ta_plot_data %>%
  group_by(treatment_clean) %>%
  summarise(Mean_Delta_Ta = round(mean(temperature_effect), 3), .groups = "drop") %>%
  arrange(desc(Mean_Delta_Ta))

tn_ranking <- tn_plot_data %>%
  group_by(treatment_clean) %>%
  summarise(Mean_Delta_Tn = round(mean(temperature_effect), 3), .groups = "drop") %>%
  arrange(desc(Mean_Delta_Tn))

comparison <- merge(ta_ranking, tn_ranking, by = "treatment_clean")
print(comparison)

cat("\n=== TIME BLOCK PATTERNS ===\n")