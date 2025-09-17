# Temperature Differential Analysis Pipeline

**Author:** Worasit Sangjan  
**Date:** May 2025

## Overview

This repository contains a complete analysis pipeline for studying grow tube temperature effects during dormancy periods. The analysis combines weather station data with grow tube temperature measurements to determine the thermal benefits of different protection materials and configurations.

## Complete Pipeline Structure

### Phase 0: Threshold Determination and Data Preparation

**Step 0.1 - Thermal Phase Analysis:**
- `0_1_thermal_phase_fig.py` - Analyzes weather data to determine data-driven temperature and solar radiation thresholds
- **Output:** Thermal phase distribution plots and threshold values (e.g., 3.3°C temperature threshold, 100 W/m² solar threshold)

**Step 0.2 - Time Block Analysis:**
- `0_2_time_block_fig.py` - Analyzes diurnal patterns to establish optimal time block boundaries
- **Output:** Time block analysis figures and validation of 4-period structure (07-11, 11-15, 15-19, 19-07)

**Step 0.3 - Integrated Data Preparation:**
- `0_3_data_preparation.py` - Uses thresholds from Steps 0.1 and 0.2 to prepare final analysis dataset
- **Output:** `grow_tube_categorized.csv` with temperature deltas and categorical assignments

### Phase 1: Statistical Analysis (R)

**Linear Mixed-Effects Modeling:**
- `1_1_lmm_delta_tn.r` - Analysis of ΔTn (tube - no-tube control)
- `1_2_lmm_delta_ta.r` - Analysis of ΔTa (tube - ambient air)

**Features:**
- Automatic transformation assessment (original, log, sqrt, cube root)
- Three modeling approaches: thermal phase, time block, and full factorial
- Comprehensive model diagnostics and selection

### Phase 2: Back-transformation and Interpretation

**Data Transformation:**
- `2_lmm_back_transform.r` - Converts transformed statistical results back to temperature scale (°C)

### Phase 3: Visualization and Results

**Figure Generation:**
- `3_1_lmm_heatmap_fig.r` - Full factorial heatmaps showing treatment effects across all conditions
- `3_2_lmm_thermalphase_fig.r` - Treatment comparisons by thermal phase
- `3_3_lmm_timeblock_fig.r` - Treatment comparisons by time block

## Data Grouping Strategy

The pipeline uses data-driven thresholds determined in Phase 0:

**Thermal Phases** (from `0_1_thermal_phase_fig.py`):
- **Cold Cloudy**: Temperature < 3.3°C, Solar radiation < 100 W/m²
- **Cold Sunny**: Temperature < 3.3°C, Solar radiation ≥ 100 W/m²
- **Warm Cloudy**: Temperature ≥ 3.3°C, Solar radiation < 100 W/m²
- **Warm Sunny**: Temperature ≥ 3.3°C, Solar radiation ≥ 100 W/m²

**Time Blocks** (from `0_2_time_block_fig.py`):
- **Morning Warming** (07:00-11:00)
- **Peak Heat** (11:00-15:00)
- **Evening Cooling** (15:00-19:00)
- **Night Period** (19:00-07:00)

## Temperature Metrics

**ΔTa - Ambient Temperature Difference:**
```
ΔTa = T_tube - T_ambient
```

**ΔTn - No-tube Control Difference:**
```
ΔTn = T_tube - T_no-tube
```

## Usage Instructions

## Requirements

### Python Environment
**Version:** Python >= 3.8

**Required Packages:**
```bash
pip install pandas>=1.3.0 numpy>=1.21.0 matplotlib>=3.4.0 seaborn
```

**Core Libraries Used:**
- `pandas>=1.3.0` - Data manipulation and analysis
- `numpy>=1.21.0` - Numerical computations  
- `matplotlib>=3.4.0` - Plotting and visualization
- `seaborn` - Statistical data visualization
- `warnings` - Warning control (built-in)

### R Environment
**Version:** R 4.0.0 or higher

**Required Packages:**
```r
install.packages(c("lme4", "emmeans", "car", "multcomp", "dplyr", "ggplot2"))
```

**Package Details:**
- `lme4` - Linear mixed-effects models
- `emmeans` - Estimated marginal means and post-hoc comparisons  
- `car` - Companion to Applied Regression (Type III ANOVA)
- `multcomp` - Multiple comparisons procedures
- `dplyr` - Data manipulation
- `ggplot2` - Grammar of graphics plotting

### Input Data Files Required
- `cleaned_weather_2023_2024.csv` - Weather station data (2023-2024 season)
- `cleaned_weather_2024_2025.csv` - Weather station data (2024-2025 season)
- `2024_data_grow_tube.csv` - Grow tube temperature measurements (2024)
- `2025_data_grow_tube.csv` - Grow tube temperature measurements (2025)

### Step-by-Step Execution

1. **Determine Thresholds:**
   ```bash
   python 0_1_thermal_phase_fig.py    # Generate thermal phase thresholds
   python 0_2_time_block_fig.py       # Validate time block structure
   ```

2. **Prepare Analysis Dataset:**
   ```bash
   python 0_3_data_preparation.py     # Create grow_tube_categorized.csv
   ```

3. **Statistical Analysis:**
   ```r
   source("1_1_lmm_delta_tn.r")       # ΔTn analysis
   source("1_2_lmm_delta_ta.r")       # ΔTa analysis
   source("2_lmm_back_transform.r")   # Back-transform results
   ```

4. **Generate Figures:**
   ```r
   source("3_1_lmm_heatmap_fig.r")    # Comprehensive heatmaps
   source("3_2_lmm_thermalphase_fig.r") # Thermal phase figures
   source("3_3_lmm_timeblock_fig.r")  # Time block figures
   ```

## Key Outputs

### Data Files
- `grow_tube_categorized.csv` - Final analysis dataset
- `Combined_LMM_Delta_Ta_Complete_Results.csv` - ΔTa statistical results
- `Combined_LMM_Delta_Tn_Complete_Results.csv` - ΔTn statistical results
- `Back_Transformed_Delta_Ta_Results.csv` - ΔTa results in °C
- `Back_Transformed_Delta_Tn_Results.csv` - ΔTn results in °C

### Figures
- Thermal phase and time block distribution plots
- Treatment effect heatmaps
- Statistical comparison figures by thermal phase and time block

## Statistical Analysis Details

**Software:** R with lme4 package (Bates et al., 2015)
**Model Selection:** AIC-based with residual diagnostics
**Transformations:** Cube root when needed for normality
**Post-hoc Analysis:** Estimated marginal means with Tukey adjustment (α = 0.05)
**Random Effects:** Treatment replicates and measurement dates

## Key Findings

**Plastic Low Treatment:**
- Maximum warming: +4.02°C during Warm Sunny conditions
- Peak diurnal effect: +2.22°C during Peak Heat period
- Consistently positive temperature effects across all conditions

**Environmental Interactions:**
- Sunny conditions amplify treatment effects
- Peak Heat period shows maximum differentiation
- Temperature effects vary significantly by both thermal phase and time block

## Repository Structure

```
Temperature_Differential_Analysis/
├── 0_1_thermal_phase_fig.py      # Threshold determination
├── 0_2_time_block_fig.py         # Time block validation  
├── 0_3_data_preparation.py       # Integrated data prep
├── 1_1_lmm_delta_tn.r           # ΔTn statistical analysis
├── 1_2_lmm_delta_ta.r           # ΔTa statistical analysis
├── 2_lmm_back_transform.r       # Result back-transformation
├── 3_1_lmm_heatmap_fig.r        # Heatmap visualizations
├── 3_2_lmm_thermalphase_fig.r   # Thermal phase figures
├── 3_3_lmm_timeblock_fig.r      # Time block figures
└── README.md                     # This documentation
```

This pipeline ensures that all categorical assignments and thresholds are data-driven and empirically validated before statistical analysis.