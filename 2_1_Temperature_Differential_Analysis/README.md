# Temperature Differential Analysis

## Overview

This section contains an analysis pipeline for studying grow tube temperature effects during dormancy periods. This analysis evaluates the thermal effects of grow tube treatments relative to both internal no-tube controls (ΔTn) and external ambient conditions (ΔTa). Environmental variability during dormancy is addressed using two grouping strategies: thermal phases and time blocks.

## Repository Structure

```
Temperature_Differential_Analysis/
├── 0_1_thermal_phase_fig.py      # Threshold determination
├── 0_2_time_block_fig.py         # Time block validation  
├── 0_3_data_preparation.py       # Integrated data prep
├── 1_1_lmm_delta_tn.r            # ΔTn statistical analysis
├── 1_2_lmm_delta_ta.r            # ΔTa statistical analysis
├── 2_lmm_back_transform.r        # Result back-transformation
├── 3_1_lmm_heatmap_fig.r         # Heatmap visualizations
├── 3_2_lmm_thermalphase_fig.r    # Thermal phase figures
├── 3_3_lmm_timeblock_fig.r       # Time block figures
└── README.md                     # This documentation
```

## Pipeline Structure

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

<img src="https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/2_1_Temperature_Differential_Analysis/images/img.png" alt="Thermal Conditions Grouping" width="70%">

*Grouping of thermal conditions during dormancy across two seasons. (A) Yearly comparisons of daily mean air temperature, (B) diurnal temperature range, and (C) solar radiation during the study period. (D–F) Distributions of thermal variables used to define thresholds for classifying cold vs. warm days and cloudy vs. sunny conditions. (G) Diurnal patterns of mean air temperature and solar radiation for defining four-time blocks.*

This data-driven approach ensures that categorical assignments (Cold/Warm at 3.3°C, Cloudy/Sunny at 100 W/m²) are grounded in the actual environmental conditions of the study period, rather than relying on literature-based assumptions.

### Phase 1: Statistical Analysis (R)

**Linear Mixed-Effects Modeling:**
- `1_1_lmm_delta_tn.r` - Analysis of ΔTn (grow tube - no-tube control)
- `1_2_lmm_delta_ta.r` - Analysis of ΔTa (grow tube - ambient air)

**Features:**
- Automatic transformation assessment (original, log, sqrt, cube root)
- Three modeling approaches: thermal phase, time block, and full factorial
- Comprehensive model diagnostics and selection

### Phase 2: Back-transformation and Interpretation

**Data Transformation:**
- `2_lmm_back_transform.r` - Converts transformed statistical results (Logarithmic) back to temperature scale (°C)

### Phase 3: Visualization and Results

**Figure Generation:**
- `3_1_lmm_heatmap_fig.r` - Full factorial heatmaps showing treatment effects across all conditions
- `3_2_lmm_thermalphase_fig.r` - Treatment comparisons by thermal phase
- `3_3_lmm_timeblock_fig.r` - Treatment comparisons by time block

## Requirements

### Python Environment
**Version:** Python >= 3.8

**Required Packages:**
```bash
pip install pandas>=1.3.0 numpy>=1.21.0 matplotlib>=3.4.0 seaborn
```

### R Environment
**Version:** R >= 4.0.0 

**Required Packages:**
```r
install.packages(c("lme4", "emmeans", "car", "multcomp", "dplyr", "ggplot2"))
```

### Input Data Files Required
- `cleaned_weather_2023_2024.csv` - Weather station data (2023-2024 season)
- `cleaned_weather_2024_2025.csv` - Weather station data (2024-2025 season)
- `2024_data_grow_tube.csv` - Grow tube temperature measurements (2024)
- `2025_data_grow_tube.csv` - Grow tube temperature measurements (2025)

## Key Outputs

### Data Files
- `grow_tube_categorized.csv` - Dataset for statistical analysis
- `Combined_LMM_Delta_Ta_Complete_Results.csv` - ΔTa statistical results
- `Combined_LMM_Delta_Tn_Complete_Results.csv` - ΔTn statistical results
- `Back_Transformed_Delta_Ta_Results.csv` - ΔTa results in °C
- `Back_Transformed_Delta_Tn_Results.csv` - ΔTn results in °C

### Figures
- Thermal phase and time block distribution plots
- Treatment effect heatmaps
- Statistical comparison figures by thermal phase and time block

## Statistical Analysis Details
- **Software:** R with lme4 package (Bates et al., 2015)
- **Model Selection:** AIC-based with residual diagnostics
- **Transformations:** log, sqrt, cube root
- **Post-hoc Analysis:** Estimated marginal means with Tukey adjustment (α = 0.05)
- **Random Effects:** Treatment replicates and measurement dates

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/LICENSE) file for details.

---

*This research was conducted to provide evidence-based recommendations for grow tube deployment in cool-climate viticulture, balancing protective benefits against potential cold-related risks during dormancy.*
