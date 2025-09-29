# Section 2.2: Cold Hardiness Simulation

## Overview

This section analyzes cold hardiness dynamics in grape cultivars under different grow tube treatments across two dormancy seasons (2023-2024 and 2024-2025). The analysis implements a dynamic physiological model to predict cultivar-specific cold hardiness responses throughout dormancy, followed by statistical comparisons to evaluate both treatment effects and genetic variation.

## Repository Structure

```
2_2_Cold_Hardiness_Simulation/
├── 1_cold_hardiness_simulation.py    # Model implementation and prediction
├── 2_data_preparation_stat.py        # Phase-specific metric extraction
├── 3_lmm_anova_cultivar.r            # Treatment comparison (by cultivar)
├── 3_lmm_anova_treatment.r           # Cultivar comparison (by treatment)
└── README.md                         # This documentation
```

## Pipeline Structure

### Phase 1: Cold Hardiness Prediction

**Step 1 - Model Implementation:**
- `1_cold_hardiness_simulation.py` - Processes experimental data and generates cold hardiness predictions using a dynamic physiological model
- **Inputs:** 
  - `weather_YYYY_YYYY_binned.csv` - Ambient weather data (15-min intervals, Sep 1 - Oct 24)
  - `grow_tube_categorized.csv` - Treatment microclimate data (Oct 25 - Mar 31)
- **Outputs:**
  - `combined_detailed_results_[Variety].csv` - Daily predictions with model state variables
  - `combined_summary_[Variety].csv` - Summary statistics per treatment-replicate
  - Individual variety comparison plots (predicted hardiness vs. ambient temperature)
  - Treatment-based comparison plots (all varieties within each treatment)

<img src="https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/2_2_Cold_Hardiness_Simulation/images/example_plot.png" alt="Cold Hardiness Dynamics" width="70%">

*Example cold hardiness dynamics showing predicted hardiness trajectories for different treatments alongside ambient temperature extremes across the dormancy period.*

**Model Framework:**

The model simulates cold hardiness through two distinct dormancy phases:

1. **Endodormancy (Sep - Dec):** Internal growth cessation driven by photoperiod and temperature
2. **Ecodormancy (Jan - Mar):** External growth suppression controlled by temperature

Key model components:
- **Acclimation:** Temperature-driven increase in cold hardiness during cold periods
- **Deacclimation:** Temperature-driven decrease in cold hardiness during warm periods
- **Phase transition:** Base-10 chilling accumulation determines endodormancy → ecodormancy shift

**Variety-Specific Parameters:**
- Initial hardiness (Hc_initial)
- Minimum/maximum hardiness thresholds (Hc_min, Hc_max)
- Temperature thresholds for each dormancy phase
- Acclimation and deacclimation rates
- Ecodormancy boundary (base-10 chilling requirement)
- Theta parameter (controls deacclimation curve shape in ecodormancy)

### Phase 2: Statistical Data Preparation

**Step 2 - Phase-Specific Metric Extraction:**
- `2_data_preparation_stat.py` - Extracts metrics from simulation results for statistical analysis
- **Output:** Six CSV files containing slope and mean hardiness data organized by analysis type

**Seasonal Phases:**
- **Fall acclimation:** Julian days 250-330 (Sep 7 - Nov 26)
- **Winter maintenance:** Days 330-365 + 1-35/50 (Nov 27 - Feb 4/19)
- **Spring deacclimation:** Days 35/50-90 (Feb 5/20 - Mar 31)

**Calculated Metrics:**
- **Dynamic phases (fall/spring):** Linear regression slope (°C/day), R², sample size
- **Winter phase:** Mean hardiness, standard deviation, min/max values

### Phase 3: Statistical Analysis and Visualization

**Linear Mixed-Effects Modeling:**

**Step 3a - By Cultivar Analysis:**
- `3_lmm_anova_cultivar.r` - Compare treatments within each cultivar to assess grow tube effects
- **Statistical Models:**
  - Dynamic phases: `Slope ~ Treatment * Phase + (1|Replicate)`
  - Winter phase: `Mean_Hardiness ~ Treatment`
- **Outputs:**
  - `results_YYYY-YYYY_cultivar/LMM_Diagnostics_[Variety].txt` - Dynamic phase diagnostics
  - `results_YYYY-YYYY_cultivar/Winter_Diagnostics_[Variety].txt` - Winter phase diagnostics
  - `plots_YYYY-YYYY_cultivar/Cold_Hardiness_YYYY-YYYY.png` - Three-panel visualization with statistical groupings

**Step 3b - By Treatment Analysis:**
- `3_lmm_anova_treatment.r` - Compare cultivars within each treatment to assess genetic variation
- **Statistical Models:**
  - Dynamic phases: `Slope ~ Variety * Phase + (1|Replicate)`
  - Winter phase: `Mean_Hardiness ~ Variety`
- **Outputs:**
  - `results_YYYY-YYYY_by_treatment/LMM_Diagnostics_[Treatment].txt` - Dynamic phase diagnostics
  - `results_YYYY-YYYY_by_treatment/Winter_Diagnostics_[Treatment].txt` - Winter phase diagnostics
  - `plots_YYYY-YYYY_by_treatment/Cold_Hardiness_YYYY-YYYY.png` - Three-panel visualization with statistical groupings

**Enhanced Statistical Features:**
- Automatic transformation assessment (original, log, sqrt, cube root)
- AIC-based model selection with residual diagnostics
- Type III ANOVA for main effects and interactions
- Post-hoc comparisons with Tukey HSD adjustment
- Dual statistical grouping validation (emmeans CLD + agricolae HSD.test)
- Comprehensive variance homogeneity testing (Levene's and Bartlett's tests)
- Descriptive statistics with standard errors
- Statistical letters displayed on visualization plots

## Requirements

### Python Environment
**Version:** Python >= 3.8

**Required Packages:**
```bash
pip install pandas>=1.3.0 numpy>=1.21.0 matplotlib>=3.4.0 scipy>=1.7.0
```

### R Environment
**Version:** R >= 4.0.0

**Required Packages:**
```r
install.packages(c("lme4", "emmeans", "car", "ggplot2", "dplyr", 
                   "gridExtra", "RColorBrewer", "multcomp", "agricolae"))
```

### Input Data Files Required
- `weather_2023_2024_binned.csv` - Ambient weather data (2023-2024 season)
- `weather_2024_2025_binned.csv` - Ambient weather data (2024-2025 season)
- `grow_tube_categorized.csv` - Treatment microclimate data

## Usage

### Complete Analysis Pipeline
```bash
# 1. Run cold hardiness simulation
python 1_cold_hardiness_simulation.py

# 2. Prepare statistical datasets
python 2_data_preparation_stat.py

# 3. Statistical analysis (choose one or both)
Rscript 3_lmm_anova_cultivar.r
Rscript 3_lmm_anova_treatment.r
```

### Season-Specific Configuration

**Python (`1_cold_hardiness_simulation.py`):**
```python
SEASON = "2023_2024"  # Change to "2024_2025" for second season
VARIETIES = ['Chardonnay', 'Concord', 'Cabernet Sauvignon', 'Mourvedre']
```

**R (both analysis scripts):**
```r
SEASON_TO_ANALYZE <- "2023-2024"  # Change to "2024-2025"
```

## Treatments
- **Paper Low:** Paper tube at soil level
- **Paper High:** Paper tube raised 10cm above soil
- **Plastic Low:** Plastic tube at soil level
- **Plastic High:** Plastic tube raised 10cm above soil
- **No-Tube:** Uncovered control

## Key Outputs

### Data Files
- `dynamic_slopes_database.csv` - All slope calculations across seasons
- `winter_means_database.csv` - All winter hardiness means
- `dynamic_slopes_by_variety.csv` / `dynamic_slopes_by_treatment.csv` - Analysis-specific datasets
- `winter_means_by_variety.csv` / `winter_means_by_treatment.csv` - Analysis-specific datasets
- `combined_detailed_results_[Variety].csv` - Daily model predictions
- `combined_summary_[Variety].csv` - Treatment-replicate summaries

### Diagnostic Files
**Dynamic Phase Analysis:**
- `LMM_Diagnostics_[Variety/Treatment].txt` - Transformation comparisons, model summaries, ANOVA tables, assumption checks

**Winter Phase Analysis:**
- `Winter_Diagnostics_[Variety/Treatment].txt` - Comprehensive ANOVA diagnostics including:
  - Descriptive statistics (n, mean, SD, SE)
  - Normality tests (Shapiro-Wilk)
  - Homogeneity tests (Levene's and Bartlett's)
  - Residual analysis
  - Tukey HSD p-values
  - Multiple grouping methods (emmeans CLD and agricolae HSD)
  - HSD statistics (alpha, HSD value, F-statistic, p-value)
  - Raw data verification

### Figures
- Individual variety comparison plots showing predicted hardiness vs. ambient temperature
- Treatment-based comparison plots showing all varieties within treatments
- Three-panel statistical comparison figures (fall/winter/spring) with:
  - Boxplots with individual data points
  - Statistical grouping letters overlaid
  - Treatment or variety color coding

## Statistical Analysis Details
- **Software:** R with lme4 package (v2.2, September 2025)
- **Model Selection:** AIC-based comparison with residual normality assessment
- **Transformations:** Original, log, sqrt, cube root (original used for biological interpretation)
- **Post-hoc Analysis:** 
  - Estimated marginal means with Tukey HSD adjustment (α = 0.05)
  - Dual validation using emmeans and agricolae packages
- **Random Effects:** Replicate as random intercept
- **Assumption Checks:** 
  - Shapiro-Wilk tests for normality of residuals and random effects
  - Levene's test and Bartlett's test for homogeneity of variance
  - Residual diagnostic plots and statistics

## Key Assumptions
1. Model assumes homogeneous microclimate within treatment zones
2. Linear relationships adequate for acclimation/deacclimation processes
3. Replicate structure captures experimental variability
4. Original (untransformed) data used for biological interpretation even when transformations improve statistical model fit
5. Phase boundaries adequately capture distinct physiological periods (may differ slightly between seasons)
6. Statistical grouping letters provide visual summary of post-hoc comparisons

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/LICENSE) file for details.

---

*This research provides evidence-based insights into how grow tube treatments affect cold hardiness dynamics across genetically diverse grape cultivars, informing deployment strategies that balance protection with dormancy-related risks.*