# Thermal Differential Analysis

## Overview

This section contains the data analysis of thermal responses in different grow tube treatments under varying environmental conditions. The study examines how different tube materials and configurations affect internal temperature relative to ambient conditions and control treatments.

## Repository Structure

```
Grow-Tube-Microclimate-Analysis/
├── 2_1_Temperature_Differential_Analysis/
│   ├── 1_1_lmm_delta_ta.r
│   ├── 1_2_lmm_delta_tn.r
│   ├── 2_lmm_back_transform.r
│   ├── 3_1_lmm_thermalphase_fig.r
│   ├── 3_2_lmm_timeblock_fig.r
│   └── 3_3_lmm_heatmap_fig.r
```
## Files Description

### Temperature Differential Analysis Scripts (`2_1_Temperature_Differential_Analysis/`)
- `1_1_lmm_delta_ta.r`: Linear mixed model analysis for ΔT(a) - ambient temperature differences
- `1_2_lmm_delta_tn.r`: Linear mixed model analysis for ΔT(n) - no-tube control differences  
- `2_lmm_back_transform.r`: Back-transformation of cube root transformed data
- `3_1_lmm_thermalphase_fig.r`: Generates thermal phase comparison figures
- `3_2_lmm_timeblock_fig.r`: Creates time block analysis visualizations
- `3_3_lmm_heatmap_fig.r`: Produces comprehensive heatmap summaries




## Temperature Differential Analysis

This analysis evaluates the thermal effects of grow tube treatments relative to both external ambient conditions and internal no-tube controls. Environmental variability during dormancy is addressed using two grouping strategies: thermal phases and time blocks.

### Data Grouping

**Thermal Phases**: Classify days into four thermal phases based on daily mean air temperature and solar radiation thresholds:
- **Cold Cloudy**: Temperature < 3.5°C, Solar radiation < 100 W/m²
- **Cold Sunny**: Temperature < 3.5°C, Solar radiation ≥ 100 W/m²
- **Warm Cloudy**: Temperature ≥ 3.5°C, Solar radiation < 100 W/m²
- **Warm Sunny**: Temperature ≥ 3.5°C, Solar radiation ≥ 100 W/m²

**Time Blocks**: Divide each day into four time blocks based on solar radiation patterns:
- **Morning Warming** (07:00-11:00)
- **Peak Heat** (11:00-15:00)
- **Evening Cooling** (15:00-19:00)
- **Night Period** (19:00-07:00)

### Thermal Metrics

**ΔT(a) - Ambient Temperature Difference:**
```
ΔT(a) = T_tube - T_ambient
```
The difference between air temperature inside a grow tube and ambient air temperature at the same bin.

**ΔT(n) - No-tube Control Difference:**
```
ΔT(n) = T_tube - T_no-tube
```
The difference between air temperature inside a grow tube and the average no-tube control at the same bin.

### Statistical Analysis

Both metrics were analyzed using Linear Mixed-Effects Models (LMMs) in R via the lme4 package. The analysis included:

- **Software**: R (https://www.r-project.org/)
- **Package**: lme4 (Bates et al., 2015)
- **Model Selection**: Based on residual diagnostics and Akaike Information Criterion (AIC)
- **Transformations**: Cube root transformations applied when needed based on residual diagnostics
- **Post-hoc Analysis**: Estimated marginal means (emmeans package) with Tukey adjustment (α = 0.05)

### Modeling Frameworks

Three modeling frameworks were applied:
1. **Treatment × Thermal Phase**: Examining treatment effects across different weather conditions
2. **Treatment × Time Block**: Analyzing treatment effects across different times of day
3. **Full Factorial Interactions**: Complete interaction analysis

## Key Findings (Temperature Differential Analysis)

### Treatment Effects

**Plastic Low Treatment** showed the strongest warming effects:
- Maximum warming during **Warm Sunny** phase: +4.02°C above ambient
- Peak warming during **Peak Heat** period: +2.22°C above ambient
- Highest combined effect during **Warm Sunny/Peak Heat**: +4.06°C above ambient
- Consistently warmer across all conditions, including cooler periods

**Paper-based Treatments** remained consistently cooler:
- **Paper Raise** showed significantly lower temperatures only during warmest conditions
- Generally maintained temperatures closer to ambient conditions

### Environmental Response Patterns

- **Sunny conditions** (≥100 W/m²) produced greater temperature differentials than cloudy conditions
- **Peak Heat period** (11:00-15:00) showed maximum treatment effects
- **Plastic treatments** responded more dramatically to solar radiation
- **Time-of-day effects** were most pronounced during high solar radiation periods



## Usage

**Data Preparation**: See `1_Data_Preparation/README.md` for detailed instructions on data preprocessing.

**Temperature Differential Analysis**: 
1. Run `1_1_lmm_delta_ta.r` and `1_2_lmm_delta_tn.r` for statistical modeling
2. Execute `2_lmm_back_transform.r` if data transformations were applied
3. Generate figures using scripts `3_1_lmm_thermalphase_fig.r`, `3_2_lmm_timeblock_fig.r`, and `3_3_lmm_heatmap_fig.r`

## Requirements

- **R** (≥4.0.0) for statistical analysis and visualization
- Required R packages:
  - `lme4`: Linear mixed-effects models
  - `emmeans`: Estimated marginal means and post-hoc comparisons
  - Additional visualization packages as specified in individual scripts
  - `emmeans`: Estimated marginal means and post-hoc comparisons
  - Additional visualization packages as specified in individual scripts
