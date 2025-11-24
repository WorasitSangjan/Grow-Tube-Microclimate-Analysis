# Grow Tube Microclimate Analysis

![License](https://img.shields.io/badge/License-MIT-blue.svg)
![Open Source](https://img.shields.io/badge/Open%20Source-Yes-brightgreen.svg)
![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![R](https://img.shields.io/badge/R-4.0+-blue.svg)
[![Platform](https://img.shields.io/badge/Platform-HOBO%20Data%20Logger-red.svg)](https://raspberrypi.org)
![Agriculture](https://img.shields.io/badge/Agriculture-Precision%20Viticulture-green.svg)
![Research](https://img.shields.io/badge/Research-USDA--ARS-navy.svg)

> **Repository for the research paper: "Grow Tube Material Selection Influences Thermal Microenvironments and Cold Hardiness Risk in Dormant Grapevines"**

## Overview

This repository contains the complete analytical workflow for studying the thermal effects of different grow tube materials and installation configurations on young grapevine microclimates during dormancy periods. The research addresses a critical knowledge gap in viticulture: understanding how grow tubes modify thermal conditions during dormancy and the potential implications for cold hardiness dynamics.

## Experimental Setup

<img src="https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/images/img_1.png" alt="Experimental Setup" width="70%">

*Field experimental setup showing (A) different grow tube treatments and (B) sensor installation.*

## Research Objectives

1. **Quantify thermal modifications** - Assess how different grow tube configurations alter internal temperatures relative to ambient conditions and no-tube controls
2. **Model physiological risk** - Simulate cold hardiness dynamics to assess deacclimation risk across grape cultivars
3. **Evaluate dormancy implications** - Analyze impacts on chilling unit and growing degree day accumulation

## Repository Structure

```
Grow-Tube-Microclimate-Analysis/
├── 1_Data_Preparation/                     # Data preprocessing and quality control
├── 2_1_Temperature_Differential_Analysis/  # Thermal effect quantification
├── 2_2_Cold_Hardiness_Simulation/          # Physiological modeling
├── 2_3_Chilling_GDD_Analysis/              # Dormancy satisfaction metrics
├── images/                                 # Repository figures
├── LICENSE                                 # License
└── README.md                               # This documentation
```

## Analysis Pipeline

<img src="https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/images/img_2.png" alt="Analysis Workflow" width="80%">

*Data analysis workflow.*


### Section 1: Data Preparation

Preprocessing of HOBO temperature logger data and weather station integration:
- Missing value imputation using cross-replicate means
- 15-minute temporal alignment between sensors and weather data
- Quality validation and visualization

**Key Output:** `merged_data_YYYY.csv` - Clean, integrated temperature dataset

### Section 2.1: Temperature Differential Analysis

Statistical evaluation of thermal effects using two complementary metrics:
- **ΔT(n)** - Deviation from internal no-tube control
- **ΔT(a)** - Deviation from external ambient conditions
- 
Environmental variability addressed through data-driven grouping:
- **Thermal phases** - Cold/Warm × Cloudy/Sunny (based on 3.3°C and 100 W/m² thresholds)
- **Time blocks** - Four diurnal periods (Morning Warming, Peak Heat, Afternoon Cooling, Night)

**Methods:** AIC-based selection, Linear Mixed-Effects Models with transformation assessment 

**Key Outputs:** 
- Treatment effect estimates across environmental conditions
- Heatmaps and statistical comparison figures
- Back-transformed results in °C

### Section 2.2: Cold Hardiness Simulation

Dynamic physiological modeling of cultivar-specific cold hardiness:
- Implementation of Ferguson et al. (2014) cold hardiness model
- Simulation across four grape cultivars with varying cold sensitivity
- Phase-specific analysis (fall acclimation, winter maintenance, spring deacclimation)

**Methods:** 
- Python-based model implementation with variety-specific parameters
- Linear regression for acclimation/deacclimation rates
- LMM and ANOVA for treatment and cultivar comparisons

**Key Outputs:**
- Daily cold hardiness predictions by treatment-cultivar
- Statistical comparisons of dynamic hardiness rates
- Three-panel visualization with treatment or cultivar groupings

### Section 2.3: Chilling and Growing Degree Day Analysis

Integration of physiological and traditional thermal metrics:
- **Ferguson DD Chilling** - Variety-specific chilling accumulation
- **Ferguson DD Heating** - Ecodormancy heat accumulation  
- **Traditional GDD** - Standard base-10°C calculation

**Methods:** One-way ANOVA 

**Key Outputs:**
- Cumulative thermal metric comparisons
- Three-panel figure (Chilling | Heating | Traditional GDD)

## Experimental Design

**Study Site:** Washington Soil Health Initiative Research Vineyard, Prosser, WA (46°15'15.8" N, 119°43'42.6" W)

**Experimental Design:** Randomized complete block design (4 replicates)

**Study Periods:** 
- 2023-2024 dormancy (Oct 25 - Mar 30) - with vines
- 2024-2025 dormancy (Oct 25 - Mar 30) - without vines

**Treatments:**

| Material | Installation | Code |
|----------|-------------|------|
| White paper | Buried  | Paper Low |
| White paper | Raised | Paper High |
| Beige plastic (double-wall) | Buried  | Plastic Low |
| Beige plastic (double-wall) | Raised | Plastic High |
| None | - | No-Tube |

**Cultivars:** Cabernet Sauvignon, Chardonnay, Concord, Mourvèdre

**Measurements:** HOBO Pro v2 data loggers at 5-minute intervals, aligned to 15-minute bins with weather data

## Key Findings

### Temperature Modifications
- **Plastic Low** tubes produced the strongest warming (+1.5-2.0°C during sunny conditions)
- **Paper tubes** showed minimal thermal modification or slight cooling
- **Installation height** significantly affected thermal performance (buried > raised)
- **Environmental conditions** modulated treatment effects (greatest differences during sunny, warm periods)

### Physiological Implications
- Warming effects varied by cultivar sensitivity
- Treatment impacts on fall acclimation and spring deacclimation rates
- Modified chilling and heat accumulation patterns
- Trade-offs between thermal protection and deacclimation risk

## Requirements

### Python Environment
```bash
pip install pandas>=1.3.0 numpy>=1.21.0 matplotlib>=3.4.0 scipy>=1.7.0 seaborn
```

### R Environment
```r
install.packages(c("lme4", "emmeans", "car", "ggplot2", "dplyr", 
                   "gridExtra", "RColorBrewer", "multcomp", "agricolae"))
```

## Usage

Each section contains detailed README files with specific instructions. General workflow:

```bash
# 1. Data Preparation
cd 1_Data_Preparation
python 1_hobo_data_imputation.py
python 2_merge_data.py

# 2. Temperature Differential Analysis
cd ../2_1_Temperature_Differential_Analysis
python 0_1_thermal_phase_fig.py
python 0_2_time_block_fig.py
python 0_3_data_preparation.py
Rscript 1_1_lmm_delta_tn.r
Rscript 1_2_lmm_delta_ta.r

# 3. Cold Hardiness Simulation
cd ../2_2_Cold_Hardiness_Simulation
python 1_cold_hardiness_simulation.py
python 2_data_preparation_stat.py
Rscript 3_lmm_anova_cultivar.r

# 4. Chilling and GDD Analysis
cd ../2_3_Chilling_GDD_Analysis
Rscript 1_chilling_gdd_analysis.r
```

## Data Availability

Raw temperature data and weather station records are available on ... (Coming Soon!)

## Citation

If you use this code or methodology, please cite:

```bibtex
@article{sangjan2025growtube,
  title={Grow Tube Material Selection Influences Thermal Microenvironments and Cold Hardiness Risk in Dormant Grapevines},
  author={Sangjan, Worasit and Gillispie, Elizabeth C. and Schrader, Mark Jake and Shaw, Madison and Moyer, Michelle M.and Rippner, Devin A.},
  journal={In preparation},
  year={2025}
}
```

## Contact

- **Corresponding Author:** Devin A. Rippner - devin.rippner@usda.gov
- **Lead Author:** Worasit Sangjan - worasitsangjan.ws@gmail.com

## Acknowledgments

- USDA-ARS Horticultural Crops Production and Genetic Improvement Research Unit, Prosser, WA
- Washington State University Irrigated Agriculture Research and Extension Center
- Washington Soil Health Initiative Research Vineyard

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

*This research provides evidence-based recommendations for grow tube deployment in cool-climate viticulture, balancing protective benefits against potential cold-related risks during dormancy.*
