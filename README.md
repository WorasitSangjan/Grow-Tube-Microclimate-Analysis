# Grow Tube Microclimate Analysis

Repository for the research paper: **"Grow Tubes During Dormancy: Modeling Thermal Protection vs. Risk in Vineyard Microclimates"**

## Key Features

- **Temperature Data Processing**: Tools for handling HOBO Pro v2 logger data and weather station integration
- **Temperature Differences Analysis**: 
- **Cold Hardiness Modeling**: 
- **Chilling Unit & Growing Degree Day Analysis**: 

## Treatments Analyzed

- **Paper Tubes**: White paper grow tubes (Low and Raised installation)
- **Plastic Tubes**: Beige double-walled plastic grow tubes (Low and Raised installation)  
- **Control**: No-tube control for baseline comparisons

## Installation
- **Low**: 
- **Raise**: 

## Key Results

- **Plastic Low** tubes consistently produced the strongest warming effects (+1.5-2.0°C during sunny conditions)
- **Paper tubes** generally showed minimal thermal modification or slight cooling effects
- **Installation height** significantly affected thermal performance, with buried tubes showing enhanced warming
- **Thermal phases** revealed treatment-specific responses to environmental conditions

## Data Sources

- **HOBO Pro v2 Data Loggers**: 5-minute interval temperature measurements
- **WSU AgWeatherNet**: Weather station data for 2023-2024
- **ATMOS 41 Weather Station**: On-site measurements for 2024-2025

## Statistical Methods

- Linear Mixed-Effects Models (LMMs) using R's lme4 package
- Temperature differential metrics: ΔT(a) and ΔT(n)
- Post-hoc comparisons with Tukey adjustment (α = 0.05)
- Cold hardiness simulation for 4 grapevine cultivars

## Citation

If you use this code or data, please cite:

```bibtex
Coming Soon!
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Corresponding Author**: Devin A. Rippner (devin.rippner@usda.gov)
- **Lead Author**: Worasit Sangjan

## Acknowledgments

- USDA-ARS Horticultural Crops Production and Genetic Improvement Research Unit
- Washington State University Irrigated Agriculture Research and Extension Center
- Washington Soil Health Initiative Research Vineyard

---

*This research was conducted to provide evidence-based recommendations for grow tube deployment in cool-climate viticulture, balancing protective benefits against potential cold-related risks during dormancy.*
