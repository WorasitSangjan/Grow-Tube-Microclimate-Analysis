# Data Preparation

This section contains Python scripts for preparing and preprocessing HOBO temperature data for grow tube microclimate analysis.

## Overview

The data preparation workflow consists of three main steps:
1. **Data Imputation** - Handle missing values in HOBO temperature data
2. **Data Merging** - Combine HOBO data with weather station data
3. **Data Visualization** - Create comprehensive plots for data validation

## Scripts

### 1. HOBO Data Imputation (`1_hobo_data_imputation.py`)

**Purpose**: Handles missing data imputation for grow tube temperature data using cross-replicate means.

**Key Features**:
- Loads HOBO data with MM/DD/YY H:MM timestamp format
- Creates a complete 5-minute grid with fuzzy timestamp matching
- Imputes missing values using replicate means within treatment groups
- Generates validation plots showing original vs imputed data

**Input**: Raw HOBO CSV file with columns: `Timestamp`, `Material`, `Condition`, `Replicate`, `Value`

**Output**: `hobo_data_imputed_complete.csv` - Complete dataset with imputed values

**Usage**:
```python
python 1_hobo_data_imputation.py
```

**Parameters**:
- `tolerance_minutes`: Tolerance for timestamp matching (default: 3 minutes)
- `create_plots`: Whether to generate visualization plots (default: True)

### 2. Data Merging (`2_merge_data.py`)

**Purpose**: Merges HOBO temperature data with weather station data using 15-minute time bins.

**Key Features**:
- Aligns HOBO and weather data to common 15-minute bins
- Handles timestamp misalignment with ±7.5-minute tolerance
- Aggregates HOBO data by treatment and replicate
- Selects the closest weather records for each time bin

**Input**: 
- HOBO data: `hobo_data_Year_complete.csv`
- Weather data: `weather_data_Year.csv`

**Output**: `merged_data_Year.csv` - Integrated dataset ready for analysis

**Usage**:
```python
python 2_merge_data.py
```

**Parameters**:
- `hobo_file`: Path to HOBO data file
- `weather_file`: Path to weather data file
- `output_file`: Output filename for merged data

### 3. Data Preview (`3_data_preview_fig.py`)

**Purpose**: Creates comprehensive visualizations of merged temperature data for validation.

**Key Features**:
- Filters data from specified start date (default: YYYY-MM-DD)
- Creates monthly detailed plots for each treatment combination
- Generates summary plots with all treatments
- Maps treatment codes to readable names

**Input**: `merged_data_Year.csv`

**Output**: Interactive matplotlib plots showing temperature patterns

**Usage**:
```python
python 3_data_preview_fig.py
```

**Parameters**:
- `input_file`: Path to merged data file
- `start_date`: Start date for filtering (default: "YYYY-MM-DD")
- `create_monthly`: Generate monthly plots (default: True)
- `create_summary`: Generate summary plots (default: True)

## Example Output

The data preparation pipeline produces high-quality temperature datasets ready for analysis:

![Temperature Data Example](https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/1_Data_Preparation/images/img.png)

*Temperature patterns for Plastic Low and Plastic Raise treatments showing grow tube effects from October 2024 to March 2025. Individual replicates (colored lines) track closely with some treatment-specific differences, while ambient air temperature (dashed black line) provides reference conditions.*

This visualization demonstrates successful data imputation, merging, and quality validation across the complete workflow.

## Data Flow

```
Raw HOBO Data
       ↓
1. Data Imputation
   - Handle missing values
   - Create a 5-minute grid
   - Cross-replicate imputation
       ↓
Complete HOBO Data
       ↓
2. Data Merging
   - Align with weather data
   - 15-minute binning
   - Timestamp matching
       ↓
Merged Dataset
       ↓
3. Data Visualization
   - Quality validation
   - Treatment comparisons
   - Temporal patterns
```

## Requirements

**Python Version:**
- Python >= 3.8

```
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
```

## Quality Checks

Each script includes validation steps:
- **Imputation**: Missing data statistics and temperature range validation
- **Merging**: Data completeness and timestamp alignment verification
- **Visualization**: Treatment combination counts and date range validation

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/WorasitSangjan/Grow-Tube-Microclimate-Analysis/blob/main/LICENSE) file for details.

---

*This research was conducted to provide evidence-based recommendations for grow tube deployment in cool-climate viticulture, balancing protective benefits against potential cold-related risks during dormancy.*
