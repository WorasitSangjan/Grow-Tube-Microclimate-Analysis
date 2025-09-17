# Author: Worasit Sangjan
# Date: May 9, 2025

# ============================================================================
# INTEGRATED and CATEGORIZED DATA 
# Combines grow tube data with weather data, calculates deltas, and categorizes
# Final output: grow_tube_categorized.csv
# ============================================================================

import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def load_and_clean_data() -> tuple:
    """Load grow tube and weather data, clean and standardize."""
    print("Loading and cleaning datasets...")
    
    # Load grow tube data
    grow_2024 = pd.read_csv("2024_data_grow_tube.csv")
    grow_2025 = pd.read_csv("2025_data_grow_tube.csv")
    
    # Load weather data
    weather_2023_24 = pd.read_csv("cleaned_weather_2023_2024.csv")
    weather_2024_25 = pd.read_csv("cleaned_weather_2024_2025.csv")
    
    print(f"Grow tube 2024: {len(grow_2024):,} records")
    print(f"Grow tube 2025: {len(grow_2025):,} records")
    print(f"Weather 2023-24: {len(weather_2023_24):,} records")
    print(f"Weather 2024-25: {len(weather_2024_25):,} records")
    
    # Clean grow tube data - remove old weather columns
    weather_columns_to_drop = ['air_temp_C', 'precip_mm', 'solar_radiation_Wpm2', 'wind_speed_mps']
    grow_2024_clean = grow_2024.drop(columns=weather_columns_to_drop, errors='ignore')
    grow_2025_clean = grow_2025.drop(columns=weather_columns_to_drop, errors='ignore')
    
    # Standardize timestamps
    grow_2024_clean['timestamp'] = pd.to_datetime(grow_2024_clean['timestamp'])
    grow_2025_clean['timestamp'] = pd.to_datetime(grow_2025_clean['timestamp'])
    weather_2023_24['timestamp'] = pd.to_datetime(weather_2023_24['timestamp'])
    weather_2024_25['timestamp'] = pd.to_datetime(weather_2024_25['timestamp'])
    
    # Combine weather data
    weather_combined = pd.concat([weather_2023_24, weather_2024_25], ignore_index=True)
    
    print(f"Combined weather: {len(weather_combined):,} records")
    
    return grow_2024_clean, grow_2025_clean, weather_combined

def merge_with_weather(grow_data: pd.DataFrame, weather_data: pd.DataFrame, year: str) -> pd.DataFrame:
    """Merge grow tube data with weather data."""
    print(f"Merging {year} data with weather...")
    
    merged = pd.merge(
        grow_data,
        weather_data[['timestamp', 'air_temp_C', 'precip_mm', 'solar_radiation_Wpm2', 'wind_speed_mps']],
        on='timestamp',
        how='left'
    )
    
    merge_success = merged['air_temp_C'].notna().sum() / len(merged) * 100
    print(f"{year} merge success: {merge_success:.1f}%")
    
    return merged

def filter_dormancy_period(data: pd.DataFrame) -> pd.DataFrame:
    """Filter to dormancy period: Oct 25 - Mar 30."""
    print("Filtering to dormancy period...")
    
    data['month'] = data['timestamp'].dt.month
    data['day'] = data['timestamp'].dt.day
    data['hour'] = data['timestamp'].dt.hour
    
    def is_dormancy_period(row):
        month, day = row['month'], row['day']
        return ((month == 10 and day >= 25) or 
                month in [11, 12, 1, 2] or 
                (month == 3 and day <= 30))
    
    dormancy_data = data[data.apply(is_dormancy_period, axis=1)].copy()
    print(f"Dormancy period data: {len(dormancy_data):,} records")
    
    return dormancy_data

def calculate_temperature_deltas(data: pd.DataFrame) -> pd.DataFrame:
    """Calculate delta_Ta and delta_Tn temperature differences."""
    print("Calculating temperature deltas...")
    
    # Filter valid data
    valid_data = data.dropna(subset=['Value', 'air_temp_C']).copy()
    
    # Calculate delta_Ta (grow tube - ambient)
    valid_data['delta_Ta'] = valid_data['Value'] - valid_data['air_temp_C']
    
    # Calculate control averages for delta_Tn
    control_data = valid_data[valid_data['Condition'] == 'uncover']
    control_avg = control_data.groupby('timestamp')['Value'].mean().reset_index()
    control_avg.columns = ['timestamp', 'control_avg_temp']
    
    # Merge treatment data with control averages
    treatment_data = valid_data[valid_data['Condition'] != 'uncover'].copy()
    merged_treatment = pd.merge(treatment_data, control_avg, on='timestamp', how='left')
    merged_treatment['delta_Tn'] = merged_treatment['Value'] - merged_treatment['control_avg_temp']
    
    # Add control data back (no delta_Tn for controls)
    control_with_delta = control_data.copy()
    control_with_delta['delta_Tn'] = np.nan
    control_with_delta['control_avg_temp'] = control_with_delta['Value']
    
    # Combine all data
    final_data = pd.concat([merged_treatment, control_with_delta], ignore_index=True)
    
    print(f"Delta_Ta calculated: {final_data['delta_Ta'].notna().sum():,} records")
    print(f"Delta_Tn calculated: {final_data['delta_Tn'].notna().sum():,} records")
    
    return final_data

def assign_time_blocks(data: pd.DataFrame) -> pd.DataFrame:
    """Assign time blocks based on hour of day."""
    print("Assigning time blocks...")
    
    def get_time_block(hour):
        if pd.isna(hour):
            return "Unknown"
        elif 7 <= hour < 11:
            return "Morning Warming"
        elif 11 <= hour < 15:
            return "Peak Heat"
        elif 15 <= hour < 19:
            return "Evening Cooling"
        else:
            return "Night Period"
    
    data['time_block'] = data['hour'].apply(get_time_block)
    
    block_counts = data['time_block'].value_counts()
    print(f"Time blocks assigned: {dict(block_counts)}")
    
    return data

def assign_thermal_phases(data: pd.DataFrame) -> pd.DataFrame:
    """Assign thermal phases based on temperature and solar radiation."""
    print("Assigning thermal phases...")
    
    def get_thermal_phase(temp, solar):
        if pd.isna(temp) or pd.isna(solar):
            return "Unknown"
        
        # Thresholds from thermal analysis
        temp_threshold = 3.3  # Mean temperature from thermal analysis
        solar_threshold = 100  # Solar radiation threshold
        
        if temp < temp_threshold and solar < solar_threshold:
            return "Cold Cloudy"
        elif temp < temp_threshold and solar >= solar_threshold:
            return "Cold Sunny"
        elif temp >= temp_threshold and solar < solar_threshold:
            return "Warm Cloudy"
        else:
            return "Warm Sunny"
    
    data['thermal_phase'] = data.apply(
        lambda row: get_thermal_phase(row['air_temp_C'], row['solar_radiation_Wpm2']), 
        axis=1
    )
    
    phase_counts = data['thermal_phase'].value_counts()
    print(f"Thermal phases assigned: {dict(phase_counts)}")
    
    return data

def create_combined_categories(data: pd.DataFrame) -> pd.DataFrame:
    """Create combined time-thermal category for analysis."""
    print("Creating combined categories...")
    
    # Create combined category
    data['time_thermal_combo'] = data['time_block'] + " × " + data['thermal_phase']
    
    # Summary statistics
    valid_deltas = data['delta_Ta'].notna().sum()
    total_records = len(data)
    
    print(f"Combined categories created for {total_records:,} records")
    print(f"Valid delta calculations: {valid_deltas:,} ({valid_deltas/total_records*100:.1f}%)")
    
    return data

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def integrated_pipeline_workflow() -> pd.DataFrame:
    """Complete pipeline from raw data to categorized output."""
    print("="*80)
    print("INTEGRATED DATA PROCESSING PIPELINE")
    print("="*80)
    
    # Load and clean data
    grow_2024, grow_2025, weather_combined = load_and_clean_data()
    
    # Merge with weather data
    integrated_2024 = merge_with_weather(grow_2024, weather_combined, "2024")
    integrated_2025 = merge_with_weather(grow_2025, weather_combined, "2025")
    
    # Combine years and filter to dormancy period
    all_integrated = pd.concat([integrated_2024, integrated_2025], ignore_index=True)
    dormancy_data = filter_dormancy_period(all_integrated)
    
    # Calculate temperature deltas
    data_with_deltas = calculate_temperature_deltas(dormancy_data)
    
    # Assign categories
    data_with_blocks = assign_time_blocks(data_with_deltas)
    data_with_phases = assign_thermal_phases(data_with_blocks)
    final_data = create_combined_categories(data_with_phases)
    
    # Save final output
    output_file = "grow_tube_categorized.csv"
    final_data.to_csv(output_file, index=False)
    
    print(f"\nPipeline complete!")
    print(f"Final dataset: {len(final_data):,} records")
    print(f"Saved: {output_file}")
    
    return final_data

# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    final_data = integrated_pipeline_workflow()