# Author: Worasit Sangjan
# Date: May 7, 2025

# ============================================================================
# DATA MERGING: HOBO and Weather Data Integration
# Merges HOBO temperature data with weather station data using 15-minute bins
# ============================================================================

import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def load_data(hobo_file: str, weather_file: str) -> tuple:
    """Load and prepare HOBO and weather data."""
    print("Loading data files...")
    
    hobo_data = pd.read_csv(hobo_file)
    weather_data = pd.read_csv(weather_file)
    
    # Convert timestamps
    hobo_data['Timestamp'] = pd.to_datetime(hobo_data['Timestamp'], errors='coerce')
    weather_data['timestamp'] = pd.to_datetime(weather_data['timestamp'], errors='coerce')
    
    # Remove invalid timestamps
    hobo_data = hobo_data.dropna(subset=['Timestamp'])
    weather_data = weather_data.dropna(subset=['timestamp'])
    
    print(f"HOBO data: {len(hobo_data):,} rows")
    print(f"Weather data: {len(weather_data):,} rows")
    
    return hobo_data, weather_data

def create_time_bins(hobo_data: pd.DataFrame, weather_data: pd.DataFrame) -> pd.DatetimeIndex:
    """Create 15-minute time bins covering both datasets."""
    print("Creating 15-minute time bins...")
    
    # Find common time range
    start = max(hobo_data['Timestamp'].min(), weather_data['timestamp'].min()).floor('15min')
    end = min(hobo_data['Timestamp'].max(), weather_data['timestamp'].max()).ceil('15min')
    
    # Create bin centers
    bin_centers = pd.date_range(start=start, end=end, freq='15min')
    
    print(f"Created {len(bin_centers):,} time bins from {start} to {end}")
    return bin_centers

def assign_to_bins(timestamps: pd.Series, bin_centers: pd.DatetimeIndex, tolerance_minutes: int = 7.5) -> pd.Series:
    """Assign timestamps to nearest bin center within tolerance."""
    tolerance_seconds = tolerance_minutes * 60
    
    # Find nearest bin for each timestamp
    idx = np.searchsorted(bin_centers, timestamps)
    closest_bins = []
    
    for i, ts in enumerate(timestamps):
        candidates = []
        
        # Check previous bin
        if idx[i] > 0:
            candidates.append(bin_centers[idx[i] - 1])
        
        # Check next bin
        if idx[i] < len(bin_centers):
            candidates.append(bin_centers[idx[i]])
        
        if candidates:
            # Find nearest candidate
            nearest = min(candidates, key=lambda x: abs((x - ts).total_seconds()))
            
            # Check if within tolerance
            if abs((nearest - ts).total_seconds()) <= tolerance_seconds:
                closest_bins.append(nearest)
            else:
                closest_bins.append(pd.NaT)
        else:
            closest_bins.append(pd.NaT)
    
    return pd.Series(closest_bins)

def aggregate_hobo_data(hobo_data: pd.DataFrame, bin_centers: pd.DatetimeIndex) -> pd.DataFrame:
    """Aggregate HOBO data by 15-minute bins."""
    print("Aggregating HOBO data to 15-minute bins...")
    
    # Assign bins
    hobo_data['bin'] = assign_to_bins(hobo_data['Timestamp'], bin_centers)
    
    # Remove unassigned data
    hobo_data = hobo_data.dropna(subset=['bin'])
    
    # Aggregate by bin and treatment
    hobo_agg = hobo_data.groupby(['bin', 'Material', 'Condition', 'Replicate']).agg({
        'Value': 'mean'
    }).reset_index().rename(columns={'bin': 'timestamp'})
    
    print(f"Aggregated to {len(hobo_agg):,} records")
    return hobo_agg

def aggregate_weather_data(weather_data: pd.DataFrame, bin_centers: pd.DatetimeIndex) -> pd.DataFrame:
    """Get closest weather record for each 15-minute bin."""
    print("Processing weather data...")
    
    # Assign bins
    weather_data['bin'] = assign_to_bins(weather_data['timestamp'], bin_centers)
    
    # Remove unassigned data
    weather_data = weather_data.dropna(subset=['bin'])
    
    # For each bin, get the closest weather record
    weather_closest = weather_data.loc[
        weather_data.groupby('bin')['timestamp'].apply(
            lambda g: (g - g.name).abs().idxmin()
        ).values
    ].reset_index(drop=True)
    
    # Clean up columns
    weather_closest = weather_closest.drop(columns='timestamp').rename(columns={'bin': 'timestamp'})
    
    print(f"Selected {len(weather_closest):,} weather records")
    return weather_closest

def merge_datasets(hobo_agg: pd.DataFrame, weather_agg: pd.DataFrame) -> pd.DataFrame:
    """Merge HOBO and weather data on timestamp."""
    print("Merging datasets...")
    
    merged = pd.merge(hobo_agg, weather_agg, on='timestamp', how='inner')
    
    print(f"Merged dataset: {len(merged):,} records")
    return merged

def validate_merge(merged_data: pd.DataFrame) -> None:
    """Validate the merged dataset."""
    print("\n=== MERGE VALIDATION ===")
    
    # Check for missing values
    missing_cols = merged_data.columns[merged_data.isna().any()].tolist()
    if missing_cols:
        print(f"Columns with missing values: {missing_cols}")
    else:
        print("No missing values")
    
    # Check time range
    time_range = merged_data['timestamp'].max() - merged_data['timestamp'].min()
    print(f"Time range: {time_range.days} days")
    
    # Check treatments
    treatments = merged_data.groupby(['Material', 'Condition']).size()
    print(f"Treatment combinations: {len(treatments)}")
    
    # Sample statistics
    print(f"Temperature range: {merged_data['Value'].min():.1f}°C to {merged_data['Value'].max():.1f}°C")

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def merge_hobo_weather_data(hobo_file: str = "hobo_data_2025_complete.csv", 
                           weather_file: str = "weather_data_2025.csv",
                           output_file: str = "merged_data_2025.csv") -> pd.DataFrame:
    """Complete data merging workflow."""
    print("="*60)
    print("HOBO-WEATHER DATA MERGING")
    print("="*60)
    
    # Load data
    hobo_data, weather_data = load_data(hobo_file, weather_file)
    
    # Create time bins
    bin_centers = create_time_bins(hobo_data, weather_data)
    
    # Aggregate both datasets
    hobo_agg = aggregate_hobo_data(hobo_data, bin_centers)
    weather_agg = aggregate_weather_data(weather_data, bin_centers)
    
    # Merge datasets
    merged = merge_datasets(hobo_agg, weather_agg)
    
    # Validate results
    validate_merge(merged)
    
    # Save output
    merged.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}")
    
    # Preview
    print(f"\n=== DATA PREVIEW ===")
    print(merged.head())
    
    return merged

# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    merged_data = merge_hobo_weather_data(
        hobo_file="hobo_data_2025_complete.csv",
        weather_file="weather_data_2025.csv",
        output_file="merged_data_2025.csv"
    )
