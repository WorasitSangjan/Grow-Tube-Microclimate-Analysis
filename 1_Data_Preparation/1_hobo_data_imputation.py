# Author: Worasit Sangjan
# Date: May 5, 2025

# ============================================================================
# HOBO DATA IMPUTATION
# Handles missing data imputation for grow tube temperature data
# ============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import time

warnings.filterwarnings('ignore')

# ============================================================================
# CORE FUNCTIONS (PRESERVED FROM ORIGINAL)
# ============================================================================

def load_hobo_data(file_path):
    """Load HOBO data with proper handling for your specific format."""
    print("Loading HOBO data...")
    try:
        df = pd.read_csv(file_path)
        print(f"Loaded {file_path}: {len(df)} rows")
    except Exception as e:
        raise ValueError(f"Could not load {file_path}: {e}")

    # Convert timestamp to datetime - handle MM/DD/YY H:MM format
    print("Converting timestamps...")
    df['Timestamp'] = pd.to_datetime(df['Timestamp'], format='%m/%d/%y %H:%M', errors='coerce')
    
    # Remove any rows with invalid timestamps
    initial_rows = len(df)
    df = df.dropna(subset=['Timestamp']).reset_index(drop=True)
    if len(df) < initial_rows:
        print(f"Removed {initial_rows - len(df)} rows with invalid timestamps")

    # Create Subject_ID efficiently
    df['Subject_ID'] = (
        df['Material'].astype(str) + '_' +
        df['Condition'].astype(str) + '_Rep' +
        df['Replicate'].astype(str)
    )

    # Sort by Subject_ID and Timestamp
    print("Sorting data...")
    df = df.sort_values(['Subject_ID', 'Timestamp']).reset_index(drop=True)

    print(f"Data loaded successfully!")
    print(f"Shape: {df.shape}")
    print(f"Date range: {df['Timestamp'].min()} to {df['Timestamp'].max()}")
    print(f"Unique Subject_IDs: {df['Subject_ID'].nunique()}")
    
    return df

def round_timestamp_to_5min_grid_vectorized(timestamps):
    """
    Vectorized timestamp rounding - much faster than loops.
    Round to nearest 5-minute mark on standard grid.
    """
    # Extract components for vectorized operations
    minutes = timestamps.dt.minute.values
    rounded_minutes = np.round(minutes / 5) * 5
    
    # Handle minute overflow efficiently
    hour_overflow = rounded_minutes >= 60
    rounded_minutes[hour_overflow] = 0
    
    # Create new timestamps efficiently
    result_timestamps = timestamps.copy()
    
    # Update minutes
    result_timestamps = result_timestamps.dt.floor('H') + pd.to_timedelta(rounded_minutes.astype(int), unit='m')
    
    # Handle hour overflow
    result_timestamps.loc[hour_overflow] += pd.Timedelta(hours=1)
    
    return result_timestamps

def create_complete_5min_grid_fast(df, tolerance_minutes=3):
    """Fast creation of complete 5-minute grid that handles timing shifts."""
    print(f"\n=== CREATING COMPLETE 5-MINUTE GRID ===")
    print(f"Tolerance for fuzzy matching: ±{tolerance_minutes} minutes")
    
    # Get all unique timestamps
    unique_timestamps = pd.Series(df['Timestamp'].unique()).sort_values()
    print(f"Original unique timestamps: {len(unique_timestamps):,}")
    
    # Round all timestamps to standard 5-minute grid (vectorized)
    print("Rounding timestamps to 5-minute grid...")
    rounded_timestamps = round_timestamp_to_5min_grid_vectorized(unique_timestamps)
    
    # Create mapping from rounded to original timestamps within tolerance
    print("Creating timestamp mapping...")
    timestamp_mapping = {}
    
    for orig_ts, rounded_ts in zip(unique_timestamps, rounded_timestamps):
        diff_minutes = abs((orig_ts - rounded_ts).total_seconds() / 60)
        if diff_minutes <= tolerance_minutes:
            if rounded_ts not in timestamp_mapping:
                timestamp_mapping[rounded_ts] = []
            timestamp_mapping[rounded_ts].append(orig_ts)
    
    # Create the complete expected grid
    start_time = rounded_timestamps.min()
    end_time = rounded_timestamps.max()
    
    # Generate complete 5-minute grid efficiently
    print("Generating complete 5-minute grid...")
    expected_timestamps = pd.date_range(start=start_time, end=end_time, freq='5min')
    
    print(f"Expected 5-minute grid: {len(expected_timestamps):,} timestamps")
    print(f"From {start_time} to {end_time}")
    
    # Identify which expected timestamps have data vs missing
    present_timestamps = expected_timestamps[expected_timestamps.isin(timestamp_mapping.keys())]
    missing_timestamps = expected_timestamps[~expected_timestamps.isin(timestamp_mapping.keys())]
    
    print(f"Present timestamps: {len(present_timestamps):,}")
    print(f"Missing timestamps: {len(missing_timestamps):,}")
    print(f"Missing percentage: {len(missing_timestamps)/len(expected_timestamps)*100:.1f}%")
    
    return expected_timestamps, present_timestamps, missing_timestamps, timestamp_mapping

def create_complete_dataset_fast(df, expected_timestamps, present_timestamps, missing_timestamps, timestamp_mapping):
    """Fast creation of complete dataset with NaN values for missing periods."""
    print(f"\n=== CREATING COMPLETE DATASET ===")
    
    all_subject_ids = df['Subject_ID'].unique()
    print(f"Creating complete dataset for {len(all_subject_ids)} subjects...")
    
    # Create subject info lookup for efficiency
    subject_info = df.groupby('Subject_ID')[['Material', 'Condition', 'Replicate']].first()
    
    # Create complete template efficiently
    print("Building complete template...")
    subjects_df = pd.DataFrame({'Subject_ID': all_subject_ids})
    timestamps_df = pd.DataFrame({'Timestamp': expected_timestamps})
    
    # Efficient cross join
    subjects_df['key'] = 1
    timestamps_df['key'] = 1
    complete_template = subjects_df.merge(timestamps_df, on='key').drop('key', axis=1)
    
    # Add subject info efficiently
    complete_template = complete_template.merge(subject_info, left_on='Subject_ID', right_index=True)
    
    print(f"Complete template created: {len(complete_template):,} rows")
    
    # Prepare actual data with rounded timestamps
    print("Preparing actual data with rounded timestamps...")
    df_with_rounded = df.copy()
    df_with_rounded['Rounded_Timestamp'] = round_timestamp_to_5min_grid_vectorized(df['Timestamp'])
    
    # Filter to only data within tolerance and aggregate duplicates
    valid_data_list = []
    for rounded_ts, orig_timestamps in timestamp_mapping.items():
        for orig_ts in orig_timestamps:
            mask = df_with_rounded['Timestamp'] == orig_ts
            if mask.any():
                data_subset = df_with_rounded[mask].copy()
                data_subset['Rounded_Timestamp'] = rounded_ts
                valid_data_list.append(data_subset)
    
    if valid_data_list:
        df_valid = pd.concat(valid_data_list, ignore_index=True)
        
        # Aggregate duplicates (if multiple measurements round to same 5-min mark)
        df_aggregated = df_valid.groupby(['Subject_ID', 'Rounded_Timestamp']).agg({
            'Value': 'mean',  # Take average if multiple values per 5-min period
            'Material': 'first',
            'Condition': 'first',
            'Replicate': 'first'
        }).reset_index()
        df_aggregated = df_aggregated.rename(columns={'Rounded_Timestamp': 'Timestamp'})
    else:
        df_aggregated = pd.DataFrame(columns=['Subject_ID', 'Timestamp', 'Value', 'Material', 'Condition', 'Replicate'])
    
    # Merge template with actual data
    print("Merging template with actual data...")
    complete_df = complete_template.merge(
        df_aggregated[['Subject_ID', 'Timestamp', 'Value']], 
        on=['Subject_ID', 'Timestamp'], 
        how='left'
    )
    
    # Sort for efficiency
    complete_df = complete_df.sort_values(['Subject_ID', 'Timestamp']).reset_index(drop=True)
    
    print(f"Data preparation complete")
    
    return complete_df

def cross_replicate_imputation_fast(df):
    """Fast cross-replicate imputation using vectorized operations."""
    print("\n=== CROSS-REPLICATE IMPUTATION ===")
    
    df_imputed = df.copy()
    
    # Create treatment group efficiently
    df_imputed['Treatment_Group'] = df_imputed['Material'] + '_' + df_imputed['Condition']
    
    # Get subjects with missing data
    subjects_with_missing = df_imputed[df_imputed['Value'].isna()]['Subject_ID'].unique()
    
    if len(subjects_with_missing) == 0:
        print("No missing data found for imputation.")
        return df_imputed.drop(['Treatment_Group'], axis=1)
    
    print(f"Processing {len(subjects_with_missing)} subjects with missing data...")
    
    # Process each treatment group for efficiency
    for treatment_group in df_imputed['Treatment_Group'].unique():
        treatment_data = df_imputed[df_imputed['Treatment_Group'] == treatment_group].copy()
        
        # Calculate mean values across replicates for each timestamp (vectorized)
        replicate_means = treatment_data.groupby('Timestamp')['Value'].mean()
        
        # Process each subject in this treatment group
        treatment_subjects = treatment_data['Subject_ID'].unique()
        subjects_with_missing_in_group = [s for s in treatment_subjects if s in subjects_with_missing]
        
        for subject_id in subjects_with_missing_in_group:
            subject_mask = (df_imputed['Subject_ID'] == subject_id)
            missing_mask = subject_mask & (df_imputed['Value'].isna())
            
            missing_count = missing_mask.sum()
            if missing_count == 0:
                continue
                
            # Get timestamps where this subject has missing data
            missing_timestamps = df_imputed.loc[missing_mask, 'Timestamp']
            
            # Fill with replicate means where available (vectorized)
            imputed_values = missing_timestamps.map(replicate_means)
            filled_count = imputed_values.notna().sum()
            
            # Update the dataframe
            df_imputed.loc[missing_mask, 'Value'] = imputed_values
    
    # Clean up
    df_imputed = df_imputed.drop(['Treatment_Group'], axis=1)
    
    return df_imputed

def visualize_imputation_results_fast(df_original, df_imputed, max_points_per_plot=2000):
    """Create fast visualizations with sampling for large datasets."""
    print("\n=== CREATING IMPUTATION VISUALIZATIONS ===")
    
    subjects = df_original['Subject_ID'].unique()
    
    fig, axes = plt.subplots(len(subjects), 1, figsize=(20, 5*len(subjects)))
    if len(subjects) == 1:
        axes = [axes]
    
    for i, subject_id in enumerate(subjects):
        ax = axes[i]
        
        # Get data for this subject
        orig_data = df_original[df_original['Subject_ID'] == subject_id].sort_values('Timestamp')
        imp_data = df_imputed[df_imputed['Subject_ID'] == subject_id].sort_values('Timestamp')
        
        # Sample for speed if too many points
        if len(orig_data) > max_points_per_plot:
            sample_indices = np.linspace(0, len(orig_data)-1, max_points_per_plot, dtype=int)
            orig_data = orig_data.iloc[sample_indices]
            imp_data = imp_data.iloc[sample_indices]
        
        # Plot original data (present values)
        orig_present = orig_data[orig_data['Value'].notna()]
        if len(orig_present) > 0:
            ax.plot(orig_present['Timestamp'], orig_present['Value'], 
                   'o-', alpha=0.7, label='Original Data', markersize=1, linewidth=0.8)
        
        # Plot imputed values
        orig_missing = orig_data[orig_data['Value'].isna()]
        if len(orig_missing) > 0:
            # Find which of these missing timestamps now have imputed values
            imputed_points = imp_data[
                imp_data['Timestamp'].isin(orig_missing['Timestamp']) & 
                imp_data['Value'].notna()
            ]
            if len(imputed_points) > 0:
                ax.plot(imputed_points['Timestamp'], imputed_points['Value'], 
                       'ro', alpha=0.8, label='Imputed Values', markersize=2)
        
        ax.set_title(f'{subject_id}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Date')
        ax.set_ylabel('Temperature (°C)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Format x-axis
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main_hobo_imputation_workflow_fast(file_path, tolerance_minutes=3, create_plots=True):
    """Complete workflow for HOBO data imputation."""
    start_time = time.time()
    
    print("="*80)
    print("HOBO DATA IMPUTATION WORKFLOW")
    print("="*80)
    
    # Step 1: Load data
    df = load_hobo_data(file_path)
    
    # Step 2: Create complete 5-minute grid
    expected_timestamps, present_timestamps, missing_timestamps, timestamp_mapping = create_complete_5min_grid_fast(
        df, tolerance_minutes
    )
    
    # Step 3: Create complete dataset
    complete_df = create_complete_dataset_fast(
        df, expected_timestamps, present_timestamps, missing_timestamps, timestamp_mapping
    )
    
    # Keep original for validation
    df_original = complete_df.copy()
    
    # Step 4: Apply cross-replicate imputation
    original_missing = complete_df['Value'].isna().sum()
    if original_missing > 0:
        print(f"\nApplying cross-replicate imputation for {original_missing:,} missing values...")
        df_imputed = cross_replicate_imputation_fast(complete_df)
    else:
        print("\nNo missing values found - dataset is complete!")
        df_imputed = complete_df
    
    # Step 5: Create visualizations
    if create_plots:
        visualize_imputation_results_fast(df_original, df_imputed)
    
    # Step 6: Save results
    output_file = "hobo_data_imputed_complete.csv"
    df_final = df_imputed.copy()
    df_final.to_csv(output_file, index=False)
    print(f"Complete imputed dataset saved to: {output_file}")
    
    # Step 7: Summary report
    end_time = time.time()
    duration = end_time - start_time
    
    print(f"\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    
    total_subjects = len(df['Subject_ID'].unique())
    total_timepoints = len(expected_timestamps)
    original_data_points = len(df)
    final_data_points = len(df_imputed[df_imputed['Value'].notna()])
    
    print(f"Original data points: {original_data_points:,}")
    print(f"Final data points: {final_data_points:,}")
    print(f"Data points added through imputation: {final_data_points - original_data_points:,}")
    print(f"Final completeness: {(final_data_points/(total_subjects * total_timepoints))*100:.1f}%")
    
    return df_imputed

# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    # Run the complete workflow
    file_path = "hobo_data_2025.csv"
    
    df_imputed = main_hobo_imputation_workflow_fast(
        file_path=file_path,
        tolerance_minutes=3,
        create_plots=True
    )
    
    print("\nWORKFLOW COMPLETED SUCCESSFULLY!")
