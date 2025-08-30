import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# =================================================================
# NEW DATA PREPARATION FUNCTIONS (FIXED FOR FULL SEASON)
# =================================================================

def prepare_data_from_raw_files(weather_file, treatment_file, season, output_dir, verbose=False):
    """
    Convert your raw data format into the format expected by your existing cold hardiness code.
    
    This creates individual CSV files that match the format your existing functions expect.
    """
    if verbose:
        print(f"PREPARING DATA FOR {season}")
        print("="*60)
    
    # Load weather data (Sep 1 - Oct 24)
    weather_data = pd.read_csv(weather_file)
    weather_data['datetime'] = pd.to_datetime(weather_data['timestamp'])
    weather_data['jday'] = weather_data['datetime'].dt.dayofyear
    weather_data['year'] = weather_data['datetime'].dt.year
    
    start_year = int(season.split('_')[0])
    weather_filtered = weather_data[
        (weather_data['year'] == start_year) &
        (weather_data['jday'] >= 244) & (weather_data['jday'] <= 297)  # Sep 1 - Oct 24
    ].copy()
    
    if verbose:
        print(f"Weather data: {len(weather_filtered)} records from {weather_filtered['datetime'].min()} to {weather_filtered['datetime'].max()}")
    
    # Load treatment data (Oct 25 onwards through March of NEXT year)
    treatment_data = pd.read_csv(treatment_file)
    treatment_data['datetime'] = pd.to_datetime(treatment_data['timestamp'])
    treatment_data['jday'] = treatment_data['datetime'].dt.dayofyear
    treatment_data['year'] = treatment_data['datetime'].dt.year
    
    # FIXED: Filter for treatment data from Oct 25 of start_year through March of next year
    next_year = start_year + 1
    treatment_filtered = treatment_data[
        ((treatment_data['year'] == start_year) & (treatment_data['jday'] >= 298)) |  # Oct 25 onwards in start year
        ((treatment_data['year'] == next_year) & (treatment_data['jday'] <= 90))      # Through Julian Day 90 in next year
    ].copy()
    
    # Create treatment identifier from Material + Condition
    treatment_filtered['Treatment'] = treatment_filtered['Material'] + '_' + treatment_filtered['Condition']
    
    if verbose:
        print(f"Treatment data: {len(treatment_filtered)} records from {treatment_filtered['datetime'].min()} to {treatment_filtered['datetime'].max()}")
        print(f"Treatments: {treatment_filtered['Treatment'].unique()}")
        print(f"Replicates: {treatment_filtered['Replicate'].unique()}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get unique treatment-replicate combinations
    combinations = treatment_filtered[['Treatment', 'Replicate']].drop_duplicates()
    
    created_files = []
    
    for _, row in combinations.iterrows():
        treatment = row['Treatment']
        replicate = row['Replicate']
        
        if verbose:
            print(f"Creating file for {treatment} replicate {replicate}")
        
        # Combine weather + treatment data for this combination
        combined_data = create_combined_data_for_treatment_replicate(
            weather_filtered, treatment_filtered, treatment, replicate, verbose=False
        )
        
        if combined_data is not None:
            # Create filename in the format your existing code expects
            filename = f"{treatment}_rep_{replicate}_complete.csv"
            filepath = os.path.join(output_dir, filename)
            
            # Save in the format your existing code expects
            combined_data.to_csv(filepath, index=False)
            created_files.append(filepath)
            
            if verbose:
                print(f"  ✓ Created: {filename}")
    
    if verbose:
        print(f"\nCreated {len(created_files)} files in {output_dir}")
        print("Files ready for your existing cold hardiness analysis!")
    
    return created_files

def create_combined_data_for_treatment_replicate(weather_data, treatment_data, treatment, replicate, verbose=False):
    """
    Combine weather data (early season) + treatment data (late season) for one treatment-replicate.
    Output format matches what your existing cold hardiness code expects.
    """
    
    # Get treatment data for this specific treatment and replicate
    treatment_subset = treatment_data[
        (treatment_data['Treatment'] == treatment) & 
        (treatment_data['Replicate'] == replicate)
    ].copy()
    
    if len(treatment_subset) == 0:
        if verbose:
            print(f"No data for {treatment} replicate {replicate}")
        return None
    
    # Prepare weather data (use air temperature as the temperature)
    weather_prepared = weather_data[['datetime', 'air_temp_C']].copy()
    weather_prepared.rename(columns={'air_temp_C': 'temperature'}, inplace=True)
    
    # Prepare treatment data (use Value as the temperature)
    treatment_prepared = treatment_subset[['datetime', 'Value']].copy()
    treatment_prepared.rename(columns={'Value': 'temperature'}, inplace=True)
    
    # Combine weather + treatment data
    combined = pd.concat([weather_prepared, treatment_prepared], ignore_index=True)
    combined = combined.sort_values('datetime').reset_index(drop=True)
    
    # Remove any duplicates (if there's overlap)
    combined = combined.drop_duplicates(subset=['datetime'], keep='last')
    
    if verbose:
        print(f"Combined data: {len(weather_prepared)} weather + {len(treatment_prepared)} treatment = {len(combined)} total")
    
    return combined

def run_cold_hardiness_analysis_with_data_prep(weather_file, treatment_file, season, varieties, output_base_dir):
    """
    Complete workflow: data preparation + cold hardiness analysis for ONE SEASON
    """
    print(f"COLD HARDINESS ANALYSIS WITH DATA PREPARATION")
    print(f"Season: {season}")
    print(f"Varieties: {varieties}")
    print("="*80)
    
    # Step 1: Prepare data (convert raw format to format expected by existing code)
    prepared_data_dir = os.path.join(output_base_dir, f"{season}_prepared_data")
    print(f"\nStep 1: Preparing data...")
    
    created_files = prepare_data_from_raw_files(
        weather_file=weather_file,
        treatment_file=treatment_file, 
        season=season,
        output_dir=prepared_data_dir,
        verbose=True
    )
    
    if len(created_files) == 0:
        print("No data files created. Check your input files.")
        return
    
    # Step 2: Run your existing cold hardiness analysis
    print(f"\nStep 2: Running cold hardiness analysis...")
    
    # Initialize your existing processor
    processor = GrowTubeColdHardnessProcessor()
    
    # Process each variety using your existing functions
    for variety in varieties:
        print(f"\n{'='*100}")
        print(f"PROCESSING VARIETY: {variety}")
        print(f"{'='*100}")
        
        # Create variety-specific output directory
        variety_output_dir = os.path.join(output_base_dir, f"{season}_{variety.replace(' ', '_')}")
        
        # Use your existing process_all_files function
        batch_results = processor.process_all_files(
            input_dir=prepared_data_dir,  # Use the prepared data
            variety=variety,
            output_dir=variety_output_dir,
            file_pattern="*_complete.csv"  # Files we just created
        )
        
        if batch_results and batch_results['successful_files'] > 0:
            # Use your existing combine_all_results function
            combined_data, summary_df = processor.combine_all_results(variety_output_dir, variety)
            
            # Use your existing plotting function
            if combined_data is not None:
                try:
                    plot_variety_treatment_comparison_with_weather(
                        combined_data, variety, variety_output_dir, 
                        weather_file=weather_file,  # Use original weather file
                        save_plots=True
                    )
                    print(f"   ✓ {variety} processing complete!")
                    print(f"   Individual files: {batch_results['successful_files']}")
                    print(f"   Combined data: {len(combined_data)} daily observations")
                    if summary_df is not None:
                        print(f"   Treatments: {summary_df['Treatment'].unique()}")
                except Exception as plot_error:
                    print(f"   Warning: Could not create plots - {plot_error}")
        else:
            print(f"   No successful processing for {variety}")
    
    print(f"\nALL PROCESSING COMPLETE!")
    print(f"Results saved in: {output_base_dir}")

# =================================================================
# FIXED FUNCTIONS - EXTENDED TO JULIAN DAY 90 (455)
# =================================================================

def find_julian_day_250(data, verbose=False):
    """Find Julian Day 250 in the data."""
    data['datetime'] = pd.to_datetime(data['datetime'])
    data['jday'] = data['datetime'].dt.dayofyear
    
    # Look for Julian Day 250 
    target_jday = 250
    
    # Find closest julian day to 250 in the data
    available_jdays = data['jday'].unique()
    closest_jday = min(available_jdays, key=lambda x: abs(x - target_jday))
    
    if verbose:
        print(f"    Using Julian Day: {closest_jday}")
    
    return closest_jday

def aggregate_to_daily(data, start_jday, verbose=False):
    """
    Aggregate 15-minute data to daily temperature values.
    Filter to start from Julian Day 250 and end at Julian Day 455 (90+365).
    
    Input: DataFrame with 15-min interval data  
    Output: DataFrame with daily T_mean, T_max, T_min, Acc_jday (Julian Day 250-455)
    """
    # Convert datetime and calculate julian day
    data['datetime'] = pd.to_datetime(data['datetime'])
    data['date'] = data['datetime'].dt.date
    data['jday'] = data['datetime'].dt.dayofyear
    
    if verbose:
        print(f"    Original data: {len(data)} rows, Julian Days {data['jday'].min()} to {data['jday'].max()}")
    
    # FIXED: Filter to dormant season using Julian Day reference
    # Julian Day 250-365 (Sep-Dec) + Julian Day 1-90 (Jan-Mar up to day 90)
    dormant_season_data = data[
        (data['jday'] >= start_jday) |  # From Julian Day 250 onwards (fall)
        (data['jday'] <= 90)           # Through Julian Day 90 (end at day 90)
    ].copy()
    
    if len(dormant_season_data) == 0:
        # Fall back to using all data but still apply the 90-day limit
        dormant_season_data = data[
            (data['jday'] >= 250) | (data['jday'] <= 90)
        ].copy()
    
    # DAILY AGGREGATION: Convert 15-minute intervals to daily values
    daily_data = dormant_season_data.groupby('date').agg({
        'temperature': ['mean', 'min', 'max', 'count'],  # Daily T_mean, T_min, T_max
        'datetime': 'first',  # Keep first datetime of the day
        'jday': 'first'      # Keep julian day
    }).reset_index()
    
    # Flatten column names
    daily_data.columns = ['date', 'T_mean', 'T_min', 'T_max', 'n_observations', 'datetime', 'jday']
    
    # Create Accumulate Julian Day
    # This handles year transition: 250, 251, ..., 365, 366, 367, ..., 455
    daily_data = daily_data.sort_values('datetime').reset_index(drop=True)
    daily_data['acc_jday'] = daily_data['jday'].copy()
    
    # Handle year transition in Julian Days
    for i in range(1, len(daily_data)):
        prev_jday = daily_data.loc[i-1, 'jday']
        curr_jday = daily_data.loc[i, 'jday']
        
        # If Julian Day jumped from high (Dec) to low (Jan), add 365
        if prev_jday > 300 and curr_jday < 100:
            daily_data.loc[i:, 'acc_jday'] = daily_data.loc[i:, 'jday'] + 365
    
    # FIXED: Filter accumulated Julian days to stop at 455 (day 90 + 365)
    daily_data = daily_data[daily_data['acc_jday'] <= 455].copy()
    
    # Reorder columns
    daily_data = daily_data[['datetime', 'T_mean', 'T_min', 'T_max', 'date', 'jday', 'acc_jday', 'n_observations']]
    
    if verbose:
        print(f"    Daily aggregation: {len(dormant_season_data)} 15-min obs → {len(daily_data)} daily values")
        print(f"    Accumulated Julian Day range: {daily_data['acc_jday'].min()} to {daily_data['acc_jday'].max()}")
        print(f"    Date range: {daily_data['date'].min()} to {daily_data['date'].max()}")
        print(f"    Temperature range: {daily_data['T_mean'].min():.1f} to {daily_data['T_mean'].max():.1f}°C")
    
    return daily_data

def apply_model(daily_data, variety, file_info, variety_params, verbose=False):
    """
    Apply cold hardiness model to daily temperature data.
    1. Calculate daily degree days
    2. Calculate acclimation/deacclimation  
    3. Update model_Hc
    4. Update DD_chilling_sum
    5. Update base10_chilling_sum
    6. Update DD_heating_sum (if dormancy_period == 2) - BEFORE transition check
    7. Check dormancy transition
    """
    if variety not in variety_params:
        available = list(variety_params.keys())
        raise ValueError(f"Variety '{variety}' not available. Choose from: {available}")
    
    params = variety_params[variety]
    
    # Initialize results
    results = daily_data.copy()
    results['Variety'] = variety
    results['Treatment'] = file_info['treatment']
    results['Replicate'] = file_info['replicate']
    results['File'] = file_info['filename']
    results['Predicted_Hc'] = np.nan
    results['Dormancy_period'] = 1
    results['DD_heating_sum'] = 0.0
    results['DD_chilling_sum'] = 0.0
    results['Base10_chilling_sum'] = 0.0
    
    # Model parameters
    Hc_min = params['Hc_min']
    Hc_max = params['Hc_max']
    Hc_range = Hc_min - Hc_max
    
    # Parameter arrays
    T_threshold = [None, params['T_threshold_endo'], params['T_threshold_eco']]
    acclimation_rate = [None, params['Acclimation_rate_endo'], params['Acclimation_rate_eco']]
    deacclimation_rate = [None, params['Deacclimation_rate_endo'], params['Deacclimation_rate_eco']]
    theta = [None, 1, params['Theta']]  # theta(1)=1, theta(2)=variety_theta
    
    # Initialize variables
    model_Hc_yesterday = params['Hc_initial']
    dormancy_period = 1
    DD_heating_sum = 0.0
    DD_chilling_sum = 0.0
    base10_chilling_sum = 0.0
    ecodormancy_boundary = params['Ecodormancy_boundary']
    budbreak_date = None
    
    if verbose:
        print(f"    Running a model for {variety}")
        print(f"    Starting Hc: {model_Hc_yesterday}°C, Theta: {params['Theta']}")
    
    # Main model loop - process each daily value
    for i, row in results.iterrows():
        T_mean = row['T_mean']
        acc_jday = row['acc_jday']  # Use Julian Day as reference
        
        # Skip if missing temperature data (If T_mean = "" Then End)
        if pd.isna(T_mean):
            if verbose and i < 5:  # Only show first few warnings
                print(f"      Skipping Julian Day {acc_jday}: missing temperature")
            continue
        
        # === VBA CALCULATION ORDER ===
        
        # 1. Calculate degree days for today
        if T_mean > T_threshold[dormancy_period]:
            DD_heating_today = T_mean - T_threshold[dormancy_period]
            DD_chilling_today = 0.0
        else:
            DD_heating_today = 0.0
            DD_chilling_today = T_mean - T_threshold[dormancy_period]
        
        # Calculate base-10 chilling
        if T_mean <= 10:
            base10_chilling_today = T_mean - 10
        else:
            base10_chilling_today = 0.0
        
        # 2. Calculate acclimation and deacclimation
        
        # DEACCLIMATION 
        deacclimation = (DD_heating_today * deacclimation_rate[dormancy_period] * 
                       (1 - ((model_Hc_yesterday - Hc_max) / Hc_range) ** theta[dormancy_period]))
        
        # Deacclimation restriction: do not allow deacclimation unless some chilling has occurred
        if DD_chilling_sum == 0:
            deacclimation = 0.0
        
        # ACCLIMATION
        acclimation = (DD_chilling_today * acclimation_rate[dormancy_period] * 
                      (1 - (Hc_min - model_Hc_yesterday) / Hc_range))
        
        # 3. Update model_Hc
        Delta_Hc = acclimation + deacclimation
        model_Hc = model_Hc_yesterday + Delta_Hc
        
        # Bounds checking
        if model_Hc <= Hc_max:
            model_Hc = Hc_max
        if model_Hc > Hc_min:
            model_Hc = Hc_min
        
        # 4. Update cumulative chilling sums
        DD_chilling_sum += DD_chilling_today
        base10_chilling_sum += base10_chilling_today
        
        # 5. Update heating sum ONLY if in ecodormancy AND BEFORE dormancy transition check
        if dormancy_period == 2:
            DD_heating_sum += DD_heating_today
        
        # 6. Check dormancy transition - AFTER heating sum update
        if base10_chilling_sum <= ecodormancy_boundary:
            dormancy_period = 2
        
        # 7. Budbreak detection
        if budbreak_date is None:
            if Hc_min == -1.2:  # Vinifera
                if model_Hc_yesterday < -2.2 and model_Hc >= -2.2:
                    budbreak_date = row['datetime']
            elif Hc_min == -2.5:  # Labruscana
                if model_Hc_yesterday < -6.4 and model_Hc >= -6.4:
                    budbreak_date = row['datetime']
        
        # Store results
        results.loc[i, 'Predicted_Hc'] = model_Hc
        results.loc[i, 'Dormancy_period'] = dormancy_period
        results.loc[i, 'DD_heating_sum'] = DD_heating_sum
        results.loc[i, 'DD_chilling_sum'] = DD_chilling_sum
        results.loc[i, 'Base10_chilling_sum'] = base10_chilling_sum
        results.loc[i, 'ACC_jday'] = acc_jday
        
        # Remember today's hardiness for tomorrow
        model_Hc_yesterday = model_Hc
    
    results['Budbreak_date'] = budbreak_date
    
    if verbose:
        final_hc = results['Predicted_Hc'].dropna().iloc[-1] if len(results['Predicted_Hc'].dropna()) > 0 else np.nan
        min_hc = results['Predicted_Hc'].min()
        n_days = len(results['Predicted_Hc'].dropna())
        print(f"    Processed {n_days} days: Final Hc: {final_hc:.2f}°C, Min Hc: {min_hc:.2f}°C")
        if budbreak_date:
            print(f"    Budbreak predicted: {budbreak_date.strftime('%Y-%m-%d')}")
        
        # Show dormancy transition info
        if 'Dormancy_period' in results.columns:
            endo_days = len(results[results['Dormancy_period'] == 1])
            eco_days = len(results[results['Dormancy_period'] == 2])
            print(f"    Dormancy periods: {endo_days} days endodormancy, {eco_days} days ecodormancy")
    
    return results

def parse_filename(filename):
    """Parse grow tube experiment filename."""
    basename = os.path.basename(filename).replace('_complete.csv', '')
    parts = basename.split('_')
    
    if len(parts) >= 4:
        material = parts[0]
        height = parts[1]
        replicate = int(parts[3])
        treatment = f"{material}_{height}"
    else:
        # Fallback parsing
        material = 'Unknown'
        height = 'Unknown'
        replicate = 1
        treatment = basename
    
    return {
        'filename': filename,
        'material': material,
        'height': height,
        'replicate': replicate,
        'treatment': treatment,
        'basename': basename
    }

def calculate_file_summary(results, file_info):
    """Calculate summary statistics for a single file."""
    valid_hc = results['Predicted_Hc'].dropna()
    
    if len(valid_hc) == 0:
        return {
            'filename': file_info['filename'],
            'treatment': file_info['treatment'],
            'replicate': file_info['replicate'],
            'final_hc': np.nan,
            'min_hc': np.nan,
            'max_hc': np.nan,
            'n_days': 0,
            'budbreak_date': None,
            'dormancy_transition_date': None
        }
    
    # Find dormancy transition
    dormancy_transition = None
    if 'Dormancy_period' in results.columns:
        transition_idx = results[results['Dormancy_period'] == 2].index
        if len(transition_idx) > 0:
            dormancy_transition = results.loc[transition_idx[0], 'datetime']
    
    return {
        'filename': file_info['filename'],
        'treatment': file_info['treatment'],
        'replicate': file_info['replicate'],
        'final_hc': valid_hc.iloc[-1],
        'min_hc': valid_hc.min(),
        'max_hc': valid_hc.max(),
        'n_days': len(valid_hc),
        'budbreak_date': results['Budbreak_date'].iloc[0] if 'Budbreak_date' in results.columns else None,
        'dormancy_transition_date': dormancy_transition,
        'date_range_start': results['datetime'].min(),
        'date_range_end': results['datetime'].max()
    }

def plot_cold_hardiness_results(results, variety, file_info, output_dir, save_plots=True):
    """Create individual file cold hardiness plots with single y-axis for temperature (Julian Day 250-455)."""
    # Prepare data
    results['datetime'] = pd.to_datetime(results['datetime'])
    results = results.sort_values('datetime').reset_index(drop=True)
    
    # Create single figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Temperature lines
    ax.plot(results['acc_jday'], results['T_max'], 'r-', linewidth=1, 
            marker='o', markersize=2, label='T_max', alpha=0.7)
    ax.plot(results['acc_jday'], results['T_min'], 'b-', linewidth=1, 
            marker='o', markersize=2, label='T_min', alpha=0.7)
    
    # Cold hardiness line
    ax.plot(results['acc_jday'], results['Predicted_Hc'], 'k-', linewidth=1, label='Predicted Hc')
    
    # Add dormancy period markers
    if 'Dormancy_period' in results.columns:
        endo_data = results[results['Dormancy_period'] == 1]
        eco_data = results[results['Dormancy_period'] == 2]
        
        if len(endo_data) > 0:
            ax.scatter(endo_data['acc_jday'], endo_data['Predicted_Hc'], 
                      c='blue', s=15, alpha=0.6, marker='x', label='Endodormancy')
        if len(eco_data) > 0:
            ax.scatter(eco_data['acc_jday'], eco_data['Predicted_Hc'], 
                      c='red', s=15, alpha=0.6, marker='x', label='Ecodormancy')
    
    # Mark budbreak with green star
    budbreak_date = results['Budbreak_date'].iloc[0]
    if pd.notna(budbreak_date):
        budbreak_row = results[results['datetime'] == budbreak_date]
        if len(budbreak_row) > 0:
            ax.scatter(budbreak_row['acc_jday'].iloc[0], budbreak_row['Predicted_Hc'].iloc[0], 
                      c='green', s=150, marker='*', label='Budbreak', 
                      zorder=5, edgecolors='black', linewidth=1)
    
    # Set up axis - FIXED to show full range
    ax.set_xlabel('Julian Day', fontsize=14, fontweight='bold')
    ax.set_ylabel('Temperature (°C)', fontsize=14, fontweight='bold')
    ax.set_xlim(250, 455)
    ax.set_ylim(-30, 40)
    ax.grid(True, alpha=0.3)
    
    # FIXED: Custom x-axis labels to include Julian Day 90
    ax.set_xticks([250, 275, 300, 325, 350, 365, 390, 415, 440, 455])
    ax.set_xticklabels(['250', '275', '300', '325', '350', '365', '25', '50', '75', '90'])
    
    # Title
    budbreak_str = budbreak_date.strftime('%m/%d/%Y') if pd.notna(budbreak_date) else 'N/A'
    title = f"{variety}, {file_info['treatment']}, Budbreak predicted on {budbreak_str}"
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add legend
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    # Save plot
    if save_plots:
        plot_file = os.path.join(output_dir, f"{file_info['basename']}_{variety}_plot.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"    Individual plot saved: {plot_file}")
    
    return fig

def plot_variety_treatment_comparison_with_weather(combined_data, variety, output_dir, 
                                                   weather_file='weather_2024_2025_binned.csv', 
                                                   save_plots=True):
    """Create variety-treatment comparison plot with ambient weather Tmax/Tmin lines."""
    if combined_data is None or len(combined_data) == 0:
        print("No data available for variety treatment comparison plot")
        return None
    
    # Define colors matching the second graph
    treatment_colors = {
        'Paper_low': '#1f77b4',      # Blue (like RDVI)
        'Paper_raise': '#ff7f0e',    # Orange (like NDVI) 
        'Plastic_low': '#2ca02c',    # Green (like CI_GREEN)
        'Plastic_raise': '#d62728',  # Red (like GNDVI)
        'Uncover_uncover': '#9467bd' # Purple (like EVI2)
    }
    
    # Load and process ambient weather data
    try:
        print(f"Loading ambient weather data from {weather_file}...")
        weather_data = pd.read_csv(weather_file)
        
        # Parse weather timestamps 
        weather_data['datetime'] = pd.to_datetime(weather_data['timestamp'])
        weather_data['date'] = weather_data['datetime'].dt.date
        weather_data['jday'] = weather_data['datetime'].dt.dayofyear
        
        # Create accumulated julian day (handle year transition)
        weather_data = weather_data.sort_values('datetime').reset_index(drop=True)
        weather_data['acc_jday'] = weather_data['jday'].copy()
        
        for i in range(1, len(weather_data)):
            prev_jday = weather_data.loc[i-1, 'jday']
            curr_jday = weather_data.loc[i, 'jday']
            
            # If Julian Day jumped from high (Dec) to low (Jan), add 365
            if prev_jday > 300 and curr_jday < 100:
                weather_data.loc[i:, 'acc_jday'] = weather_data.loc[i:, 'jday'] + 365
        
        # Calculate daily Tmax and Tmin from ambient air temperature
        daily_weather = weather_data.groupby('acc_jday')['air_temp_C'].agg(['min', 'max']).reset_index()
        daily_weather.columns = ['acc_jday', 'ambient_tmin', 'ambient_tmax']
        
        # FIXED: Filter to Julian Day 250-455 range
        daily_weather = daily_weather[
            (daily_weather['acc_jday'] >= 250) & 
            (daily_weather['acc_jday'] <= 455)
        ].copy()
        
        print(f"Processed ambient weather: {len(daily_weather)} days")
        
    except Exception as e:
        print(f"Warning: Could not load weather data - {e}")
        daily_weather = None
    
    # Create figure with single y-axis
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # === PLOT COLD HARDINESS LINES ===
    treatments = sorted(combined_data['Treatment'].unique())
    
    for treatment in treatments:
        treatment_data = combined_data[combined_data['Treatment'] == treatment]
        
        # Calculate treatment means for cold hardiness
        treatment_stats = treatment_data.groupby('acc_jday')['Predicted_Hc'].mean().reset_index()
        treatment_stats = treatment_stats.sort_values('acc_jday')
        
        # Get color for this treatment
        color = treatment_colors.get(treatment, '#000000')  # Default to black if not found
        
        ax.plot(treatment_stats['acc_jday'], treatment_stats['Predicted_Hc'], 
               color=color, linewidth=3, label=f"{treatment}", 
               marker='o', markersize=4, alpha=0.9)
    
    # === PLOT AMBIENT TEMPERATURE LINES ===
    if daily_weather is not None and len(daily_weather) > 0:
        ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmax'], 
               color='red', linestyle='--', linewidth=2, alpha=0.7, 
               label='Ambient Tmax')
        
        ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmin'], 
               color='blue', linestyle='--', linewidth=2, alpha=0.7, 
               label='Ambient Tmin')
        
        print(f"Added ambient temperature lines: {daily_weather['ambient_tmin'].min():.1f}°C to {daily_weather['ambient_tmax'].max():.1f}°C")
    
    # === FORMATTING ===
    
    # Single y-axis for all temperature data
    ax.set_xlabel('Julian Day', fontsize=30, fontweight='bold')
    ax.set_ylabel('Temperature (°C)', fontsize=30, fontweight='bold')
    
    # FIXED: Title to show correct range
    ax.set_title(f'{variety}', 
                fontsize=37, fontweight='bold')
    
    # FIXED: X-axis setup to show full range
    ax.set_xlim(240, 455)
    ax.set_xticks([250, 275, 300, 325, 350, 365, 390, 415, 440, 455])
    ax.set_xticklabels(['250', '275', '300', '325', '350', '365', '25', '50', '75', '90'])
    ax.tick_params(axis='both', labelsize=28)

    # Y-axis limits - combine both cold hardiness and ambient temperature ranges
    if len(combined_data) > 0:
        hc_min = combined_data['Predicted_Hc'].min()
        hc_max = combined_data['Predicted_Hc'].max()
        
        if daily_weather is not None and len(daily_weather) > 0:
            temp_min = min(hc_min, daily_weather['ambient_tmin'].min())
            temp_max = max(hc_max, daily_weather['ambient_tmax'].max())
        else:
            temp_min = hc_min
            temp_max = hc_max
        
        ax.set_ylim(temp_min - 5, temp_max + 5)
    
    # Grid
    ax.grid(True, alpha=0.3)
    
    # Single legend combining all lines
    ax.legend(loc='best', ncol=4, fontsize=18)
    
    plt.tight_layout()
    
    if save_plots:
        plot_file = os.path.join(output_dir, f"{variety}_variety_treatment_comparison.png")
        plt.savefig(plot_file, dpi=400, bbox_inches='tight')
        print(f"Variety treatment comparison plot saved: {plot_file}")
    
    return fig

def plot_all_varieties_comparison(varieties_data, output_dir, season, 
                                 weather_file='weather_2024_2025_binned.csv', 
                                 save_plots=True):
    """
    Create a single plot comparing cold hardiness across all varieties.
    """
    # Define colors for each variety
    variety_colors = {
        'Chardonnay': '#1f77b4',      # Blue
        'Concord': '#ff7f0e',         # Orange
        'Cabernet Sauvignon': '#2ca02c',  # Green
        'Mourvedre': '#d62728',       # Red
        'Barbera': '#9467bd',         # Purple
        'Cabernet franc': '#8c564b',  # Brown
        'Chenin blanc': '#e377c2',    # Pink
        'Dolcetto': '#7f7f7f',        # Gray
        'Gewürztraminer': '#bcbd22',  # Olive
        'Grenache': '#17becf'         # Cyan
    }
    
    # Load and process ambient weather data
    try:
        print(f"Loading ambient weather data from {weather_file}...")
        weather_data = pd.read_csv(weather_file)
        
        # Parse weather timestamps 
        weather_data['datetime'] = pd.to_datetime(weather_data['timestamp'])
        weather_data['date'] = weather_data['datetime'].dt.date
        weather_data['jday'] = weather_data['datetime'].dt.dayofyear
        
        # Create accumulated julian day (handle year transition)
        weather_data = weather_data.sort_values('datetime').reset_index(drop=True)
        weather_data['acc_jday'] = weather_data['jday'].copy()
        
        for i in range(1, len(weather_data)):
            prev_jday = weather_data.loc[i-1, 'jday']
            curr_jday = weather_data.loc[i, 'jday']
            
            # If Julian Day jumped from high (Dec) to low (Jan), add 365
            if prev_jday > 300 and curr_jday < 100:
                weather_data.loc[i:, 'acc_jday'] = weather_data.loc[i:, 'jday'] + 365
        
        # Calculate daily Tmax and Tmin from ambient air temperature
        daily_weather = weather_data.groupby('acc_jday')['air_temp_C'].agg(['min', 'max']).reset_index()
        daily_weather.columns = ['acc_jday', 'ambient_tmin', 'ambient_tmax']
        
        # Filter to Julian Day 250-455 range
        daily_weather = daily_weather[
            (daily_weather['acc_jday'] >= 250) & 
            (daily_weather['acc_jday'] <= 455)
        ].copy()
        
        print(f"Processed ambient weather: {len(daily_weather)} days")
        
    except Exception as e:
        print(f"Warning: Could not load weather data - {e}")
        daily_weather = None
    
    # Create figure with high resolution
    fig, ax = plt.subplots(figsize=(16, 10), dpi=150)
    
    # Plot cold hardiness for each variety
    for variety_name, combined_data in varieties_data.items():
        if combined_data is None or len(combined_data) == 0:
            print(f"Warning: No data for {variety_name}")
            continue
            
        # Get all treatments for this variety and calculate overall mean
        variety_mean = combined_data.groupby('acc_jday')['Predicted_Hc'].mean().reset_index()
        variety_mean = variety_mean.sort_values('acc_jday')
        
        # Get color for this variety
        color = variety_colors.get(variety_name, '#000000')  # Default to black if not found
        
        # Plot variety mean cold hardiness
        ax.plot(variety_mean['acc_jday'], variety_mean['Predicted_Hc'], 
               color=color, linewidth=3, label=f"{variety_name}", 
               marker='o', markersize=3, alpha=0.9)
        
        print(f"Added {variety_name}: {len(variety_mean)} data points")
    
    # Plot ambient temperature lines
    if daily_weather is not None and len(daily_weather) > 0:
        ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmax'], 
               color='red', linestyle='--', linewidth=2, alpha=0.7, 
               label='Ambient Tmax')
        
        ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmin'], 
               color='blue', linestyle='--', linewidth=2, alpha=0.7, 
               label='Ambient Tmin')
        
        print(f"Added ambient temperature lines: {daily_weather['ambient_tmin'].min():.1f}°C to {daily_weather['ambient_tmax'].max():.1f}°C")
    
    # Formatting
    ax.set_xlabel('Julian Day', fontsize=16, fontweight='bold')
    ax.set_ylabel('Temperature (°C)', fontsize=16, fontweight='bold')
    
    # Title
    ax.set_title(f'Cold Hardiness Comparison - All Varieties ({season})\nJulian Day 250-455', 
                fontsize=18, fontweight='bold', pad=20)
    
    # X-axis setup
    ax.set_xlim(240, 460)
    ax.set_xticks([250, 275, 300, 325, 350, 365, 390, 415, 440, 455])
    ax.set_xticklabels(['250\n(Sep 7)', '275\n(Oct 2)', '300\n(Oct 27)', 
                       '325\n(Nov 21)', '350\n(Dec 16)', '365\n(Dec 31)',
                       '25\n(Jan 25)', '50\n(Feb 19)', '75\n(Mar 16)', '90\n(Mar 31)'],
                       fontsize=12)
    
    # Y-axis limits
    if len(varieties_data) > 0:
        all_hc_values = []
        for combined_data in varieties_data.values():
            if combined_data is not None and len(combined_data) > 0:
                all_hc_values.extend(combined_data['Predicted_Hc'].dropna().values)
        
        if len(all_hc_values) > 0:
            hc_min = min(all_hc_values)
            hc_max = max(all_hc_values)
            
            if daily_weather is not None and len(daily_weather) > 0:
                temp_min = min(hc_min, daily_weather['ambient_tmin'].min())
                temp_max = max(hc_max, daily_weather['ambient_tmax'].max())
            else:
                temp_min = hc_min
                temp_max = hc_max
            
            ax.set_ylim(temp_min - 5, temp_max + 5)
    
    # Grid and legend
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=12, ncol=2)
    
    # Increase tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    plt.tight_layout()
    
    # Save plot
    if save_plots:
        plot_file = os.path.join(output_dir, f"all_varieties_comparison_{season}.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"All varieties comparison plot saved: {plot_file}")
    
    return fig

def plot_varieties_by_treatment(varieties_data, output_dir, season, 
                                weather_file='weather_2024_2025_binned.csv', 
                                save_plots=True):
    """
    Create 4 separate plots, one for each treatment, showing all varieties under that treatment.
    """
    # Define colors for each variety (consistent across all plots)
    variety_colors = {
        'Chardonnay': '#1f77b4',      # Blue
        'Concord': '#ff7f0e',         # Orange
        'Cabernet Sauvignon': '#2ca02c',  # Green
        'Mourvedre': '#d62728',       # Red
        'Barbera': '#9467bd',         # Purple
        'Cabernet franc': '#8c564b',  # Brown
        'Chenin blanc': '#e377c2',    # Pink
        'Dolcetto': '#7f7f7f',        # Gray
        'Gewürztraminer': '#bcbd22',  # Olive
        'Grenache': '#17becf'         # Cyan
    }
    
    # Get list of treatments from the data
    all_treatments = set()
    for variety_name, combined_data in varieties_data.items():
        if combined_data is not None and len(combined_data) > 0:
            all_treatments.update(combined_data['Treatment'].unique())
    
    treatments = sorted(list(all_treatments))
    print(f"Found treatments: {treatments}")
    
    # Load and process ambient weather data
    try:
        print(f"Loading ambient weather data from {weather_file}...")
        weather_data = pd.read_csv(weather_file)
        
        # Parse weather timestamps 
        weather_data['datetime'] = pd.to_datetime(weather_data['timestamp'])
        weather_data['date'] = weather_data['datetime'].dt.date
        weather_data['jday'] = weather_data['datetime'].dt.dayofyear
        
        # Create accumulated julian day (handle year transition)
        weather_data = weather_data.sort_values('datetime').reset_index(drop=True)
        weather_data['acc_jday'] = weather_data['jday'].copy()
        
        for i in range(1, len(weather_data)):
            prev_jday = weather_data.loc[i-1, 'jday']
            curr_jday = weather_data.loc[i, 'jday']
            
            # If Julian Day jumped from high (Dec) to low (Jan), add 365
            if prev_jday > 300 and curr_jday < 100:
                weather_data.loc[i:, 'acc_jday'] = weather_data.loc[i:, 'jday'] + 365
        
        # Calculate daily Tmax and Tmin from ambient air temperature
        daily_weather = weather_data.groupby('acc_jday')['air_temp_C'].agg(['min', 'max']).reset_index()
        daily_weather.columns = ['acc_jday', 'ambient_tmin', 'ambient_tmax']
        
        # Filter to Julian Day 250-455 range
        daily_weather = daily_weather[
            (daily_weather['acc_jday'] >= 250) & 
            (daily_weather['acc_jday'] <= 455)
        ].copy()
        
        print(f"Processed ambient weather: {len(daily_weather)} days")
        
    except Exception as e:
        print(f"Warning: Could not load weather data - {e}")
        daily_weather = None
    
    # Create one plot for each treatment
    created_plots = []
    
    for treatment in treatments:
        print(f"\nCreating plot for treatment: {treatment}")
        
        # Create figure with high resolution
        fig, ax = plt.subplots(figsize=(16, 10), dpi=150)
        
        varieties_plotted = 0
        
        # Plot each variety's performance under this treatment
        for variety_name, combined_data in varieties_data.items():
            if combined_data is None or len(combined_data) == 0:
                continue
                
            # Filter data for this specific treatment
            treatment_data = combined_data[combined_data['Treatment'] == treatment]
            
            if len(treatment_data) == 0:
                print(f"  No data for {variety_name} under {treatment}")
                continue
            
            # Calculate mean cold hardiness for this variety under this treatment
            variety_treatment_mean = treatment_data.groupby('acc_jday')['Predicted_Hc'].mean().reset_index()
            variety_treatment_mean = variety_treatment_mean.sort_values('acc_jday')
            
            # Get color for this variety
            color = variety_colors.get(variety_name, '#000000')
            
            # Plot variety performance under this treatment
            ax.plot(variety_treatment_mean['acc_jday'], variety_treatment_mean['Predicted_Hc'], 
                   color=color, linewidth=3, label=f"{variety_name}", 
                   marker='o', markersize=3, alpha=0.9)
            
            varieties_plotted += 1
            print(f"  Added {variety_name}: {len(variety_treatment_mean)} data points")
        
        # Only create plot if we have data for this treatment
        if varieties_plotted == 0:
            plt.close(fig)
            print(f"  Skipping {treatment} - no variety data found")
            continue
        
        # Plot ambient temperature lines
        if daily_weather is not None and len(daily_weather) > 0:
            ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmax'], 
                   color='red', linestyle='--', linewidth=2, alpha=0.7, 
                   label='Ambient Tmax')
            
            ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmin'], 
                   color='blue', linestyle='--', linewidth=2, alpha=0.7, 
                   label='Ambient Tmin')
        
        # Formatting
        ax.set_xlabel('Julian Day', fontsize=35, fontweight='bold')
        ax.set_ylabel('Temperature (°C)', fontsize=35, fontweight='bold')
        
        # Title
        treatment_title = treatment.replace('_', ' ').title()
        ax.set_title(f'{treatment_title}', 
                    fontsize=40, fontweight='bold', pad=20)
        
        # X-axis setup
        ax.set_xlim(240, 460)
        ax.set_xticks([250, 275, 300, 325, 350, 365, 390, 415, 440, 455])
        ax.set_xticklabels(['250)', '275', '300', '325', '350', '365', '25', '50', '75', '90'])


        # Y-axis limits - calculate from all varieties under this treatment
        all_hc_values = []
        for variety_name, combined_data in varieties_data.items():
            if combined_data is not None and len(combined_data) > 0:
                treatment_data = combined_data[combined_data['Treatment'] == treatment]
                if len(treatment_data) > 0:
                    all_hc_values.extend(treatment_data['Predicted_Hc'].dropna().values)
        
        if len(all_hc_values) > 0:
            hc_min = min(all_hc_values)
            hc_max = max(all_hc_values)
            
            if daily_weather is not None and len(daily_weather) > 0:
                temp_min = min(hc_min, daily_weather['ambient_tmin'].min())
                temp_max = max(hc_max, daily_weather['ambient_tmax'].max())
            else:
                temp_min = hc_min
                temp_max = hc_max
            
            ax.set_ylim(temp_min - 5, temp_max + 5)
        
        # Grid and legend
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=25, ncol=3)
        
        # Increase tick label sizes
        ax.tick_params(axis='both', which='major', labelsize=32)
        
        plt.tight_layout()
        
        # Save plot
        if save_plots:
            plot_file = os.path.join(output_dir, f"varieties_comparison_{treatment}_{season}.png")
            plt.savefig(plot_file, dpi=400, bbox_inches='tight', facecolor='white')
            print(f"  Treatment comparison plot saved: {plot_file}")
            created_plots.append(plot_file)
        
        plt.close(fig)
    
    print(f"\nCreated {len(created_plots)} treatment comparison plots")
    return created_plots


def run_analysis_with_treatment_comparisons(weather_file, treatment_file, season, 
                                          varieties, output_base_dir):
    """
    Complete workflow: analyze all varieties and create treatment-based comparison plots.
    """
    print(f"COMPLETE ANALYSIS WITH TREATMENT-BASED COMPARISONS FOR {season}")
    print(f"Varieties: {varieties}")
    print("="*80)
    
    # Store results for each variety
    varieties_data = {}
    
    # Run analysis for each variety
    for variety in varieties:
        print(f"\n{'='*100}")
        print(f"PROCESSING VARIETY: {variety}")
        print(f"{'='*100}")
        
        # Step 1: Prepare data for this variety
        prepared_data_dir = os.path.join(output_base_dir, f"{season}_prepared_data")
        
        if not os.path.exists(prepared_data_dir):
            print(f"Preparing data for {season}...")
            prepare_data_from_raw_files(
                weather_file=weather_file,
                treatment_file=treatment_file, 
                season=season,
                output_dir=prepared_data_dir,
                verbose=True
            )
        
        # Step 2: Process this variety
        processor = GrowTubeColdHardnessProcessor()
        variety_output_dir = os.path.join(output_base_dir, f"{season}_{variety.replace(' ', '_')}")
        
        batch_results = processor.process_all_files(
            input_dir=prepared_data_dir,
            variety=variety,
            output_dir=variety_output_dir,
            file_pattern="*_complete.csv"
        )
        
        # Step 3: Combine results for this variety
        if batch_results and batch_results['successful_files'] > 0:
            combined_data, summary_df = processor.combine_all_results(variety_output_dir, variety)
            
            # Store for comparison plots
            varieties_data[variety] = combined_data
            
            # Create individual variety plot
            if combined_data is not None:
                try:
                    plot_variety_treatment_comparison_with_weather(
                        combined_data, variety, variety_output_dir, 
                        weather_file=weather_file,
                        save_plots=True
                    )
                    print(f"   ✓ {variety} processing complete!")
                except Exception as plot_error:
                    print(f"   Warning: Could not create individual plot - {plot_error}")
        else:
            print(f"   No successful processing for {variety}")
            varieties_data[variety] = None
    
    # Step 4: Create treatment-based comparison plots
    print(f"\n{'='*100}")
    print("CREATING TREATMENT-BASED VARIETY COMPARISON PLOTS")
    print(f"{'='*100}")
    
    # Filter out varieties with no data
    valid_varieties_data = {k: v for k, v in varieties_data.items() if v is not None and len(v) > 0}
    
    if len(valid_varieties_data) > 0:
        try:
            created_plots = plot_varieties_by_treatment(
                varieties_data=valid_varieties_data,
                output_dir=output_base_dir,
                season=season,
                weather_file=weather_file,
                save_plots=True
            )
            print(f"✓ Treatment-based comparison plots created!")
            print(f"   Varieties included: {list(valid_varieties_data.keys())}")
            print(f"   Plots created: {len(created_plots)}")
        except Exception as e:
            print(f"Error creating treatment comparison plots: {e}")
    else:
        print("No valid variety data available for treatment comparison plots")
    
    print(f"\n🎉 COMPLETE ANALYSIS FINISHED!")
    print(f"Results saved in: {output_base_dir}")
    
    return varieties_data


def run_all_varieties_analysis_and_comparison(weather_file, treatment_file, season, 
                                            varieties, output_base_dir):
    """
    Complete workflow: analyze all varieties and create comparison plot.
    """
    print(f"COMPLETE MULTI-VARIETY ANALYSIS FOR {season}")
    print(f"Varieties: {varieties}")
    print("="*80)
    
    # Store results for each variety
    varieties_data = {}
    
    # Run analysis for each variety
    for variety in varieties:
        print(f"\n{'='*100}")
        print(f"PROCESSING VARIETY: {variety}")
        print(f"{'='*100}")
        
        # Step 1: Prepare data for this variety
        prepared_data_dir = os.path.join(output_base_dir, f"{season}_prepared_data")
        
        if not os.path.exists(prepared_data_dir):
            print(f"Preparing data for {season}...")
            prepare_data_from_raw_files(
                weather_file=weather_file,
                treatment_file=treatment_file, 
                season=season,
                output_dir=prepared_data_dir,
                verbose=True
            )
        
        # Step 2: Process this variety
        processor = GrowTubeColdHardnessProcessor()
        variety_output_dir = os.path.join(output_base_dir, f"{season}_{variety.replace(' ', '_')}")
        
        batch_results = processor.process_all_files(
            input_dir=prepared_data_dir,
            variety=variety,
            output_dir=variety_output_dir,
            file_pattern="*_complete.csv"
        )
        
        # Step 3: Combine results for this variety
        if batch_results and batch_results['successful_files'] > 0:
            combined_data, summary_df = processor.combine_all_results(variety_output_dir, variety)
            
            # Store for comparison plot
            varieties_data[variety] = combined_data
            
            # Create individual variety plot
            if combined_data is not None:
                try:
                    plot_variety_treatment_comparison_with_weather(
                        combined_data, variety, variety_output_dir, 
                        weather_file=weather_file,
                        save_plots=True
                    )
                    print(f"   ✓ {variety} processing complete!")
                except Exception as plot_error:
                    print(f"   Warning: Could not create individual plot - {plot_error}")
        else:
            print(f"   No successful processing for {variety}")
            varieties_data[variety] = None
    
    # Step 4: Create comparison plot with all varieties
    print(f"\n{'='*100}")
    print("CREATING MULTI-VARIETY COMPARISON PLOT")
    print(f"{'='*100}")
    
    # Filter out varieties with no data
    valid_varieties_data = {k: v for k, v in varieties_data.items() if v is not None and len(v) > 0}
    
    if len(valid_varieties_data) > 0:
        try:
            plot_all_varieties_comparison(
                varieties_data=valid_varieties_data,
                output_dir=output_base_dir,
                season=season,
                weather_file=weather_file,
                save_plots=True
            )
            print(f"✓ Multi-variety comparison plot created!")
            print(f"   Varieties included: {list(valid_varieties_data.keys())}")
        except Exception as e:
            print(f"Error creating comparison plot: {e}")
    else:
        print("No valid variety data available for comparison plot")
    
    print(f"\n🎉 COMPLETE ANALYSIS FINISHED!")
    print(f"Results saved in: {output_base_dir}")
    
    return varieties_data

def combine_all_results(output_dir, variety):
    """Combine all individual result files into comprehensive summaries."""
    print(f"\nCOMBINING ALL RESULTS FOR {variety}")
    
    # Find all result files
    result_files = glob.glob(os.path.join(output_dir, f"*_{variety}_results.csv"))
    
    if len(result_files) == 0:
        print(f"No result files found for {variety}")
        return None, None
    
    print(f"Found {len(result_files)} result files to combine")
    
    # Load and combine all results
    all_data = []
    summary_data = []
    
    for result_file in result_files:
        try:
            data = pd.read_csv(result_file)
            all_data.append(data)
            
            # Extract summary for this file
            if len(data) > 0:
                valid_hc = data['Predicted_Hc'].dropna()
                if len(valid_hc) > 0:
                    summary_data.append({
                        'Variety': variety,
                        'Treatment': data['Treatment'].iloc[0],
                        'Replicate': data['Replicate'].iloc[0],
                        'File': data['File'].iloc[0],
                        'Final_Hc': valid_hc.iloc[-1],
                        'Min_Hc': valid_hc.min(),
                        'Max_Hc': valid_hc.max(),
                        'N_Days': len(valid_hc),
                        'Budbreak_Date': data['Budbreak_date'].iloc[0],
                        'Start_Date': data['datetime'].min(),
                        'End_Date': data['datetime'].max()
                    })
        
        except Exception as e:
            print(f"Error reading {result_file}: {e}")
    
    # Combine all data
    if len(all_data) > 0:
        combined_data = pd.concat(all_data, ignore_index=True)
        
        # Save combined detailed results
        combined_file = os.path.join(output_dir, f"combined_detailed_results_{variety}.csv")
        combined_data.to_csv(combined_file, index=False)
        print(f"    Saved combined detailed results: {combined_file}")
        
        # Create and save summary
        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(output_dir, f"combined_summary_{variety}.csv")
        summary_df.to_csv(summary_file, index=False)
        print(f"    Saved combined summary: {summary_file}")
        
        # Print summary statistics
        print(f"\nSUMMARY FOR {variety}:")
        if len(summary_df) > 0:
            by_treatment = summary_df.groupby('Treatment').agg({
                'Final_Hc': ['mean', 'std', 'count'],
                'Min_Hc': ['mean', 'std']
            }).round(2)
            
            print("Treatment Summary:")
            print(by_treatment)
        
        return combined_data, summary_df
    
    else:
        print("No data could be combined")
        return None, None

def process_all_files(input_dir, variety, output_dir, variety_params, file_pattern="*_complete.csv"):
    """Process all grow tube experiment files for a specific variety."""
    print(f"PROCESSING ALL FILES FOR {variety.upper()}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"File pattern: {file_pattern}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all files
    file_pattern_path = os.path.join(input_dir, file_pattern)
    files = glob.glob(file_pattern_path)
    files.sort()
    
    if len(files) == 0:
        print(f"No files found matching pattern: {file_pattern_path}")
        return None
    
    print(f"Found {len(files)} files to process")
    
    # Process each file
    all_results = []
    successful_files = 0
    failed_files = 0
    
    for filepath in files:
        result = process_single_file(filepath, variety, output_dir, variety_params, verbose=True)
        all_results.append(result)
        
        if result['success']:
            successful_files += 1
        else:
            failed_files += 1
    
    print(f"\n{'='*80}")
    print(f"PROCESSING COMPLETE FOR {variety}")
    print(f"Successful: {successful_files}, Failed: {failed_files}")
    print(f"{'='*80}")
    
    return {
        'variety': variety,
        'input_dir': input_dir,
        'output_dir': output_dir,
        'total_files': len(files),
        'successful_files': successful_files,
        'failed_files': failed_files,
        'results': all_results
    }

def process_single_file(filepath, variety, output_dir, variety_params, verbose=True):
    """Process a single grow tube experiment file."""
    filename = os.path.basename(filepath)
    file_info = parse_filename(filename)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Processing: {filename}")
        print(f"Treatment: {file_info['treatment']}, Replicate: {file_info['replicate']}")
        print(f"Variety: {variety}")
        print(f"{'='*60}")
    
    try:
        # Step 1: Load 15-minute data
        if verbose:
            print("1. Loading 15-minute data...")
        data = pd.read_csv(filepath)
        
        if verbose:
            print(f"   Loaded {len(data)} rows, columns: {list(data.columns)}")
        
        # Step 2: Aggregate to daily values
        if verbose:
            print("2. Aggregating to daily values...")
        start_jday = find_julian_day_250(data, verbose=verbose)
        daily_data = aggregate_to_daily(data, start_jday, verbose=verbose)
        
        # Step 3: Apply cold hardiness model
        if verbose:
            print("3. Applying cold hardiness model...")
        results = apply_model(daily_data, variety, file_info, variety_params, verbose=verbose)
        
        # Step 4: Calculate summary statistics
        summary = calculate_file_summary(results, file_info)
        
        # Step 5: Save individual results and create plots
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{file_info['basename']}_{variety}_results.csv")
        results.to_csv(output_file, index=False)
        
        # Create plots
        try:
            plot_cold_hardiness_results(results, variety, file_info, output_dir, save_plots=True)
        except Exception as plot_error:
            print(f"Warning: Could not create plot - {plot_error}")
        
        if verbose:
            print(f"4. Results saved: {output_file}")
            print(f"   Summary: Final Hc = {summary['final_hc']:.2f}°C, Min Hc = {summary['min_hc']:.2f}°C")
        
        return {
            'success': True,
            'filepath': filepath,
            'output_file': output_file,
            'file_info': file_info,
            'results': results,
            'summary': summary,
            'variety': variety
        }
        
    except Exception as e:
        if verbose:
            print(f"Error processing {filename}: {e}")
        return {
            'success': False,
            'filepath': filepath,
            'file_info': file_info,
            'error': str(e),
            'variety': variety
        }

# =================================================================
# CLASS DEFINITION (PRESERVED EXACTLY)
# =================================================================

class GrowTubeColdHardnessProcessor:
    """
    Process individual grow tube experiment files for cold hardiness analysis.
    
    Workflow:
    1. Load individual CSV file (e.g., Paper_low_rep_1_complete.csv)
    2. Aggregate 15-min data to daily T_mean, T_max, T_min
    3. Apply cold hardiness model
    4. Save individual results
    5. Combine all files for comprehensive analysis
    """
    
    def __init__(self, start_date_target="Julian Day 250", end_date_flexible=True):
        """Initialize with the model parameters and processing settings."""
        
        # PROCESSING SETTINGS
        self.start_date_target = start_date_target  # Model starts Julian Day 250
        self.end_date_flexible = end_date_flexible  # OK to end at end of March
        self.acc_start_julian_target = 250  #  Julian Day 250
        
        # MODEL PARAMETERS (by variety) - PRESERVED EXACTLY FROM YOUR ORIGINAL CODE
        self.variety_params = {
            'Barbera': {
                'Hc_initial': -10.1, 'Hc_min': -1.2, 'Hc_max': -23.5,
                'T_threshold_endo': 15.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.06, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.10, 'Deacclimation_rate_eco': 0.08, 
                'Theta': 7  # Used as theta(2) in the model, theta(1) always = 1
            },
            'Cabernet franc': {
                'Hc_initial': -9.9, 'Hc_min': -1.2, 'Hc_max': -25.4,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 4.0, 'Ecodormancy_boundary': -500,
                'Acclimation_rate_endo': 0.12, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.04, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 7
            },
            'Cabernet Sauvignon': {
                'Hc_initial': -10.3, 'Hc_min': -1.2, 'Hc_max': -25.1,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 5.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.12, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.08, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 7
            },
            'Chardonnay': {
                'Hc_initial': -11.8, 'Hc_min': -1.2, 'Hc_max': -25.7,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -600,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.10, 'Deacclimation_rate_eco': 0.08, 
                'Theta': 7
            },
            'Chenin blanc': {
                'Hc_initial': -12.1, 'Hc_min': -1.2, 'Hc_max': -24.1,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 4.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.04, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 7
            },
            'Concord': {
                'Hc_initial': -12.8, 'Hc_min': -2.5, 'Hc_max': -29.5,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -600,
                'Acclimation_rate_endo': 0.12, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 3  # Lower theta = less curved deacclimation in ecodormancy
            },
            'Dolcetto': {
                'Hc_initial': -10.1, 'Hc_min': -1.2, 'Hc_max': -23.2,
                'T_threshold_endo': 12.0, 'T_threshold_eco': 4.0, 'Ecodormancy_boundary': -600,
                'Acclimation_rate_endo': 0.16, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.10, 'Deacclimation_rate_eco': 0.12, 
                'Theta': 3
            },
            'Gewürztraminer': {
                'Hc_initial': -11.6, 'Hc_min': -1.2, 'Hc_max': -24.9,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 6.0, 'Ecodormancy_boundary': -400,
                'Acclimation_rate_endo': 0.12, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.06, 'Deacclimation_rate_eco': 0.18, 
                'Theta': 5
            },
            'Grenache': {
                'Hc_initial': -10.0, 'Hc_min': -1.2, 'Hc_max': -22.7,
                'T_threshold_endo': 12.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -500,
                'Acclimation_rate_endo': 0.16, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.06, 
                'Theta': 5
            },
            'Lemberger': {
                'Hc_initial': -13.0, 'Hc_min': -1.2, 'Hc_max': -25.6,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 5.0, 'Ecodormancy_boundary': -800,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.18, 
                'Theta': 7
            },
            'Malbec': {
                'Hc_initial': -11.5, 'Hc_min': -1.2, 'Hc_max': -25.1,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 4.0, 'Ecodormancy_boundary': -400,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.08,
                'Deacclimation_rate_endo': 0.06, 'Deacclimation_rate_eco': 0.08, 
                'Theta': 7
            },
            'Merlot': {
                'Hc_initial': -10.3, 'Hc_min': -1.2, 'Hc_max': -25.0,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 5.0, 'Ecodormancy_boundary': -500,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.04, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 7
            },
            'Mourvedre': {
                'Hc_initial': -9.5, 'Hc_min': -1.2, 'Hc_max': -22.1,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 6.0, 'Ecodormancy_boundary': -600,
                'Acclimation_rate_endo': 0.12, 'Acclimation_rate_eco': 0.06,
                'Deacclimation_rate_endo': 0.08, 'Deacclimation_rate_eco': 0.14, 
                'Theta': 5  # Moderate curve in ecodormancy deacclimation
            },
            'Nebbiolo': {
                'Hc_initial': -11.1, 'Hc_min': -1.2, 'Hc_max': -24.4,
                'T_threshold_endo': 11.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.16, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 3
            },
            'Pinot gris': {
                'Hc_initial': -12.0, 'Hc_min': -1.2, 'Hc_max': -24.1,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 6.0, 'Ecodormancy_boundary': -400,
                'Acclimation_rate_endo': 0.12, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.20, 
                'Theta': 3
            },
            'Riesling': {
                'Hc_initial': -12.6, 'Hc_min': -1.2, 'Hc_max': -26.1,
                'T_threshold_endo': 12.0, 'T_threshold_eco': 5.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.14, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.12, 
                'Theta': 7
            },
            'Sangiovese': {
                'Hc_initial': -10.7, 'Hc_min': -1.2, 'Hc_max': -21.9,
                'T_threshold_endo': 11.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.14, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.06, 
                'Theta': 7
            },
            'Sauvignon blanc': {
                'Hc_initial': -10.6, 'Hc_min': -1.2, 'Hc_max': -24.9,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 5.0, 'Ecodormancy_boundary': -300,
                'Acclimation_rate_endo': 0.08, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.06, 'Deacclimation_rate_eco': 0.12, 
                'Theta': 7
            },
            'Sémillon': {
                'Hc_initial': -10.4, 'Hc_min': -1.2, 'Hc_max': -22.4,
                'T_threshold_endo': 13.0, 'T_threshold_eco': 7.0, 'Ecodormancy_boundary': -300,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.02,
                'Deacclimation_rate_endo': 0.08, 'Deacclimation_rate_eco': 0.20, 
                'Theta': 5
            },
            'Sunbelt': {
                'Hc_initial': -11.8, 'Hc_min': -2.5, 'Hc_max': -29.1,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -400,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.06, 'Deacclimation_rate_eco': 0.12, 
                'Theta': 1.5  # Very low theta = almost linear deacclimation
            },
            'Syrah': {
                'Hc_initial': -10.3, 'Hc_min': -1.2, 'Hc_max': -24.2,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 4.0, 'Ecodormancy_boundary': -700,
                'Acclimation_rate_endo': 0.08, 'Acclimation_rate_eco': 0.04,
                'Deacclimation_rate_endo': 0.06, 'Deacclimation_rate_eco': 0.08, 
                'Theta': 7
            },
            'Viognier': {
                'Hc_initial': -11.2, 'Hc_min': -1.2, 'Hc_max': -24.0,
                'T_threshold_endo': 14.0, 'T_threshold_eco': 5.0, 'Ecodormancy_boundary': -300,
                'Acclimation_rate_endo': 0.10, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.08, 'Deacclimation_rate_eco': 0.10, 
                'Theta': 7
            },
            'Zinfandel': {
                'Hc_initial': -10.4, 'Hc_min': -1.2, 'Hc_max': -24.4,
                'T_threshold_endo': 12.0, 'T_threshold_eco': 3.0, 'Ecodormancy_boundary': -500,
                'Acclimation_rate_endo': 0.16, 'Acclimation_rate_eco': 0.10,
                'Deacclimation_rate_endo': 0.02, 'Deacclimation_rate_eco': 0.06, 
                'Theta': 7
            }
        }
        
        # MODEL CALCULATION SETTINGS
        # Theta(1) = 1 (endodormancy), theta(2) = variety Theta (ecodormancy)
        self.theta_endodormancy = 1 
        
        # Budbreak detection thresholds
        self.budbreak_threshold_vinifera = -2.2  # Most varieties
        self.budbreak_threshold_labruscana = -6.4  # Concord, Sunbelt
        
        print(f"   Grow Tube Cold Hardiness Processor Initialized")
        print(f"   Target start: {self.start_date_target}")
        print(f"   End flexibility: {self.end_date_flexible}")
        print(f"   Available varieties: {len(self.variety_params)}")
        print(f"   Model settings: theta(1)=1, theta(2)=variety-specific")
    
    def get_available_varieties(self):
        """Get list of available grape varieties."""
        return list(self.variety_params.keys())
    
    def process_single_file(self, filepath, variety, output_dir, verbose=True):
        """Process a single grow tube experiment file using class methods."""
        return process_single_file(filepath, variety, output_dir, self.variety_params, verbose)
    
    def process_all_files(self, input_dir, variety, output_dir, file_pattern="*_complete.csv"):
        """Process all grow tube experiment files for a specific variety using class methods."""
        return process_all_files(input_dir, variety, output_dir, self.variety_params, file_pattern)
    
    def combine_all_results(self, output_dir, variety):
        """Combine all individual result files into comprehensive summaries."""
        return combine_all_results(output_dir, variety)

# =================================================================
# MAIN EXECUTION SCRIPT - ONE SEASON AT A TIME
# =================================================================

if __name__ == "__main__":
    
    print("GROW TUBE COLD HARDINESS ANALYSIS - ONE SEASON AT A TIME")
    print("Weather Data (Sep 1 - Oct 24) + Treatment Data (Oct 25 onwards through March)")
    print("="*80)
    
    # =================================================================
    # CONFIGURATION - UPDATE THESE FOR YOUR FILES
    # =================================================================
    
    # Your file paths
    TREATMENT_FILE = "grow_tube_categorized.csv"  # Your combined treatment data file
    
    # Choose which season to run (run one at a time)
    SEASON = "2023_2024"  # Change to "2024_2025" for the other season
    
    # Weather file for the chosen season
    if SEASON == "2023_2024":
        WEATHER_FILE = "weather_2023_2024_binned.csv"
    elif SEASON == "2024_2025":
        WEATHER_FILE = "weather_2024_2025_binned.csv"
    else:
        print(f"Unknown season: {SEASON}")
        exit(1)
    
    # Analysis settings
    VARIETIES = ['Chardonnay', 'Concord', 'Cabernet Sauvignon', 'Mourvedre']  # Start with just a couple, add more as needed
    OUTPUT_DIR = f"cold_hardiness_results_{SEASON}"
    
    # =================================================================
    # FILE VALIDATION
    # =================================================================
    
    print(f"Season: {SEASON}")
    print(f"Weather file: {WEATHER_FILE}")
    print(f"Treatment file: {TREATMENT_FILE}")
    print(f"Varieties: {VARIETIES}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("-" * 60)
    
    # Check if files exist
    if not os.path.exists(WEATHER_FILE):
        print(f"❌ Weather file not found: {WEATHER_FILE}")
        print("Please check the file name and location")
        exit(1)
    else:
        print(f"✅ Weather file found: {WEATHER_FILE}")
    
    if not os.path.exists(TREATMENT_FILE):
        print(f"❌ Treatment file not found: {TREATMENT_FILE}")
        print("Please check the file name and location")
        exit(1)
    else:
        print(f"✅ Treatment file found: {TREATMENT_FILE}")
    
    # =================================================================
    # RUN ANALYSIS
    # =================================================================
    
    print(f"\nStarting analysis for {SEASON}...")
    
    try:
        run_analysis_with_treatment_comparisons(
            weather_file=WEATHER_FILE,
            treatment_file=TREATMENT_FILE,
            season=SEASON,
            varieties=VARIETIES,
            output_base_dir=OUTPUT_DIR
        )
        
        print(f"\n🎉 ANALYSIS COMPLETE FOR {SEASON}!")
        print(f"Results saved in: {OUTPUT_DIR}")
        print(f"\nTo run the other season:")
        print(f"1. Change SEASON to '2024_2025' if you ran '2023_2024' (or vice versa)")
        print(f"2. Run the script again")
        
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        print("\nPlease check:")
        print("1. File paths are correct")
        print("2. Data files have the expected format")
        print("3. The season exists in your data")

    # =================================================================
    # OPTIONAL: QUICK DATA CHECK
    # =================================================================
    
    print(f"\n" + "="*60)
    print("QUICK DATA CHECK")
    print("="*60)
    
    try:
        # Quick check of treatment data
        treatment_data = pd.read_csv(TREATMENT_FILE)
        print(f"Treatment data: {len(treatment_data)} total records")
        
        # Check seasons in data
        treatment_data['datetime'] = pd.to_datetime(treatment_data['timestamp'])
        treatment_data['year'] = treatment_data['datetime'].dt.year
        seasons_in_data = treatment_data['year'].unique()
        print(f"Years in treatment data: {seasons_in_data}")
        
        # Check current season
        start_year = int(SEASON.split('_')[0])
        season_data = treatment_data[treatment_data['year'] == start_year]
        print(f"Records for {SEASON}: {len(season_data)}")
        
        if len(season_data) > 0:
            # Check treatments and replicates
            season_data['Treatment'] = season_data['Material'] + '_' + season_data['Condition']
            treatments = season_data['Treatment'].unique()
            replicates = season_data['Replicate'].unique()
            print(f"Treatments in {SEASON}: {treatments}")
            print(f"Replicates in {SEASON}: {replicates}")
            
            # Check date range
            print(f"Date range in {SEASON}: {season_data['datetime'].min()} to {season_data['datetime'].max()}")
        else:
            print(f"⚠️  No data found for {SEASON}")
    
    except Exception as e:
        print(f"Error in data check: {e}")