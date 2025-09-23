# Author: Worasit Sangjan
# Date: July 5, 2025

# ============================================================================
# COLD HARDINESS ANALYSIS FOR GRAPE CULTIVAS
# Processes grow tube experiment data and creates comparison visualizations
# ============================================================================

import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

# Set plotting style for consistency
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 30
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 3.0

# ============================================================================
# CORE DATA PROCESSING FUNCTIONS
# ============================================================================

def prepare_data_from_raw_files(weather_file: str, treatment_file: str, season: str, 
                               output_dir: str, verbose: bool = False) -> list:
    """Convert raw data format into format expected by cold hardiness analysis."""
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
        (weather_data['jday'] >= 244) & (weather_data['jday'] <= 297)
    ].copy()
    
    # Load treatment data (Oct 25 onwards through March of next year)
    treatment_data = pd.read_csv(treatment_file)
    treatment_data['datetime'] = pd.to_datetime(treatment_data['timestamp'])
    treatment_data['jday'] = treatment_data['datetime'].dt.dayofyear
    treatment_data['year'] = treatment_data['datetime'].dt.year
    
    next_year = start_year + 1
    treatment_filtered = treatment_data[
        ((treatment_data['year'] == start_year) & (treatment_data['jday'] >= 298)) |
        ((treatment_data['year'] == next_year) & (treatment_data['jday'] <= 90))
    ].copy()
    
    # Create treatment identifier and apply mapping
    treatment_filtered['Treatment'] = treatment_filtered['Material'] + '_' + treatment_filtered['Condition']
    treatment_name_mapping = {
        'Paper_low': 'Paper Low',
        'Paper_raise': 'Paper High', 
        'Plastic_low': 'Plastic Low',
        'Plastic_raise': 'Plastic High',
        'Uncover_uncover': 'No-Tube'
    }
    treatment_filtered['Treatment'] = treatment_filtered['Treatment'].map(treatment_name_mapping).fillna(treatment_filtered['Treatment'])
    
    if verbose:
        print(f"Weather data: {len(weather_filtered)} records")
        print(f"Treatment data: {len(treatment_filtered)} records")
        print(f"Treatments: {treatment_filtered['Treatment'].unique()}")
    
    # Create output files
    os.makedirs(output_dir, exist_ok=True)
    combinations = treatment_filtered[['Treatment', 'Replicate']].drop_duplicates()
    created_files = []
    
    for _, row in combinations.iterrows():
        treatment, replicate = row['Treatment'], row['Replicate']
        combined_data = create_combined_data_for_treatment_replicate(
            weather_filtered, treatment_filtered, treatment, replicate
        )
        
        if combined_data is not None:
            filename = f"{treatment.replace(' ', '_').replace('-', '_')}_rep_{replicate}_complete.csv"
            filepath = os.path.join(output_dir, filename)
            combined_data.to_csv(filepath, index=False)
            created_files.append(filepath)
    
    if verbose:
        print(f"Created {len(created_files)} files in {output_dir}")
    
    return created_files

def create_combined_data_for_treatment_replicate(weather_data: pd.DataFrame, treatment_data: pd.DataFrame, 
                                               treatment: str, replicate: int) -> pd.DataFrame:
    """Combine weather and treatment data for one treatment-replicate combination."""
    treatment_subset = treatment_data[
        (treatment_data['Treatment'] == treatment) & 
        (treatment_data['Replicate'] == replicate)
    ].copy()
    
    if len(treatment_subset) == 0:
        return None
    
    # Prepare data with consistent column names
    weather_prepared = weather_data[['datetime', 'air_temp_C']].copy()
    weather_prepared.rename(columns={'air_temp_C': 'temperature'}, inplace=True)
    
    treatment_prepared = treatment_subset[['datetime', 'Value']].copy()
    treatment_prepared.rename(columns={'Value': 'temperature'}, inplace=True)
    
    # Combine and clean
    combined = pd.concat([weather_prepared, treatment_prepared], ignore_index=True)
    combined = combined.sort_values('datetime').reset_index(drop=True)
    combined = combined.drop_duplicates(subset=['datetime'], keep='last')
    
    return combined

def aggregate_to_daily(data: pd.DataFrame, start_jday: int, verbose: bool = False) -> pd.DataFrame:
    """Aggregate 15-minute data to daily temperature values."""
    data['datetime'] = pd.to_datetime(data['datetime'])
    data['date'] = data['datetime'].dt.date
    data['jday'] = data['datetime'].dt.dayofyear
    
    # Filter to dormant season
    dormant_season_data = data[
        (data['jday'] >= start_jday) | (data['jday'] <= 90)
    ].copy()
    
    # Daily aggregation
    daily_data = dormant_season_data.groupby('date').agg({
        'temperature': ['mean', 'min', 'max', 'count'],
        'datetime': 'first',
        'jday': 'first'
    }).reset_index()
    
    # Flatten column names
    daily_data.columns = ['date', 'T_mean', 'T_min', 'T_max', 'n_observations', 'datetime', 'jday']
    
    # Create accumulated Julian day for year transition
    daily_data = daily_data.sort_values('datetime').reset_index(drop=True)
    daily_data['acc_jday'] = daily_data['jday'].copy()
    
    for i in range(1, len(daily_data)):
        prev_jday = daily_data.loc[i-1, 'jday']
        curr_jday = daily_data.loc[i, 'jday']
        if prev_jday > 300 and curr_jday < 100:
            daily_data.loc[i:, 'acc_jday'] = daily_data.loc[i:, 'jday'] + 365
    
    daily_data = daily_data[daily_data['acc_jday'] <= 455].copy()
    daily_data = daily_data[['datetime', 'T_mean', 'T_min', 'T_max', 'date', 'jday', 'acc_jday', 'n_observations']]
    
    if verbose:
        print(f"Daily aggregation: {len(dormant_season_data)} 15-min obs → {len(daily_data)} daily values")
    
    return daily_data

def apply_model(daily_data: pd.DataFrame, variety: str, file_info: dict, 
               variety_params: dict, verbose: bool = False) -> pd.DataFrame:
    """Apply cold hardiness model to daily temperature data."""
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
    
    T_threshold = [None, params['T_threshold_endo'], params['T_threshold_eco']]
    acclimation_rate = [None, params['Acclimation_rate_endo'], params['Acclimation_rate_eco']]
    deacclimation_rate = [None, params['Deacclimation_rate_endo'], params['Deacclimation_rate_eco']]
    theta = [None, 1, params['Theta']]
    
    # Initialize variables
    model_Hc_yesterday = params['Hc_initial']
    dormancy_period = 1
    DD_heating_sum = 0.0
    DD_chilling_sum = 0.0
    base10_chilling_sum = 0.0
    ecodormancy_boundary = params['Ecodormancy_boundary']
    budbreak_date = None
    
    if verbose:
        print(f"Running model for {variety}")
    
    # Main model loop
    for i, row in results.iterrows():
        T_mean = row['T_mean']
        acc_jday = row['acc_jday']
        
        if pd.isna(T_mean):
            continue
        
        # Calculate degree days
        if T_mean > T_threshold[dormancy_period]:
            DD_heating_today = T_mean - T_threshold[dormancy_period]
            DD_chilling_today = 0.0
        else:
            DD_heating_today = 0.0
            DD_chilling_today = T_mean - T_threshold[dormancy_period]
        
        # Base-10 chilling
        base10_chilling_today = T_mean - 10 if T_mean <= 10 else 0.0
        
        # Calculate acclimation and deacclimation
        deacclimation = (DD_heating_today * deacclimation_rate[dormancy_period] * 
                        (1 - ((model_Hc_yesterday - Hc_max) / Hc_range) ** theta[dormancy_period]))
        
        if DD_chilling_sum == 0:
            deacclimation = 0.0
        
        acclimation = (DD_chilling_today * acclimation_rate[dormancy_period] * 
                      (1 - (Hc_min - model_Hc_yesterday) / Hc_range))
        
        # Update model
        model_Hc = model_Hc_yesterday + acclimation + deacclimation
        model_Hc = max(Hc_max, min(Hc_min, model_Hc))
        
        # Update sums
        DD_chilling_sum += DD_chilling_today
        base10_chilling_sum += base10_chilling_today
        
        if dormancy_period == 2:
            DD_heating_sum += DD_heating_today
        
        # Check dormancy transition
        if base10_chilling_sum <= ecodormancy_boundary:
            dormancy_period = 2
        
        # Store results
        results.loc[i, 'Predicted_Hc'] = model_Hc
        results.loc[i, 'Dormancy_period'] = dormancy_period
        results.loc[i, 'DD_heating_sum'] = DD_heating_sum
        results.loc[i, 'DD_chilling_sum'] = DD_chilling_sum
        results.loc[i, 'Base10_chilling_sum'] = base10_chilling_sum
        results.loc[i, 'ACC_jday'] = acc_jday
        
        model_Hc_yesterday = model_Hc
    
    results['Budbreak_date'] = budbreak_date
    return results

def parse_filename(filename: str) -> dict:
    """Parse filename and extract treatment information."""
    basename = os.path.basename(filename).replace('_complete.csv', '')
    parts = basename.split('_')
    
    if len(parts) >= 4:
        if len(parts) >= 5:
            treatment_parts = parts[:-2]
            original_treatment = '_'.join(treatment_parts)
            replicate = int(parts[-1])
        else:
            original_treatment = f"{parts[0]}_{parts[1]}"
            replicate = int(parts[3])
    else:
        original_treatment = basename
        replicate = 1
    
    # Treatment name mapping
    treatment_name_mapping = {
        'Paper_low': 'Paper Low',
        'Paper_Low': 'Paper Low',
        'Paper_raise': 'Paper High', 
        'Paper_Raise': 'Paper High',
        'Plastic_low': 'Plastic Low',
        'Plastic_Low': 'Plastic Low',
        'Plastic_raise': 'Plastic High',
        'Plastic_Raise': 'Plastic High',
        'Uncover_uncover': 'No-Tube',
        'No_Tube': 'No-Tube',
    }
    
    mapped_treatment = treatment_name_mapping.get(original_treatment, original_treatment)
    
    return {
        'filename': filename,
        'replicate': replicate,
        'treatment': mapped_treatment,  
        'original_treatment': original_treatment, 
        'basename': basename
    }

# ============================================================================    
# VISUALIZATION FUNCTIONS
# ============================================================================

def plot_variety_treatment_comparison_with_weather(combined_data: pd.DataFrame, variety: str, 
                                                  output_dir: str, weather_file: str = 'weather_2024_2025_binned.csv', 
                                                  save_plots: bool = True):
    """Create variety-treatment comparison plot with ambient weather lines."""
    if combined_data is None or len(combined_data) == 0:
        print("No data available for variety treatment comparison plot")
        return None
    
    # Treatment styles
    treatment_styles = {
        'Paper Low': {'color': '#d62728', 'marker': 'o'},
        'Paper High': {'color': '#ff7f0e', 'marker': 's'},
        'Plastic Low': {'color': '#1f77b4', 'marker': 'v'},
        'Plastic High': {'color': '#9467bd', 'marker': 'D'},
        'No-Tube': {'color': '#2ca02c', 'marker': '^'},
    }
    
    # Load ambient weather data
    try:
        weather_data = pd.read_csv(weather_file)
        weather_data['datetime'] = pd.to_datetime(weather_data['timestamp'])
        weather_data['jday'] = weather_data['datetime'].dt.dayofyear
        
        # Handle year transition
        weather_data = weather_data.sort_values('datetime').reset_index(drop=True)
        weather_data['acc_jday'] = weather_data['jday'].copy()
        
        for i in range(1, len(weather_data)):
            prev_jday = weather_data.loc[i-1, 'jday']
            curr_jday = weather_data.loc[i, 'jday']
            if prev_jday > 300 and curr_jday < 100:
                weather_data.loc[i:, 'acc_jday'] = weather_data.loc[i:, 'jday'] + 365
        
        daily_weather = weather_data.groupby('acc_jday')['air_temp_C'].agg(['min', 'max']).reset_index()
        daily_weather.columns = ['acc_jday', 'ambient_tmin', 'ambient_tmax']
        daily_weather = daily_weather[
            (daily_weather['acc_jday'] >= 250) & (daily_weather['acc_jday'] <= 455)
        ].copy()
        
    except Exception as e:
        print(f"Warning: Could not load weather data - {e}")
        daily_weather = None
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    treatments = sorted(combined_data['Treatment'].unique())
    plotted_treatments = set()
    
    for treatment in treatments:
        treatment_data = combined_data[combined_data['Treatment'] == treatment]
        
        if len(treatment_data) == 0:
            continue
            
        treatment_stats = treatment_data.groupby('acc_jday')['Predicted_Hc'].mean().reset_index()
        treatment_stats = treatment_stats.sort_values('acc_jday')
        
        # Determine display name and style
        if "Paper" in treatment and ("raise" in treatment.lower() or "high" in treatment.lower()):
            display_name = "Paper High"
            style = {'color': '#ff7f0e', 'marker': 's'}
        elif "Paper" in treatment and "low" in treatment.lower():
            display_name = "Paper Low"
            style = {'color': '#d62728', 'marker': 'o'}
        elif "Plastic" in treatment and ("raise" in treatment.lower() or "high" in treatment.lower()):
            display_name = "Plastic High"
            style = {'color': '#9467bd', 'marker': 'D'}
        elif "Plastic" in treatment and "low" in treatment.lower():
            display_name = "Plastic Low"
            style = {'color': '#1f77b4', 'marker': 'v'}
        elif "No" in treatment or "tube" in treatment.lower() or "Uncover" in treatment:
            display_name = "No-Tube"
            style = {'color': '#2ca02c', 'marker': '^'}
        else:
            display_name = treatment.replace("Raise", "High").replace("raise", "high").replace("_", " ")
            style = treatment_styles.get(display_name, {'color': '#000000', 'marker': 'o'})
        
        if display_name in plotted_treatments:
            continue
            
        plotted_treatments.add(display_name)
        
        ax.plot(treatment_stats['acc_jday'], treatment_stats['Predicted_Hc'], 
               color=style['color'], linewidth=3, label=display_name, 
               marker=style['marker'], markersize=4, markeredgewidth=0, alpha=0.9)
    
    # Plot ambient temperature lines
    if daily_weather is not None and len(daily_weather) > 0:
        ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmax'], 
               color='red', linestyle='--', linewidth=2, alpha=0.7, label='Ambient Tmax')
        ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmin'], 
               color='blue', linestyle='--', linewidth=2, alpha=0.7, label='Ambient Tmin')
    
    # Formatting
    ax.set_xlabel('Julian Day', fontsize=30)
    ax.set_ylabel('Temperature (°C)', fontsize=30)
    ax.set_title(f'{variety}', fontsize=33)
    ax.set_xlim(240, 455)
    ax.set_xticks([250, 275, 300, 325, 350, 365, 390, 415, 440, 455])
    ax.set_xticklabels(['250', '275', '300', '325', '350', '365', '25', '50', '75', '90'])
    ax.set_ylim(-40, 40)
    ax.tick_params(axis='both', labelsize=28)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', ncol=3, fontsize=21, frameon=False)
    
    plt.tight_layout()
    
    if save_plots:
        plot_file = os.path.join(output_dir, f"{variety}_variety_treatment_comparison.png")
        plt.savefig(plot_file, dpi=400, bbox_inches='tight', facecolor='white')
        print(f"Variety treatment comparison plot saved: {plot_file}")
    
    return fig

def plot_varieties_by_treatment(varieties_data: dict, output_dir: str, season: str, 
                               weather_file: str = 'weather_2024_2025_binned.csv', 
                               save_plots: bool = True) -> list:
    """Create separate plots for each treatment showing all varieties."""
    variety_colors = {
        'Chardonnay': '#1f77b4',
        'Concord': '#ff7f0e',
        'Cabernet Sauvignon': '#2ca02c',
        'Mourvedre': '#d62728',
        'Barbera': '#9467bd',
        'Cabernet franc': '#8c564b',
        'Chenin blanc': '#e377c2',
        'Dolcetto': '#7f7f7f',
        'Gewürztraminer': '#bcbd22',
        'Grenache': '#17becf'
    }
    
    # Get treatments
    all_treatments = set()
    for variety_name, combined_data in varieties_data.items():
        if combined_data is not None and len(combined_data) > 0:
            all_treatments.update(combined_data['Treatment'].unique())
    
    treatments = sorted(list(all_treatments))
    
    # Load weather data (same as above function)
    try:
        weather_data = pd.read_csv(weather_file)
        weather_data['datetime'] = pd.to_datetime(weather_data['timestamp'])
        weather_data['jday'] = weather_data['datetime'].dt.dayofyear
        
        weather_data = weather_data.sort_values('datetime').reset_index(drop=True)
        weather_data['acc_jday'] = weather_data['jday'].copy()
        
        for i in range(1, len(weather_data)):
            prev_jday = weather_data.loc[i-1, 'jday']
            curr_jday = weather_data.loc[i, 'jday']
            if prev_jday > 300 and curr_jday < 100:
                weather_data.loc[i:, 'acc_jday'] = weather_data.loc[i:, 'jday'] + 365
        
        daily_weather = weather_data.groupby('acc_jday')['air_temp_C'].agg(['min', 'max']).reset_index()
        daily_weather.columns = ['acc_jday', 'ambient_tmin', 'ambient_tmax']
        daily_weather = daily_weather[
            (daily_weather['acc_jday'] >= 250) & (daily_weather['acc_jday'] <= 455)
        ].copy()
        
    except Exception as e:
        print(f"Warning: Could not load weather data - {e}")
        daily_weather = None
    
    created_plots = []
    
    for treatment in treatments:
        fig, ax = plt.subplots(figsize=(16, 10), dpi=150)
        varieties_plotted = 0
        
        for variety_name, combined_data in varieties_data.items():
            if combined_data is None or len(combined_data) == 0:
                continue
                
            treatment_data = combined_data[combined_data['Treatment'] == treatment]
            if len(treatment_data) == 0:
                continue
            
            variety_treatment_mean = treatment_data.groupby('acc_jday')['Predicted_Hc'].mean().reset_index()
            variety_treatment_mean = variety_treatment_mean.sort_values('acc_jday')
            
            color = variety_colors.get(variety_name, '#000000')
            
            ax.plot(variety_treatment_mean['acc_jday'], variety_treatment_mean['Predicted_Hc'], 
                   color=color, linewidth=3, label=variety_name, 
                   marker='o', markersize=1, markeredgewidth=0, alpha=0.9)
            
            varieties_plotted += 1
        
        if varieties_plotted == 0:
            plt.close(fig)
            continue
        
        # Plot ambient temperature
        if daily_weather is not None and len(daily_weather) > 0:
            ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmax'], 
                   color='red', linestyle='--', linewidth=2, alpha=0.7, label='Ambient Tmax')
            ax.plot(daily_weather['acc_jday'], daily_weather['ambient_tmin'], 
                   color='blue', linestyle='--', linewidth=2, alpha=0.7, label='Ambient Tmin')
        
        # Formatting
        ax.set_xlabel('Julian Day', fontsize=35)
        ax.set_ylabel('Temperature (°C)', fontsize=35)
        ax.set_title(f'{treatment}', fontsize=40, pad=20)
        ax.set_xlim(240, 460)
        ax.set_xticks([250, 275, 300, 325, 350, 365, 390, 415, 440, 455])
        ax.set_xticklabels(['250', '275', '300', '325', '350', '365', '25', '50', '75', '90'])
        ax.set_ylim(-40, 40)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=25, ncol=3, frameon=False)
        ax.tick_params(axis='both', which='major', labelsize=32)
        
        plt.tight_layout()
        
        if save_plots:
            plot_file = os.path.join(output_dir, f"varieties_comparison_{treatment.replace(' ', '_').replace('-', '_')}_{season}.png")
            plt.savefig(plot_file, dpi=400, bbox_inches='tight', facecolor='white')
            created_plots.append(plot_file)
        
        plt.close(fig)
    
    print(f"Created {len(created_plots)} treatment comparison plots")
    return created_plots

# ============================================================================
# MAIN PROCESSOR CLASS
# ============================================================================

class GrowTubeColdHardnessProcessor:
    """Process grow tube experiment files for cold hardiness analysis."""
    
    def __init__(self):
        """Initialize with model parameters for grape varieties."""
        
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

        print(f"Cold Hardiness Processor Initialized")
        print(f"Available varieties: {len(self.variety_params)}")
    
    def find_julian_day_250(self, data: pd.DataFrame, verbose: bool = False) -> int:
        """Find Julian Day 250 in the data."""
        data['datetime'] = pd.to_datetime(data['datetime'])
        data['jday'] = data['datetime'].dt.dayofyear
        
        target_jday = 250
        available_jdays = data['jday'].unique()
        closest_jday = min(available_jdays, key=lambda x: abs(x - target_jday))
        
        if verbose:
            print(f"Using Julian Day: {closest_jday}")
        
        return closest_jday
    
    def process_single_file(self, filepath: str, variety: str, verbose: bool = True) -> dict:
        """Process a single file and return results."""
        filename = os.path.basename(filepath)
        file_info = parse_filename(filename)
        
        if verbose:
            print(f"Processing: {filename} for {variety}")
        
        try:
            data = pd.read_csv(filepath)
            start_jday = self.find_julian_day_250(data, verbose=False)
            daily_data = aggregate_to_daily(data, start_jday, verbose=False)
            results = apply_model(daily_data, variety, file_info, self.variety_params, verbose=False)
            
            return {
                'success': True,
                'filepath': filepath,
                'results': results,
                'variety': variety
            }
            
        except Exception as e:
            if verbose:
                print(f"Error processing {filename}: {e}")
            return {
                'success': False,
                'filepath': filepath,
                'error': str(e),
                'variety': variety
            }
    
    def process_variety(self, input_dir: str, variety: str, verbose: bool = True) -> pd.DataFrame:
        """Process all files for a variety and return combined data."""
        if verbose:
            print(f"Processing variety: {variety}")
        
        file_pattern_path = os.path.join(input_dir, "*_complete.csv")
        files = glob.glob(file_pattern_path)
        
        if len(files) == 0:
            print(f"No files found for {variety}")
            return None
        
        all_results = []
        successful_files = 0
        
        for filepath in files:
            result = self.process_single_file(filepath, variety, verbose=False)
            
            if result['success']:
                all_results.append(result['results'])
                successful_files += 1
        
        if len(all_results) == 0:
            print(f"No successful processing for {variety}")
            return None
        
        combined_data = pd.concat(all_results, ignore_index=True)
        
        if verbose:
            print(f"Successfully processed {successful_files} files")
            print(f"Combined data: {len(combined_data)} daily observations")
            print(f"Treatments: {combined_data['Treatment'].unique()}")
        
        return combined_data
    
    def save_combined_csv(self, combined_data: pd.DataFrame, variety: str, output_dir: str) -> None:
        """Save combined CSV file for a variety."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save detailed results
        combined_file = os.path.join(output_dir, f"combined_detailed_results_{variety.replace(' ', '_')}.csv")
        combined_data.to_csv(combined_file, index=False)
        print(f"Combined CSV saved: {combined_file}")
        
        # Create summary
        summary_data = []
        for treatment in combined_data['Treatment'].unique():
            treatment_data = combined_data[combined_data['Treatment'] == treatment]
            for replicate in treatment_data['Replicate'].unique():
                rep_data = treatment_data[treatment_data['Replicate'] == replicate]
                valid_hc = rep_data['Predicted_Hc'].dropna()
                
                if len(valid_hc) > 0:
                    summary_data.append({
                        'Variety': variety,
                        'Treatment': treatment,
                        'Replicate': replicate,
                        'Final_Hc': valid_hc.iloc[-1],
                        'Min_Hc': valid_hc.min(),
                        'Max_Hc': valid_hc.max(),
                        'N_Days': len(valid_hc),
                        'Start_Date': rep_data['datetime'].min(),
                        'End_Date': rep_data['datetime'].max()
                    })
        
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_file = os.path.join(output_dir, f"combined_summary_{variety.replace(' ', '_')}.csv")
            summary_df.to_csv(summary_file, index=False)
            print(f"Summary CSV saved: {summary_file}")

# ============================================================================
# MAIN WORKFLOW FUNCTION
# ============================================================================

def run_streamlined_analysis(weather_file: str, treatment_file: str, season: str, 
                            varieties: list, output_base_dir: str) -> dict:
    """
    Streamlined workflow: data preparation + analysis + comparison plots.
    """
    print(f"STREAMLINED COLD HARDINESS ANALYSIS FOR {season}")
    print(f"Varieties: {varieties}")
    print("="*80)
    
    # Step 1: Prepare data
    prepared_data_dir = os.path.join(output_base_dir, f"{season}_prepared_data")
    
    if not os.path.exists(prepared_data_dir):
        print(f"Step 1: Preparing data for {season}...")
        
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
        
        print(f"Created {len(created_files)} data files")
    else:
        print(f"Using existing prepared data in {prepared_data_dir}")
    
    # Step 2: Process varieties
    print(f"Step 2: Processing varieties...")
    
    processor = GrowTubeColdHardnessProcessor()
    varieties_data = {}
    
    for variety in varieties:
        print(f"{'='*60}")
        print(f"PROCESSING: {variety}")
        print(f"{'='*60}")
        
        combined_data = processor.process_variety(prepared_data_dir, variety, verbose=True)
        
        if combined_data is not None:
            varieties_data[variety] = combined_data
            
            variety_output_dir = os.path.join(output_base_dir, f"{season}_{variety.replace(' ', '_')}")
            processor.save_combined_csv(combined_data, variety, variety_output_dir)
            
            try:
                plot_variety_treatment_comparison_with_weather(
                    combined_data, variety, variety_output_dir, 
                    weather_file=weather_file, save_plots=True
                )
                print(f"Individual comparison plot created for {variety}")
            except Exception as plot_error:
                print(f"Warning: Could not create individual plot for {variety} - {plot_error}")
        else:
            print(f"No data processed for {variety}")
            varieties_data[variety] = None
    
    # Step 3: Create treatment-based comparison plots
    print(f"{'='*60}")
    print("CREATING TREATMENT COMPARISON PLOTS")
    print(f"{'='*60}")
    
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
            print(f"Treatment comparison plots created: {len(created_plots)}")
            print(f"Varieties included: {list(valid_varieties_data.keys())}")
        except Exception as e:
            print(f"Error creating treatment comparison plots: {e}")
    else:
        print("No valid variety data available for treatment comparison plots")
    
    print(f"STREAMLINED ANALYSIS COMPLETE!")
    print(f"Results saved in: {output_base_dir}")
    
    return varieties_data

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    print("STREAMLINED GROW TUBE COLD HARDINESS ANALYSIS")
    print("Comparison plots and combined CSV files")
    print("="*80)
    
    # Configuration
    TREATMENT_FILE = "grow_tube_categorized.csv"
    SEASON = "2023_2024" # Chage the season here!
    
    if SEASON == "2023_2024":
        WEATHER_FILE = "weather_2023_2024_binned.csv"
    elif SEASON == "2024_2025":
        WEATHER_FILE = "weather_2024_2025_binned.csv"
    else:
        print(f"Unknown season: {SEASON}")
        return
    
    VARIETIES = ['Chardonnay', 'Concord', 'Cabernet Sauvignon', 'Mourvedre']
    OUTPUT_DIR = f"cold_hardiness_results_{SEASON}"
    
    # File validation
    print(f"Season: {SEASON}")
    print(f"Weather file: {WEATHER_FILE}")
    print(f"Treatment file: {TREATMENT_FILE}")
    print(f"Varieties: {VARIETIES}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("-" * 60)
    
    if not os.path.exists(WEATHER_FILE):
        print(f"Weather file not found: {WEATHER_FILE}")
        return
    else:
        print(f"Weather file found: {WEATHER_FILE}")
    
    if not os.path.exists(TREATMENT_FILE):
        print(f"Treatment file not found: {TREATMENT_FILE}")
        return
    else:
        print(f"Treatment file found: {TREATMENT_FILE}")
    
    # Run analysis
    print(f"Starting streamlined analysis for {SEASON}...")
    
    try:
        varieties_data = run_streamlined_analysis(
            weather_file=WEATHER_FILE,
            treatment_file=TREATMENT_FILE,
            season=SEASON,
            varieties=VARIETIES,
            output_base_dir=OUTPUT_DIR
        )
        
        print(f"ANALYSIS COMPLETE FOR {SEASON}!")
        print(f"Results saved in: {OUTPUT_DIR}")
        print("Files created:")
        print("- Combined CSV files for each variety")
        print("- Individual variety comparison plots")
        print("- Treatment-based comparison plots")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        print("Please check:")
        print("1. File paths are correct")
        print("2. Data files have the expected format")
        print("3. The season exists in your data")

if __name__ == "__main__":
    main()