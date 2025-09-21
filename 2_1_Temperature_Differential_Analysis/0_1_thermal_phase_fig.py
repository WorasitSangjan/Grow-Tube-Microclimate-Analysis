# Author: Worasit Sangjan
# Date: May 7, 2025

# ============================================================================
# THERMAL PHASE ANALYSIS - DATA-DRIVEN THRESHOLD DETERMINATION
# Analyzes weather data for dormancy period (October 25 - March 30)
# Determines data-driven thresholds for temperature and solar radiation
# ============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# Set plotting style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.linewidth'] = 1.25
plt.rcParams['lines.linewidth'] = 1.0

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def load_weather_datasets() -> pd.DataFrame:
    """Load both weather datasets and combine for analysis."""
    print("Loading weather datasets...")
    
    df_2023_24 = pd.read_csv("cleaned_weather_2023_2024.csv")
    df_2024_25 = pd.read_csv("cleaned_weather_2024_2025.csv")
    
    df_2023_24['season'] = '2023-2024'
    df_2024_25['season'] = '2024-2025'
    
    merged_df = pd.concat([df_2023_24, df_2024_25], ignore_index=True)
    print(f"Combined dataset: {len(merged_df):,} records")
    
    return merged_df

def filter_dormancy_period(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Filter data for dormancy period: Oct 25 - Mar 30."""
    print("Filtering for dormancy period (Oct 25 - Mar 30)...")
    
    merged_df["timestamp"] = pd.to_datetime(merged_df["timestamp"], errors="coerce")
    merged_df["month"] = merged_df["timestamp"].dt.month
    merged_df["day"] = merged_df["timestamp"].dt.day
    merged_df["date"] = merged_df["timestamp"].dt.date
    
    def is_dormancy_period(row):
        month, day = row['month'], row['day']
        return ((month == 10 and day >= 25) or 
                month in [11, 12, 1, 2] or 
                (month == 3 and day <= 30))
    
    dormancy_data = merged_df[merged_df.apply(is_dormancy_period, axis=1)].copy()
    print(f"Dormancy period data: {len(dormancy_data):,} records")
    
    return dormancy_data

def calculate_daily_summaries(dormancy_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate daily weather summaries."""
    print("Calculating daily summaries...")
    
    daily_summary = (
        dormancy_data.groupby(["date", "season"])
        .agg({
            "air_temp_C": ["min", "mean", "max"],
            "solar_radiation_Wpm2": "mean"
        })
        .reset_index()
    )
    
    # Flatten column names
    daily_summary.columns = ['date', 'season', 'air_temp_min', 'air_temp_mean', 'air_temp_max', 'solar_radiation_mean']
    daily_summary['diurnal_range'] = daily_summary['air_temp_max'] - daily_summary['air_temp_min']
    
    print(f"Daily summaries: {len(daily_summary)} days")
    return daily_summary

def determine_thresholds(daily_summary: pd.DataFrame) -> dict:
    """Determine data-driven thresholds."""
    print("Determining thresholds...")
    
    # Temperature statistics
    temp_stats = daily_summary['air_temp_mean'].describe()
    temp_mean, temp_std = temp_stats['mean'], temp_stats['std']
    temp_q1, temp_q3 = temp_stats['25%'], temp_stats['75%']
    
    # Solar statistics
    solar_stats = daily_summary['solar_radiation_mean'].describe()
    solar_q1, solar_mean, solar_q3 = solar_stats['25%'], solar_stats['mean'], solar_stats['75%']
    
    thresholds = {
        'temperature_bins': [
            (temp_mean - 2*temp_std, temp_mean - temp_std),
            (temp_mean - temp_std, temp_mean),
            (temp_mean, temp_mean + temp_std),
            (temp_mean + temp_std, temp_mean + 2*temp_std)
        ],
        'solar_thresholds': {
            'solar_start': solar_q1,
            'solar_moderate': solar_mean,
            'solar_peak': solar_q3
        }
    }
    
    print("Temperature bins (°C):")
    for i, (low, high) in enumerate(thresholds['temperature_bins']):
        print(f"  Bin {i+1}: {low:.1f} to {high:.1f}")
    
    print(f"\nSolar thresholds (W/m²):")
    print(f"  Start: {solar_q1:.1f}, Peak: {solar_q3:.1f}")
    
    return thresholds

def create_distribution_plots(daily_summary: pd.DataFrame, thresholds: dict) -> None:
    """Create comprehensive distribution plots with data-driven thresholds."""
    print("Creating plots...")

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.2])

    # Define subplot axes
    box_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]  # Top row: boxplots
    hist_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]  # Bottom row: histograms

    # Data and styling - restore original colors with accessibility features
    colors = ['steelblue', 'darkorange', 'seagreen']  # Original colors
    variables = ['air_temp_mean', 'diurnal_range', 'solar_radiation_mean']
    titles = ['Daily mean air temperature', 'Daily diurnal range', 'Daily mean solar radiation']
    units = ['°C', '°C', 'W/m²']

    for i, (var, title, unit, color) in enumerate(zip(variables, titles, units, colors)):
        stats = daily_summary[var].describe()

        # --- Histogram (bottom row) ---
        ax_left = hist_axes[i]
        ax_left.hist(daily_summary[var], bins=25, color=color,
                     alpha=0.7, edgecolor='black', linewidth=0.5)

        # Add statistics lines with different styles for accessibility
        ax_left.axvline(stats['mean'], color='red', linestyle='-', linewidth=2.0,
                        label=f"Mean: {stats['mean']:.1f} {unit}")
        ax_left.axvline(stats['50%'], color='green', linestyle='--', linewidth=2.0,
                        label=f"Median: {stats['50%']:.1f} {unit}")
        ax_left.axvline(stats['25%'], color='orange', linestyle=':', linewidth=2.0,
                        label=f"Q1: {stats['25%']:.1f} {unit}")
        ax_left.axvline(stats['75%'], color='purple', linestyle=':', linewidth=2.0,
                        label=f"Q3: {stats['75%']:.1f} {unit}")
        ax_left.axvspan(stats['mean'] - stats['std'], stats['mean'] + stats['std'],
                        alpha=0.2, color='gray', label='±1 std dev')

        # Threshold markers
        if var == 'air_temp_mean':
            ax_left.axvline(stats['25%'], color='blue', linestyle='-.', linewidth=2.0,
                            label=f"Cold dormancy: {stats['25%']:.1f}°C")
            ax_left.axvline(stats['75%'], color='red', linestyle='-.', linewidth=2.0,
                            label=f"Warm dormancy: {stats['75%']:.1f}°C")

        elif var == 'solar_radiation_mean':
            solar_thresh = thresholds['solar_thresholds']
            ax_left.axvline(solar_thresh['solar_start'], color='orange', linestyle='-.', linewidth=2.0,
                            label=f"Low solar: {solar_thresh['solar_start']:.1f} W/m²")
            ax_left.axvline(solar_thresh['solar_peak'], color='red', linestyle='-.', linewidth=2.0,
                            label=f"High solar: {solar_thresh['solar_peak']:.1f} W/m²")

        ax_left.set_xlabel(f'{title} ({unit})', fontsize=18)
        ax_left.set_ylabel('Count', fontsize=18)
        ax_left.tick_params(axis='both', width=0.5, labelsize=17) 
        
        # Legend below the graph
        ax_left.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), 
                      ncol=2, frameon=False, fontsize=16)
        
        ax_left.grid(True, alpha=0.3, linewidth=0.5)

        # --- Boxplot (top row) ---
        ax_right = box_axes[i]
        season_data = [daily_summary[daily_summary['season'] == season][var]
                       for season in ['2023-2024', '2024-2025']]

        box_plot = ax_right.boxplot(season_data, labels=['2023-2024', '2024-2025'],
                                    patch_artist=True, notch=True)

        # Use different colors for seasons with visual distinction
        for patch, season_color in zip(box_plot['boxes'], ['lightblue', 'lightcoral']):
            patch.set_facecolor(season_color)
            patch.set_alpha(0.7)

        ax_right.axhline(stats['mean'], color='red', linestyle='-', linewidth=1.0,
                         label=f"Overall mean: {stats['mean']:.1f} {unit}")

        ax_right.set_ylabel(f'{title} ({unit})', fontsize=18)
        ax_right.tick_params(axis='both', width=0.5, labelsize=17)
        ax_right.tick_params(axis='x', which='both', length=0)  # No x-axis ticks for categorical
        
        # Legend below the graph
        ax_right.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), 
                       frameon=False, fontsize=16)
        
        ax_right.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.215, wspace=0.25)
    
    plt.savefig('thermal_phase_analysis.png', dpi=600, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.show()

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def thermal_phase_analysis_workflow() -> dict:
    """Main workflow to determine thresholds and create plots."""
    print("="*60)
    print("THERMAL PHASE ANALYSIS")
    print("="*60)
    
    # Load and process data
    merged_df = load_weather_datasets()
    dormancy_data = filter_dormancy_period(merged_df)
    daily_summary = calculate_daily_summaries(dormancy_data)
    
    # Determine thresholds
    thresholds = determine_thresholds(daily_summary)
    
    # Create plots
    create_distribution_plots(daily_summary, thresholds)
    
    print("\nAnalysis complete!")
    return thresholds

# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    thresholds = thermal_phase_analysis_workflow()