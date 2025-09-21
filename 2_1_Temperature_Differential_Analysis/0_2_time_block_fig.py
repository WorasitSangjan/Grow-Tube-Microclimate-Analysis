# Author: Worasit Sangjan
# Date: May 7, 2025

# ============================================================================
# TIME BLOCK ANALYSIS - 4-PERIOD TIME BLOCK DETERMINATION
# Analyzes weather data to establish time blocks for correction factors
# Creates comprehensive visualization of diurnal patterns
# ============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

# Set plotting style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['lines.linewidth'] = 1.0

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def load_weather_datasets() -> pd.DataFrame:
    """Load both weather datasets and combine for time block analysis."""
    print("Loading weather datasets...")
    
    df_2023_24 = pd.read_csv("cleaned_weather_2023_2024.csv")
    df_2024_25 = pd.read_csv("cleaned_weather_2024_2025.csv")
    
    # Standardize column names
    if 'solar_radiation_Wpm2' in df_2023_24.columns:
        df_2023_24.rename(columns={'solar_radiation_Wpm2': 'solar_radiation'}, inplace=True)
    if 'solar_radiation_Wpm2' in df_2024_25.columns:
        df_2024_25.rename(columns={'solar_radiation_Wpm2': 'solar_radiation'}, inplace=True)
    
    df_2023_24['season'] = '2023-2024'
    df_2024_25['season'] = '2024-2025'
    
    merged_df = pd.concat([df_2023_24, df_2024_25], ignore_index=True)
    print(f"Combined dataset: {len(merged_df):,} records")
    
    return merged_df

def filter_dormancy_period(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Filter data for dormancy period: Oct 25 - Mar 31."""
    print("Filtering for dormancy period (Oct 25 - Mar 31)...")
    
    merged_df["timestamp"] = pd.to_datetime(merged_df["timestamp"], errors="coerce")
    merged_df["hour"] = merged_df["timestamp"].dt.hour
    merged_df["month"] = merged_df["timestamp"].dt.month
    merged_df["day"] = merged_df["timestamp"].dt.day
    
    def is_dormancy_period(row):
        month, day = row['month'], row['day']
        return ((month == 10 and day >= 25) or 
                month in [11, 12, 1, 2] or 
                (month == 3 and day <= 31))
    
    dormancy_data = merged_df[merged_df.apply(is_dormancy_period, axis=1)].copy()
    print(f"Dormancy period data: {len(dormancy_data):,} records")
    
    return dormancy_data

def assign_time_blocks(weather_df: pd.DataFrame) -> pd.DataFrame:
    """Assign time blocks based on hour of day."""
    print("Assigning time blocks...")
    
    def get_time_block(hour):
        if 7 <= hour < 11:
            return "Morning warming"
        elif 11 <= hour < 15:
            return "Peak heat"
        elif 15 <= hour < 19:
            return "Afternoon cooling"
        else:
            return "Night period"
    
    weather_df["time_block"] = weather_df["hour"].apply(get_time_block)
    print(f"Time blocks assigned to {len(weather_df):,} records")
    
    return weather_df

def calculate_hourly_summaries(weather_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate hourly summaries for plotting."""
    print("Calculating hourly summaries...")
    
    hourly_stats = (
        weather_df.groupby("hour")[["air_temp_C", "solar_radiation"]]
        .agg(["mean", "min", "max"])
        .round(2)
    )
    
    # Flatten column names
    hourly_stats.columns = ['_'.join(col).strip() for col in hourly_stats.columns.values]
    hourly_stats = hourly_stats.reset_index()
    
    return hourly_stats

def create_time_block_figure(hourly_stats: pd.DataFrame, weather_df: pd.DataFrame) -> None:
    """Create comprehensive time block analysis figure."""
    print("Creating time block figure...")
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7.25, 9.5))
    
    # Prepare hourly data for circular plot (add hour 24 = hour 0)
    hourly_plot = hourly_stats.copy()
    hour_0_row = hourly_stats[hourly_stats['hour'] == 0].copy()
    hour_0_row['hour'] = 24
    hourly_plot = pd.concat([hourly_plot, hour_0_row], ignore_index=True)
    
    # GRAPH 1: DIURNAL PATTERNS WITH TIME BLOCKS
    # Temperature lines
    line1 = ax1.plot(hourly_plot["hour"], hourly_plot["air_temp_C_mean"], 
                     label="Mean temp", color="steelblue", linewidth=1.0, linestyle='-')
    line2 = ax1.plot(hourly_plot["hour"], hourly_plot["air_temp_C_min"], 
                     label="Min temp", color="steelblue", linewidth=1.0, linestyle='--')
    line3 = ax1.plot(hourly_plot["hour"], hourly_plot["air_temp_C_max"], 
                     label="Max temp", color="steelblue", linewidth=1.0, linestyle=':')
    
    ax1.set_ylabel("Air temperature (°C)", color="steelblue")
    ax1.set_xlabel("Hour of day")
    ax1.set_xticks(range(0, 25, 2))
    ax1.set_xlim(0, 24)
    ax1.tick_params(axis='y', labelcolor="steelblue", width=0.5)
    ax1.tick_params(axis='x', width=0.5)
    ax1.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)
    
    # Solar radiation on secondary axis
    ax1_twin = ax1.twinx()
    line4 = ax1_twin.plot(hourly_plot["hour"], hourly_plot["solar_radiation_mean"], 
                          label="Mean radiation", color="darkorange", linewidth=1.0, linestyle='-')
    line5 = ax1_twin.plot(hourly_plot["hour"], hourly_plot["solar_radiation_min"], 
                          label="Min radiation", color="darkorange", linewidth=1.0, linestyle='--')
    line6 = ax1_twin.plot(hourly_plot["hour"], hourly_plot["solar_radiation_max"], 
                          label="Max radiation", color="darkorange", linewidth=1.0, linestyle=':')
    
    ax1_twin.set_ylabel("Solar radiation (W/m²)", color="darkorange")
    ax1_twin.tick_params(axis='y', labelcolor="darkorange", width=0.5, direction='out')
    
    # Add time block boundaries
    time_dividers = [7, 11, 15, 19]
    for divider in time_dividers:
        ax1.axvline(x=divider, color="gray", linestyle="--", linewidth=0.5, alpha=0.6)
    
    # Time block labels
    time_labels = [
        ("Morning warming | 07-11", 7),
        ("Peak heat | 11-15", 11),
        ("Afternoon cooling | 15-19", 15),
        ("Night period | 19-07", 19)
    ]
    for label, hour in time_labels:
        ax1.text(hour + 0.25, 17, label, rotation=90, ha="left", va="top",
                 fontsize=8, color="gray", alpha=0.7)
    
    # Legend
    lines = line1 + line2 + line3 + line4 + line5 + line6
    labels = ["Mean temp", "Min temp", "Max temp", "Mean radiation", "Min radiation", "Max radiation"]
    ax1.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.15), 
               ncol=6, frameon=False, fontsize=8)
    
    # GRAPH 2: SOLAR RADIATION DISTRIBUTION BY TIME BLOCK
    time_block_order = ["Morning warming", "Peak heat", "Afternoon cooling", "Night period"]
    colors = ['white', '#CCCCCC', '#808080', '#404040']
    
    solar_violin_data = []
    for block in time_block_order:
        block_solar = weather_df[weather_df["time_block"] == block]["solar_radiation"].dropna()
        solar_violin_data.append(block_solar)
    
    solar_violin_parts = ax2.violinplot(solar_violin_data, positions=range(len(time_block_order)), 
                                       showmeans=True, showmedians=True)
    
    for i, pc in enumerate(solar_violin_parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(1.0)
        if i == 0:
            pc.set_edgecolor('black')
            pc.set_linewidth(0.25)
    
    overall_solar_mean = weather_df["solar_radiation"].mean()
    ax2.axhline(overall_solar_mean, color='red', linestyle='-', linewidth=1.0,
                label=f"Overall mean: {overall_solar_mean:.0f} W/m²")
    
    ax2.set_xlabel("Time block")
    ax2.set_ylabel("Solar radiation (W/m²)")
    ax2.set_xticks(range(len(time_block_order)))
    ax2.set_xticklabels([name.replace(' ', '\n') for name in time_block_order])
    ax2.tick_params(axis='x', which='both', length=0)
    ax2.tick_params(axis='y', width=0.5, direction='out')
    
    legend2 = ax2.legend(loc='upper right', frameon=True, fontsize=8)
    legend2.get_frame().set_linewidth(0.5)
    
    ax2.grid(True, alpha=0.3, linewidth=0.5)
    
    # GRAPH 3: TEMPERATURE DISTRIBUTION BY TIME BLOCK
    temp_violin_data = []
    for block in time_block_order:
        block_temps = weather_df[weather_df["time_block"] == block]["air_temp_C"].dropna()
        temp_violin_data.append(block_temps)
    
    temp_violin_parts = ax3.violinplot(temp_violin_data, positions=range(len(time_block_order)), 
                                      showmeans=True, showmedians=True)
    
    for i, pc in enumerate(temp_violin_parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(1.0)
        if i == 0:
            pc.set_edgecolor('black')
            pc.set_linewidth(0.25)
    
    overall_temp_mean = weather_df["air_temp_C"].mean()
    ax3.axhline(overall_temp_mean, color='red', linestyle='-', linewidth=1.0, 
                label=f"Overall mean: {overall_temp_mean:.1f}°C")
    
    ax3.set_xlabel("Time block")
    ax3.set_ylabel("Air temperature (°C)")
    ax3.set_xticks(range(len(time_block_order)))
    ax3.set_xticklabels([name.replace(' ', '\n') for name in time_block_order])
    ax3.tick_params(axis='x', which='both', length=0)
    ax3.tick_params(axis='y', width=0.5, direction='out')
    
    legend3 = ax3.legend(loc='upper right', frameon=True, fontsize=8)
    legend3.get_frame().set_linewidth(0.5)
    
    ax3.grid(True, alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4)
    
    plt.savefig('time_block_analysis.png', dpi=600, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.show()

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def time_block_analysis_workflow() -> None:
    """Main workflow to create time block analysis figure."""
    print("="*60)
    print("TIME BLOCK ANALYSIS")
    print("="*60)
    
    # Load and process data
    merged_df = load_weather_datasets()
    dormancy_data = filter_dormancy_period(merged_df)
    weather_df = assign_time_blocks(dormancy_data)
    
    # Calculate summaries
    hourly_stats = calculate_hourly_summaries(weather_df)
    
    # Create figure
    create_time_block_figure(hourly_stats, weather_df)
    
    print("\nAnalysis complete!")

# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    time_block_analysis_workflow()