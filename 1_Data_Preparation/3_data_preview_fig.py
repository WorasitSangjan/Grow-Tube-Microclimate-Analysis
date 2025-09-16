# Author: Worasit Sangjan
# Date: May 7, 2025

# ============================================================================
# DATA PREVIEW: Temperature Data Visualization
# Creates comprehensive plots of merged HOBO and weather data
# ============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def load_and_validate_data(file_path: str, start_date: str = "2024-10-25") -> pd.DataFrame:
    """Load and validate temperature data."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found")
    
    print(f"Loading data from {file_path}...")
    
    df = pd.read_csv(file_path)
    df["datetime"] = pd.to_datetime(df["timestamp"])
    df = df[df["datetime"] >= start_date]
    
    combinations = df.groupby(['Material', 'Condition', 'Replicate']).size()
    
    print(f"Date range: {df['datetime'].min()} to {df['datetime'].max()}")
    print(f"Total records: {len(df):,}")
    print(f"Unique combinations: {len(combinations)}")
    
    return df

def get_treatment_info(df: pd.DataFrame) -> dict:
    """Get treatment combination information."""
    title_mapping = {
        ('Paper', 'low'): 'Paper Low',
        ('Paper', 'raise'): 'Paper Raise', 
        ('Plastic', 'low'): 'Plastic Low',
        ('Plastic', 'raise'): 'Plastic Raise',
        ('Uncover', 'uncover'): 'No-Tube'
    }
    
    combinations = df[['Material', 'Condition']].drop_duplicates()
    return combinations, title_mapping

def create_monthly_plots(df: pd.DataFrame) -> None:
    """Create detailed monthly plots for each treatment combination."""
    print("\nCreating monthly temperature plots...")
    
    combinations, title_mapping = get_treatment_info(df)
    
    # Group by month
    df["month"] = df["datetime"].dt.to_period("M")
    months = sorted(df["month"].unique())
    
    if len(months) == 0:
        print("No data to plot")
        return
    
    colors = ['red', 'green', 'blue', 'purple']
    
    for _, combo in combinations.iterrows():
        material, condition = combo['Material'], combo['Condition']
        combo_data = df[(df['Material'] == material) & (df['Condition'] == condition)]
        
        if len(combo_data) == 0:
            continue
        
        # Create subplots
        fig, axes = plt.subplots(len(months), 1, figsize=(16, 4 * len(months)))
        if len(months) == 1:
            axes = [axes]
        
        # Plot each month
        for i, month in enumerate(months):
            month_data = combo_data[combo_data["month"] == month]
            ax = axes[i]
            
            if len(month_data) == 0:
                ax.text(0.5, 0.5, f'No data for {month}', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"{month}", fontsize=12, fontweight='bold')
                continue
            
            # Plot replicates
            replicates = sorted(month_data['Replicate'].unique())
            for j, rep in enumerate(replicates):
                rep_data = month_data[month_data['Replicate'] == rep]
                if len(rep_data) > 0:
                    ax.plot(rep_data["datetime"], rep_data["Value"], 
                           label=f"Rep {rep} Tube", 
                           color=colors[j % len(colors)], 
                           linewidth=1.2, alpha=0.8)
            
            # Plot air temperature
            ax.plot(month_data["datetime"], month_data["air_temp_C"], 
                   label="Air Temp", color="black", linewidth=1, linestyle='--')
            
            # Format subplot
            ax.set_title(f"{month} - {title_mapping.get((material, condition), f'{material} {condition}')}", 
                        fontsize=12, fontweight='bold')
            ax.set_ylabel("Temperature (°C)")
            ax.grid(True, alpha=0.3)
            
            # Format x-axis
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
            
            # Set limits
            ax.set_xlim(month_data["datetime"].min(), month_data["datetime"].max())
            
            # Legend on first subplot
            if i == 0:
                ax.legend(loc='upper right', fontsize=8)
        
        # Final formatting
        treatment_title = title_mapping.get((material, condition), f"{material} - {condition}")
        plt.suptitle(f"Temperature Data: {treatment_title}", 
                    fontsize=16, fontweight='bold', y=0.98)
        plt.xlabel("Date")
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        plt.show()

def create_summary_plot(df: pd.DataFrame) -> None:
    """Create summary plot with all treatment combinations."""
    print("\nCreating summary temperature plot...")
    
    combinations, title_mapping = get_treatment_info(df)
    colors = ['red', 'green', 'blue', 'purple']
    
    fig, axes = plt.subplots(len(combinations), 1, figsize=(16, 5 * len(combinations)))
    if len(combinations) == 1:
        axes = [axes]
    
    for i, (_, combo) in enumerate(combinations.iterrows()):
        material, condition = combo['Material'], combo['Condition']
        combo_data = df[(df['Material'] == material) & (df['Condition'] == condition)]
        
        ax = axes[i]
        
        # Plot replicates
        replicates = sorted(combo_data['Replicate'].unique())
        for j, rep in enumerate(replicates):
            rep_data = combo_data[combo_data['Replicate'] == rep]
            ax.plot(rep_data["datetime"], rep_data["Value"], 
                   label=f"Rep {rep}", color=colors[j % len(colors)], 
                   linewidth=1, alpha=0.7)
        
        # Plot air temperature
        ax.plot(combo_data["datetime"], combo_data["air_temp_C"], 
               label="Air Temp", color="black", linewidth=1, linestyle='--')
        
        # Formatting
        treatment_title = title_mapping.get((material, condition), f"{material} - {condition}")
        ax.set_title(treatment_title, fontsize=23, fontweight='bold')
        ax.set_ylabel("Temperature (°C)", fontsize=21, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=17, ncol=5)
        ax.tick_params(axis='both', labelsize=15)
        
        # Format x-axis
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
        ax.xaxis.set_major_locator(mdates.MonthLocator())
    
    plt.xlabel("Date", fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.show()

def filter_data_by_date(input_file: str, output_file: str, start_date: str = "2024-10-25") -> str:
    """Filter data from start date and save to new file."""
    print(f"Filtering data from {start_date}...")
    
    df = pd.read_csv(input_file)
    
    print(f"Original: {len(df):,} rows")
    print(f"Date range: {df['timestamp'].min()} to {df['timestamp'].max()}")
    
    # Filter data
    df["datetime"] = pd.to_datetime(df["timestamp"])
    filtered_df = df[df["datetime"] >= start_date].drop('datetime', axis=1)
    
    print(f"Filtered: {len(filtered_df):,} rows")
    print(f"Date range: {filtered_df['timestamp'].min()} to {filtered_df['timestamp'].max()}")
    
    # Save
    filtered_df.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}")
    
    return output_file

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def preview_temperature_data(input_file: str = "merged_data_2025.csv", 
                           start_date: str = "2024-10-25",
                           create_monthly: bool = True, 
                           create_summary: bool = True) -> None:
    """Complete temperature data preview workflow."""
    print("="*60)
    print("TEMPERATURE DATA PREVIEW")
    print("="*60)
    
    # Filter data if needed
    filtered_file = input_file.replace('.csv', '_filtered.csv')
    
    if not os.path.exists(filtered_file):
        filter_data_by_date(input_file, filtered_file, start_date)
    
    # Load and validate
    df = load_and_validate_data(filtered_file, start_date)
    
    # Create plots
    if create_monthly:
        create_monthly_plots(df)
    
    if create_summary:
        create_summary_plot(df)
    
    print("\nPreview complete!")

# ============================================================================
# USAGE
# ============================================================================

if __name__ == "__main__":
    preview_temperature_data(
        input_file="merged_data_2025.csv",
        start_date="2024-10-25",
        create_monthly=True,
        create_summary=True
    )
