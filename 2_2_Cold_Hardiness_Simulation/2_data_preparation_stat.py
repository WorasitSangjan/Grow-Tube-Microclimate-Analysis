# Author: Worasit Sangjan
# Date: July 7, 2025

# ============================================================================
# COLD HARDINESS DATA PREPARATION
# Prepares datasets for by_variety and by_treatment analyses
# ============================================================================

import pandas as pd
import numpy as np
from scipy import stats
import os

# ============================================================================
# CONFIGURATION
# ============================================================================

PHASE_BOUNDARIES = {
    '2023-2024': {
        'fall': (250, 330),
        'winter': [(330, 365), (1, 35)],
        'spring': (35, 90)
    },
    '2024-2025': {
        'fall': (250, 330), 
        'winter': [(330, 365), (1, 50)],
        'spring': (50, 90)
    }
}

VARIETIES = ['Cabernet Sauvignon', 'Chardonnay', 'Concord', 'Mourvedre']
TREATMENTS = ['Paper_low', 'Paper_raise', 'Plastic_low', 'Plastic_raise', 'Uncover_uncover']
SEASONS = ['2023-2024', '2024-2025']

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def calculate_slope(x, y, min_points=10):
    """Calculate slope using linear regression."""
    if len(x) < min_points or len(y) < min_points:
        return None, None, None
    
    valid_indices = ~(np.isnan(x) | np.isnan(y))
    x_clean = np.array(x)[valid_indices]
    y_clean = np.array(y)[valid_indices]
    
    if len(x_clean) < min_points:
        return None, None, None
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_clean, y_clean)
    return slope, r_value**2, len(x_clean)

def get_phase_data(df, phase, season):
    """Extract data for specific phase."""
    phase_boundaries = PHASE_BOUNDARIES[season][phase]
    
    if phase == 'winter':
        mask = ((df['jday'] >= phase_boundaries[0][0]) & (df['jday'] <= phase_boundaries[0][1])) | \
               ((df['jday'] >= phase_boundaries[1][0]) & (df['jday'] <= phase_boundaries[1][1]))
    else:
        mask = (df['jday'] >= phase_boundaries[0]) & (df['jday'] <= phase_boundaries[1])
    
    return df[mask]

def load_variety_data(variety, season):
    """Load data for specific variety and season."""
    season_folder = f"{season.replace('-', '_')}_coldhardiness"
    filename = f'{season_folder}/combined_detailed_results_{variety}.csv'
    
    try:
        return pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Warning: {filename} not found")
        return None

# ============================================================================
# DATASET CREATION
# ============================================================================

def create_dynamic_slopes_database():
    """Create slope database for fall and spring phases."""
    slope_data = []
    
    for season in SEASONS:
        for variety in VARIETIES:
            df = load_variety_data(variety, season)
            if df is None:
                continue
            
            for treatment in TREATMENTS:
                treatment_data = df[df['Treatment'] == treatment]
                if treatment_data.empty:
                    continue
                
                for replicate in treatment_data['Replicate'].unique():
                    rep_data = treatment_data[treatment_data['Replicate'] == replicate]
                    if rep_data.empty:
                        continue
                    
                    for phase in ['fall', 'spring']:
                        phase_data = get_phase_data(rep_data, phase, season)
                        
                        if not phase_data.empty:
                            slope, r2, n_points = calculate_slope(
                                phase_data['jday'].values, 
                                phase_data['Predicted_Hc'].values
                            )
                            
                            slope_data.append({
                                'Season': season,
                                'Variety': variety,
                                'Treatment': treatment,
                                'Replicate': replicate,
                                'Phase': phase,
                                'Slope_C_per_day': slope,
                                'R_squared': r2,
                                'N_points': n_points
                            })
    
    slopes_df = pd.DataFrame(slope_data)
    slopes_df.to_csv('dynamic_slopes_database.csv', index=False)
    return slopes_df

def create_winter_means_database():
    """Create winter hardiness database."""
    winter_data = []
    
    for season in SEASONS:
        for variety in VARIETIES:
            df = load_variety_data(variety, season)
            if df is None:
                continue
            
            for treatment in TREATMENTS:
                treatment_data = df[df['Treatment'] == treatment]
                if treatment_data.empty:
                    continue
                
                for replicate in treatment_data['Replicate'].unique():
                    rep_data = treatment_data[treatment_data['Replicate'] == replicate]
                    if rep_data.empty:
                        continue
                    
                    winter_phase_data = get_phase_data(rep_data, 'winter', season)
                    
                    if not winter_phase_data.empty and len(winter_phase_data) >= 10:
                        winter_data.append({
                            'Season': season,
                            'Variety': variety,
                            'Treatment': treatment,
                            'Replicate': replicate,
                            'Phase': 'winter',
                            'Mean_Hardiness': winter_phase_data['Predicted_Hc'].mean(),
                            'Std_Hardiness': winter_phase_data['Predicted_Hc'].std(),
                            'Min_Hardiness': winter_phase_data['Predicted_Hc'].min(),
                            'Max_Hardiness': winter_phase_data['Predicted_Hc'].max(),
                            'N_points': len(winter_phase_data)
                        })
    
    winter_df = pd.DataFrame(winter_data)
    winter_df.to_csv('winter_means_database.csv', index=False)
    return winter_df

def create_analysis_datasets(dynamic_df, winter_df):
    """Create analysis-specific datasets."""
    # By variety datasets
    dynamic_by_variety = dynamic_df.copy()
    dynamic_by_variety['Analysis_Type'] = 'by_variety'
    dynamic_by_variety.to_csv('dynamic_slopes_by_variety.csv', index=False)
    
    winter_by_variety = winter_df.copy()  
    winter_by_variety['Analysis_Type'] = 'by_variety'
    winter_by_variety.to_csv('winter_means_by_variety.csv', index=False)
    
    # By treatment datasets
    dynamic_by_treatment = dynamic_df.copy()
    dynamic_by_treatment['Analysis_Type'] = 'by_treatment'
    dynamic_by_treatment.to_csv('dynamic_slopes_by_treatment.csv', index=False)
    
    winter_by_treatment = winter_df.copy()
    winter_by_treatment['Analysis_Type'] = 'by_treatment'  
    winter_by_treatment.to_csv('winter_means_by_treatment.csv', index=False)

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    """Main data preparation workflow."""
    print("Creating dynamic slopes database...")
    dynamic_slopes = create_dynamic_slopes_database()
    print(f"Dynamic slopes: {len(dynamic_slopes)} observations")
    
    print("Creating winter means database...")
    winter_means = create_winter_means_database()
    print(f"Winter means: {len(winter_means)} observations")
    
    print("Creating analysis-specific datasets...")
    create_analysis_datasets(dynamic_slopes, winter_means)
    
    print("Data preparation complete!")
    return dynamic_slopes, winter_means

if __name__ == "__main__":
    main()