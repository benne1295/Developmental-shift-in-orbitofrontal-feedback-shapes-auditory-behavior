"""
Simple script to plot bouton density from Excel file
Filters axons with density = 0 or < 0.05 per 10 microns, then plots
"""

import os
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.formula.api import mixedlm

# ============================================================================
# CONFIGURATION - MODIFY THESE PARAMETERS
# ============================================================================

# File paths
ADOLESCENT_RESULTS = '/Users/bennepraegel/Downloads/adolescent_filtered/consolidated_axon_peaks_interactive.xlsx'
ADULT_RESULTS = '/Users/bennepraegel/Downloads/adult_filtered/consolidated_axon_peaks_interactive.xlsx'
OUTPUT_DIR = '/Users/bennepraegel/Downloads/analysis_results'

# Density calculation: per how many microns to calculate density (e.g., 10 = per 10 microns)
DENSITY_PER_MICRONS = 5.0

# Filter threshold: minimum density required to include axon (e.g., 0.05 = 0.05 boutons per DENSITY_PER_MICRONS)
MIN_DENSITY_THRESHOLD = 0.2

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================

def calculate_density(df, microns_per_density):
    """
    Calculate bouton density per specified microns.
    Uses DENSITY_PER_MICRONS parameter from configuration.e
    
    Args:
        df: DataFrame with peak count and axon length data
        microns_per_density: Number of microns to calculate density per (from DENSITY_PER_MICRONS config)
    """
    peak_count_col = None
    for col in ['PeakCount', 'SpikeCount_Total']:
        if col in df.columns:
            peak_count_col = col
            break
    
    # Column name for density - uses DENSITY_PER_MICRONS value
    density_col_name = f'Density_per_{int(microns_per_density)}um'
    
    # Check if exact column exists
    if density_col_name in df.columns:
        df[density_col_name] = df[density_col_name]
    elif f'BoutonDensity_per_{int(microns_per_density)}um' in df.columns:
        df[density_col_name] = df[f'BoutonDensity_per_{int(microns_per_density)}um']
    elif f'SpikeDensity_per_{int(microns_per_density)}um' in df.columns:
        df[density_col_name] = df[f'SpikeDensity_per_{int(microns_per_density)}um']
    elif 'SpikeDensity_per_um' in df.columns:
        # CALCULATION: Uses DENSITY_PER_MICRONS here
        df[density_col_name] = df['SpikeDensity_per_um'] * microns_per_density
    elif 'BoutonDensity_per_um' in df.columns:
        # CALCULATION: Uses DENSITY_PER_MICRONS here
        df[density_col_name] = df['BoutonDensity_per_um'] * microns_per_density
    elif peak_count_col and 'AxonLength_um' in df.columns:
        # CALCULATION: Uses DENSITY_PER_MICRONS here
        # Formula: (peak_count / length) * DENSITY_PER_MICRONS
        df[density_col_name] = (df[peak_count_col] / df['AxonLength_um']) * microns_per_density
        print(f"    Calculated density using: ({peak_count_col} / AxonLength_um) * {microns_per_density}")
    else:
        raise ValueError(f"Could not find appropriate columns to calculate density per {microns_per_density} microns")
    
    return df, density_col_name


def extract_mouse_id_from_brain_region(brain_region):
    """
    Extract mouse number from BrainRegion string.
    Format: "1_sl1" -> mouse "1", "2_sl2" -> mouse "2"
    Handles: "1_sl1", "2_sl1", "10_sl2", etc.
    
    Args:
        brain_region: BrainRegion string (e.g., "1_sl1", "2_sl2", "10_sl3")
    
    Returns:
        Mouse ID string (e.g., "1", "2", "10") or None if extraction fails
    """
    if pd.isna(brain_region) or brain_region == '':
        return None
    
    brain_region_str = str(brain_region).strip()
    
    # Extract mouse number before first underscore (if underscore exists)
    # e.g., "2_sl1" -> "2", "10_sl2" -> "10"
    if '_' in brain_region_str:
        mouse_id = brain_region_str.split('_')[0]
        # Verify it's a number (can be multi-digit like "10")
        if mouse_id.isdigit():
            return mouse_id
        return None
    else:
        # If no underscore, try to extract first number from string
        match = re.search(r'^\d+', brain_region_str)
        if match:
            return match.group(0)
        return None


def load_and_filter(microns_per_density, min_threshold):
    """
    Load data, calculate density, filter, and save filtered Excel.
    
    Args:
        microns_per_density: Number of microns to calculate density per
        min_threshold: Minimum density threshold for filtering
    """
    print("="*70)
    print("LOADING DATA")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Density calculation: per {microns_per_density} microns")
    print(f"  Minimum density threshold: {min_threshold}")
    
    # Load
    print("\nLoading data...")
    df_adolescent = pd.read_excel(ADOLESCENT_RESULTS)
    df_adult = pd.read_excel(ADULT_RESULTS)
    print(f"  Adolescent: {len(df_adolescent)} axons")
    print(f"  Adult: {len(df_adult)} axons")
    
    # Calculate density - USES DENSITY_PER_MICRONS parameter
    print(f"\nCalculating density per {microns_per_density} microns (using DENSITY_PER_MICRONS = {microns_per_density})...")
    df_adolescent, density_col = calculate_density(df_adolescent, microns_per_density)
    df_adult, _ = calculate_density(df_adult, microns_per_density)
    print(f"    Created density column: {density_col}")
    
    # Filter: remove density = 0 or < threshold - USES MIN_DENSITY_THRESHOLD parameter
    print(f"\n{'='*70}")
    print(f"FILTERING: Removing axons with density = 0 or < {min_threshold} (using MIN_DENSITY_THRESHOLD = {min_threshold})")
    print(f"{'='*70}")
    
    before_adolescent = len(df_adolescent)
    # FILTER CONDITION: Uses MIN_DENSITY_THRESHOLD here
    df_adolescent = df_adolescent[
        (df_adolescent[density_col] > 0) & 
        (df_adolescent[density_col] >= min_threshold) &  # MIN_DENSITY_THRESHOLD applied here
        (df_adolescent[density_col].notna())
    ].copy()
    removed_adolescent = before_adolescent - len(df_adolescent)
    
    before_adult = len(df_adult)
    # FILTER CONDITION: Uses MIN_DENSITY_THRESHOLD here
    df_adult = df_adult[
        (df_adult[density_col] > 0) & 
        (df_adult[density_col] >= min_threshold) &  # MIN_DENSITY_THRESHOLD applied here
        (df_adult[density_col].notna())
    ].copy()
    removed_adult = before_adult - len(df_adult)
    
    print(f"  Adolescent: Removed {removed_adolescent} axons ({len(df_adolescent)} remaining)")
    print(f"  Adult: Removed {removed_adult} axons ({len(df_adult)} remaining)")
    print(f"  Total removed: {removed_adolescent + removed_adult} axons")
    print(f"  Total remaining: {len(df_adolescent) + len(df_adult)} axons")
    if len(df_adolescent) > 0 or len(df_adult) > 0:
        df_combined_check = pd.concat([df_adolescent, df_adult], ignore_index=True)
        print(f"  Minimum density in filtered data: {df_combined_check[density_col].min():.4f}")
        print(f"  Maximum density in filtered data: {df_combined_check[density_col].max():.4f}")
    
    # Add age group
    df_adolescent['AgeGroup'] = 'Adolescent'
    df_adult['AgeGroup'] = 'Adult'
    
    # Extract MouseID from BrainRegion - ALWAYS re-extract to ensure correctness
    print(f"\n{'='*70}")
    print("EXTRACTING MOUSE ID FROM BRAINREGION (Column B)")
    print(f"{'='*70}")
    print("Note: Extracting mouse number (e.g., '1', '2') from BrainRegion (e.g., '2_sl1' -> '2')")
    
    # Always extract MouseID from BrainRegion (even if column exists, update it)
    if 'BrainRegion' in df_adolescent.columns:
        print(f"\nAdolescent data:")
        print(f"  Total rows: {len(df_adolescent)}")
        print(f"  Unique BrainRegion values: {df_adolescent['BrainRegion'].nunique()}")
        
        df_adolescent['MouseID'] = df_adolescent['BrainRegion'].apply(extract_mouse_id_from_brain_region)
        
        # Show all unique BrainRegion -> MouseID mappings
        unique_mappings = df_adolescent[['BrainRegion', 'MouseID']].drop_duplicates().sort_values('BrainRegion')
        print(f"\n  All BrainRegion -> MouseID mappings:")
        for _, row in unique_mappings.iterrows():
            br = row['BrainRegion']
            mouse_id = row['MouseID']
            count = len(df_adolescent[df_adolescent['BrainRegion'] == br])
            print(f"    '{br}' -> MouseID '{mouse_id}' ({count} axons)")
        
        # Check for any failed extractions
        failed = df_adolescent[df_adolescent['MouseID'].isna()]
        if len(failed) > 0:
            print(f"\n  WARNING: {len(failed)} rows with failed MouseID extraction:")
            for br in failed['BrainRegion'].unique():
                print(f"    '{br}' -> Could not extract MouseID")
    
    if 'BrainRegion' in df_adult.columns:
        print(f"\nAdult data:")
        print(f"  Total rows: {len(df_adult)}")
        print(f"  Unique BrainRegion values: {df_adult['BrainRegion'].nunique()}")
        
        df_adult['MouseID'] = df_adult['BrainRegion'].apply(extract_mouse_id_from_brain_region)
        
        # Show all unique BrainRegion -> MouseID mappings
        unique_mappings = df_adult[['BrainRegion', 'MouseID']].drop_duplicates().sort_values('BrainRegion')
        print(f"\n  All BrainRegion -> MouseID mappings:")
        for _, row in unique_mappings.iterrows():
            br = row['BrainRegion']
            mouse_id = row['MouseID']
            count = len(df_adult[df_adult['BrainRegion'] == br])
            print(f"    '{br}' -> MouseID '{mouse_id}' ({count} axons)")
        
        # Check for any failed extractions
        failed = df_adult[df_adult['MouseID'].isna()]
        if len(failed) > 0:
            print(f"\n  WARNING: {len(failed)} rows with failed MouseID extraction:")
            for br in failed['BrainRegion'].unique():
                print(f"    '{br}' -> Could not extract MouseID")
    
    print(f"\n{'='*70}")
    
    # Create unique mouse identifier combining AgeGroup and MouseID
    # This ensures mice from different age groups are treated as different groups
    df_adolescent['MouseID_unique'] = 'Adolescent_' + df_adolescent['MouseID'].astype(str)
    df_adult['MouseID_unique'] = 'Adult_' + df_adult['MouseID'].astype(str)
    
    # Combine
    df_combined = pd.concat([df_adolescent, df_adult], ignore_index=True)
    
    # Show mouse statistics
    if 'MouseID' in df_combined.columns:
        mice_per_group = df_combined.groupby('AgeGroup')['MouseID'].nunique()
        print(f"\nMouse statistics:")
        print(f"  Adolescent: {mice_per_group.get('Adolescent', 0)} unique mice")
        print(f"  Adult: {mice_per_group.get('Adult', 0)} unique mice")
        axons_per_mouse = df_combined.groupby(['AgeGroup', 'MouseID']).size()
        if 'Adolescent' in axons_per_mouse.index.get_level_values(0):
            avg_adolescent = axons_per_mouse.xs('Adolescent', level='AgeGroup').mean()
            print(f"  Average axons per mouse - Adolescent: {avg_adolescent:.1f}")
        if 'Adult' in axons_per_mouse.index.get_level_values(0):
            avg_adult = axons_per_mouse.xs('Adult', level='AgeGroup').mean()
            print(f"  Average axons per mouse - Adult: {avg_adult:.1f}")
    
    # Save filtered Excel files
    print(f"\n{'='*70}")
    print("SAVING FILTERED EXCEL FILES")
    print(f"{'='*70}")
    df_adolescent.to_excel(os.path.join(OUTPUT_DIR, 'adolescent_filtered_results.xlsx'), index=False)
    df_adult.to_excel(os.path.join(OUTPUT_DIR, 'adult_filtered_results.xlsx'), index=False)
    df_combined.to_excel(os.path.join(OUTPUT_DIR, 'combined_filtered_results.xlsx'), index=False)
    print(f"  Saved filtered Excel files to {OUTPUT_DIR}")
    
    return df_combined, density_col


def fit_mixed_effects_model(df_subset, density_col):
    """
    Fit Linear Mixed Effects Model with SlideID as random effect.
    
    Model Structure:
    - Response variable: Density (boutons per X microns)
    - Fixed effect: AgeGroup (Adolescent vs Adult)
    - Random effect: MouseID (random intercept per mouse)
    - This accounts for: Multiple axons per mouse, variable sample sizes per mouse
    
    Formula: Density ~ AgeGroup + (1|MouseID)
    - AgeGroup is a fixed factor (comparison between Adolescent and Adult)
    - (1|MouseID) means random intercept varies by mouse
    - Each mouse gets its own baseline/intercept, but the AgeGroup effect is fixed
    
    Args:
        df_subset: DataFrame with density data
        density_col: Name of the density column to analyze
    
    Returns:
        p_value: p-value for AgeGroup effect
        model_result: Fitted model result object
        model_info: Dictionary with model details
    """
    # Check if MouseID or MouseID_unique column exists
    mouse_id_col = 'MouseID_unique' if 'MouseID_unique' in df_subset.columns else 'MouseID'
    
    if mouse_id_col not in df_subset.columns:
        print(f"    Warning: MouseID column not found, falling back to Mann-Whitney U test")
        adolescent = df_subset[df_subset['AgeGroup'] == 'Adolescent'][density_col].dropna()
        adult = df_subset[df_subset['AgeGroup'] == 'Adult'][density_col].dropna()
        if len(adolescent) > 0 and len(adult) > 0:
            _, p_val = stats.mannwhitneyu(adolescent.values, adult.values, alternative='two-sided')
            print(f"    Using Mann-Whitney U test, p-value: {p_val:.6f}")
            return p_val, None, None
        else:
            return np.nan, None, None
    
    print(f"    Using {mouse_id_col} for random effects")
    
    # Prepare data for model
    df_model = df_subset[[density_col, 'AgeGroup', mouse_id_col]].copy()
    df_model = df_model.dropna(subset=[density_col, 'AgeGroup', mouse_id_col])
    
    if len(df_model) == 0:
        return np.nan, None, None
    
    # Check if we have multiple mice per group
    mice_per_group = df_model.groupby('AgeGroup')[mouse_id_col].nunique()
    print(f"    Mice per group: {mice_per_group.to_dict()}")
    if mice_per_group.min() < 2:
        print(f"    Warning: Not enough mice per group for LMM (min={mice_per_group.min()}), falling back to Mann-Whitney U test")
        adolescent = df_subset[df_subset['AgeGroup'] == 'Adolescent'][density_col].dropna()
        adult = df_subset[df_subset['AgeGroup'] == 'Adult'][density_col].dropna()
        if len(adolescent) > 0 and len(adult) > 0:
            _, p_val = stats.mannwhitneyu(adolescent.values, adult.values, alternative='two-sided')
            print(f"    Using Mann-Whitney U test, p-value: {p_val:.6f}")
            return p_val, None, None
        else:
            return np.nan, None, None
    
    try:
        # Create formula: density ~ AgeGroup with random intercept per MouseID
        # Fit mixed effects model with random intercept per mouse
        # This accounts for multiple axons per mouse and variable sample sizes
        formula = f'{density_col} ~ C(AgeGroup)'
        print(f"    Fitting LMM: {formula} with random intercept per {mouse_id_col}")
        print(f"    Model structure: Density ~ AgeGroup + (1|{mouse_id_col})")
        print(f"    Data: {len(df_model)} observations, {df_model[mouse_id_col].nunique()} unique mice")
        
        model = mixedlm(formula, df_model, groups=df_model[mouse_id_col])
        result = model.fit(reml=True)
        
        # Print model summary for debugging
        print(f"    LMM fitted successfully")
        print(f"    Model coefficients: {result.params}")
        print(f"    P-values: {result.pvalues}")
        
        # Create model info dictionary
        model_info = {
            'formula': formula,
            'random_effect': mouse_id_col,
            'n_observations': len(df_model),
            'n_mice': df_model[mouse_id_col].nunique(),  # Number of unique mice
            'mice_per_group': mice_per_group.to_dict(),  # Mice per age group
            'coefficients': result.params.to_dict(),
            'pvalues': result.pvalues.to_dict(),
            'aic': result.aic,
            'bic': result.bic,
            'llf': result.llf,
            'scale': result.scale,
            'fe_params': result.fe_params.to_dict(),
            're_params': result.cov_re.to_dict() if hasattr(result, 'cov_re') else {},
            'model_summary': str(result.summary())
        }
        
        # Extract p-value for AgeGroup effect
        # The p-value is for the coefficient of AgeGroup (comparing Adult to Adolescent)
        if 'C(AgeGroup)[T.Adult]' in result.pvalues.index:
            p_val = result.pvalues['C(AgeGroup)[T.Adult]']
            print(f"    Using p-value from C(AgeGroup)[T.Adult]: {p_val:.6f}")
        elif len(result.pvalues) > 1:
            # Usually the second coefficient is AgeGroup (first is Intercept)
            p_val = result.pvalues.iloc[1]
            print(f"    Using p-value from coefficient index 1: {p_val:.6f}")
        else:
            p_val = np.nan
            print(f"    Warning: Could not extract p-value from model")
        
        return p_val, result, model_info
    
    except Exception as e:
        print(f"    Warning: LMM failed ({str(e)}), falling back to Mann-Whitney U test")
        adolescent = df_subset[df_subset['AgeGroup'] == 'Adolescent'][density_col].dropna()
        adult = df_subset[df_subset['AgeGroup'] == 'Adult'][density_col].dropna()
        if len(adolescent) > 0 and len(adult) > 0:
            _, p_val = stats.mannwhitneyu(adolescent.values, adult.values, alternative='two-sided')
            return p_val, None
        else:
            return np.nan, None


def plot_single_density(df_subset, density_col, microns_per_density, title, filename_base):
    """
    Helper function to plot a single density comparison.
    
    Args:
        df_subset: DataFrame subset to plot
        density_col: Name of the density column to plot
        microns_per_density: Number of microns for density calculation (for labels)
        title: Title for the plot
        filename_base: Base filename for saving (without extension)
    
    Returns:
        model_info: Dictionary with LME model information (or None if not available)
    """
    if len(df_subset) == 0:
        print(f"  ⚠ No data for {title} - skipping plot")
        return None
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Violin plot
    sns.violinplot(data=df_subset, x='AgeGroup', y=density_col,
                   ax=ax, palette=['lightgreen', 'darkgreen'])
    
    # Scatter overlay
    sns.stripplot(data=df_subset, x='AgeGroup', y=density_col,
                  ax=ax, color='black', alpha=0.5, size=4)
    
    # Statistics using Linear Mixed Effects Model
    adolescent = df_subset[df_subset['AgeGroup'] == 'Adolescent'][density_col].dropna()
    adult = df_subset[df_subset['AgeGroup'] == 'Adult'][density_col].dropna()
    
    model_info = None
    if len(adolescent) > 0 and len(adult) > 0:
        # Fit Linear Mixed Effects Model
        p_val, model_result, model_info = fit_mixed_effects_model(df_subset, density_col)
        
        # Calculate descriptive statistics for display
        y_max = df_subset[density_col].max()
        y_range = y_max - df_subset[density_col].min()
        
        if not np.isnan(p_val):
            if p_val < 0.001:
                sig_text, bg_color = '***', 'red'
            elif p_val < 0.01:
                sig_text, bg_color = '**', 'orange'
            elif p_val < 0.05:
                sig_text, bg_color = '*', 'yellow'
            else:
                sig_text, bg_color = 'ns', 'lightgray'
            
            # Add p-value annotation
            ax.text(0.5, y_max + y_range * 0.15, 
                   f"{sig_text}\np={p_val:.4f} (LMM)", 
                   ha='center', fontsize=12, weight='bold',
                   bbox=dict(boxstyle='round', facecolor=bg_color, alpha=0.9, 
                            edgecolor='black', linewidth=3))
        else:
            # If p-value couldn't be calculated
            ax.text(0.5, y_max + y_range * 0.15, 
                   "p=NA (LMM failed)", 
                   ha='center', fontsize=12, weight='bold',
                   bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.9, 
                            edgecolor='black', linewidth=3))
        
        # Add stats (showing number of axons and mice per group)
        for j, age in enumerate(['Adolescent', 'Adult']):
            age_data = df_subset[df_subset['AgeGroup'] == age][density_col].dropna()
            if len(age_data) > 0:
                mean_val = age_data.mean()
                sem_val = age_data.sem()
                n_axons = len(age_data)
                
                # Count number of mice in this age group
                mouse_id_col = 'MouseID_unique' if 'MouseID_unique' in df_subset.columns else 'MouseID'
                if mouse_id_col in df_subset.columns:
                    n_mice = df_subset[df_subset['AgeGroup'] == age][mouse_id_col].nunique()
                    stats_text = f"n={n_axons} axons\nm={n_mice} mice\nμ={mean_val:.2f}\n±{sem_val:.2f}"
                else:
                    stats_text = f"n={n_axons}\nμ={mean_val:.2f}\n±{sem_val:.2f}"
                
                ax.text(j, y_max + y_range * 0.05, 
                       stats_text,
                       ha='center', fontsize=9,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='gray'))
        
        ax.set_ylim(0, y_max + y_range * 0.30)
    
    # Format microns for display (show as integer if whole number, otherwise show decimal)
    if microns_per_density == int(microns_per_density):
        microns_display = int(microns_per_density)
    else:
        microns_display = microns_per_density
    ax.set_ylabel(f'Bouton Density (boutons per {microns_display} µm)', fontsize=14, weight='bold')
    ax.set_xlabel('Age Group', fontsize=13)
    ax.set_title(title, fontsize=15, weight='bold', pad=20)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    
    # Save
    png_path = os.path.join(OUTPUT_DIR, f'{filename_base}.png')
    svg_path = os.path.join(OUTPUT_DIR, f'{filename_base}.svg')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(svg_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {png_path}")
    
    # Return model_info if available
    return model_info


def save_lme_results(all_model_results, output_dir):
    """
    Save LME model results to Excel file.
    
    Args:
        all_model_results: Dictionary with keys as analysis names and values as model_info dictionaries
        output_dir: Directory to save the Excel file
    """
    if not all_model_results:
        print("\nNo LME model results to save.")
        return
    
    print(f"\n{'='*70}")
    print("SAVING LME MODEL RESULTS")
    print(f"{'='*70}")
    
    # Create a list to store flattened results
    results_list = []
    
    for analysis_name, model_info in all_model_results.items():
        if model_info is None:
            continue
        
        # Extract key information for the table
        result_row = {
            'Analysis': analysis_name,
            'Formula': model_info.get('formula', 'N/A'),
            'Random_Effect': model_info.get('random_effect', 'N/A'),
            'N_Observations': model_info.get('n_observations', np.nan),
            'N_Mice': model_info.get('n_mice', np.nan),
            'AIC': model_info.get('aic', np.nan),
            'BIC': model_info.get('bic', np.nan),
            'Log_Likelihood': model_info.get('llf', np.nan),
            'Scale': model_info.get('scale', np.nan),
        }
        
        # Add mice per group
        mice_per_group = model_info.get('mice_per_group', {})
        result_row['N_Mice_Adolescent'] = mice_per_group.get('Adolescent', np.nan)
        result_row['N_Mice_Adult'] = mice_per_group.get('Adult', np.nan)
        
        # Add coefficients
        coefficients = model_info.get('coefficients', {})
        for coef_name, coef_value in coefficients.items():
            result_row[f'Coefficient_{coef_name}'] = coef_value
        
        # Add p-values
        pvalues = model_info.get('pvalues', {})
        for pval_name, pval_value in pvalues.items():
            result_row[f'PValue_{pval_name}'] = pval_value
        
        # Add standard errors and confidence intervals if available from fe_params
        fe_params = model_info.get('fe_params', {})
        if isinstance(fe_params, dict):
            for param_name, param_value in fe_params.items():
                result_row[f'FixedEffect_{param_name}'] = param_value
        
        # Add random effect variance if available
        re_params = model_info.get('re_params', {})
        if isinstance(re_params, dict) and re_params:
            for re_name, re_value in re_params.items():
                result_row[f'RandomEffect_Variance_{re_name}'] = re_value
        
        results_list.append(result_row)
    
    if not results_list:
        print("No valid LME model results to save.")
        return
    
    # Create DataFrame
    df_results = pd.DataFrame(results_list)
    
    # Save to Excel
    excel_path = os.path.join(output_dir, 'LME_model_results.xlsx')
    df_results.to_excel(excel_path, index=False)
    print(f"  Saved LME results table to: {excel_path}")
    
    # Also save detailed model summaries as text files
    for analysis_name, model_info in all_model_results.items():
        if model_info is None:
            continue
        
        model_summary = model_info.get('model_summary', '')
        if model_summary:
            # Clean filename
            safe_name = analysis_name.lower().replace(' ', '_').replace('/', '_')
            txt_path = os.path.join(output_dir, f'LME_model_summary_{safe_name}.txt')
            with open(txt_path, 'w') as f:
                f.write(f"Linear Mixed Effects Model Results: {analysis_name}\n")
                f.write("=" * 80 + "\n\n")
                f.write("MODEL STRUCTURE:\n")
                f.write("-" * 80 + "\n")
                f.write(f"Formula: {model_info.get('formula', 'N/A')}\n")
                f.write(f"Random Effect: {model_info.get('random_effect', 'N/A')}\n")
                f.write(f"Model: Density ~ AgeGroup + (1|{model_info.get('random_effect', 'MouseID')})\n\n")
                f.write("EXPLANATION:\n")
                f.write("-" * 80 + "\n")
                f.write("Response variable: Bouton Density (boutons per X microns)\n")
                f.write("Fixed effect: AgeGroup (Adolescent vs Adult)\n")
                f.write("Random effect: Random intercept per mouse (accounts for multiple axons per mouse)\n")
                f.write("This model accounts for:\n")
                f.write("  - Multiple axons per mouse (random intercept per mouse)\n")
                f.write("  - Variable sample sizes per mouse\n")
                f.write("  - Non-independence of axons from the same mouse\n\n")
                f.write("MODEL RESULTS:\n")
                f.write("-" * 80 + "\n")
                f.write(model_summary)
            print(f"  Saved detailed summary to: {txt_path}")


def plot_density(df, density_col, microns_per_density):
    """
    Plot density column as violin plot with scatter for all axons and separately for each layer.
    
    Args:
        df: DataFrame with density data
        density_col: Name of the density column to plot
        microns_per_density: Number of microns for density calculation (for labels)
    """
    print(f"\n{'='*70}")
    print("CREATING PLOTS")
    print(f"{'='*70}")
    
    # Format microns for display (show as integer if whole number, otherwise show decimal)
    if microns_per_density == int(microns_per_density):
        microns_display = int(microns_per_density)
    else:
        microns_display = microns_per_density
    
    # Normalize layer names
    df_plot = df.copy()
    if 'Layer' in df_plot.columns:
        df_plot['Layer'] = df_plot['Layer'].replace({
            'L1': 'Supragranular',
            'L5': 'Infragranular'
        })
    
    # Filter for layers that exist
    df_layers = df_plot[df_plot['Layer'].isin(['Supragranular', 'Infragranular'])].copy() if 'Layer' in df_plot.columns else pd.DataFrame()
    
    # Collect LME model results from all plots
    all_model_results = {}
    
    # Plot 1: All axons
    print("\nPlot 1: All axons")
    model_info = plot_single_density(
        df_plot,
        density_col,
        microns_per_density,
        f'Bouton Density per {microns_display} µm - All Axons\nAdolescent vs Adult',
        'bouton_density_all'
    )
    all_model_results['All Axons'] = model_info
    
    # Plot 2: Infragranular + Supragranular combined
    if len(df_layers) > 0:
        print("\nPlot 2: Infragranular + Supragranular layers combined")
        model_info = plot_single_density(
            df_layers,
            density_col,
            microns_per_density,
            f'Bouton Density per {microns_display} µm - Infragranular + Supragranular Combined\nAdolescent vs Adult',
            'bouton_density_infra_supra_combined'
        )
        all_model_results['Infragranular + Supragranular Combined'] = model_info
    
    # Plot 3: Infragranular (L5)
    if len(df_layers) > 0:
        df_l5 = df_layers[df_layers['Layer'] == 'Infragranular'].copy()
        print("\nPlot 3: Infragranular layers (L5)")
        model_info = plot_single_density(
            df_l5,
            density_col,
            microns_per_density,
            f'Bouton Density per {microns_display} µm - Infragranular Layers (L5)\nAdolescent vs Adult',
            'bouton_density_infragranular_L5'
        )
        all_model_results['Infragranular (L5)'] = model_info
    
    # Plot 4: Supragranular (L1)
    if len(df_layers) > 0:
        df_l1 = df_layers[df_layers['Layer'] == 'Supragranular'].copy()
        print("\nPlot 4: Supragranular layers (L1)")
        model_info = plot_single_density(
            df_l1,
            density_col,
            microns_per_density,
            f'Bouton Density per {microns_display} µm - Supragranular Layers (L1)\nAdolescent vs Adult',
            'bouton_density_supragranular_L1'
        )
        all_model_results['Supragranular (L1)'] = model_info
    
    # Save all LME model results
    save_lme_results(all_model_results, OUTPUT_DIR)


def print_lme_structure_explanation():
    """
    Print explanation of the Linear Mixed Effects Model structure.
    """
    print("\n" + "="*80)
    print("LINEAR MIXED EFFECTS MODEL (LME) STRUCTURE EXPLANATION")
    print("="*80)
    print("\nMODEL FORMULA:")
    print("  Density ~ AgeGroup + (1|MouseID)")
    print("\nCOMPONENTS:")
    print("  1. Response Variable (Density):")
    print("     - Bouton density (boutons per X microns)")
    print("     - Continuous numerical value")
    print("\n  2. Fixed Effect (AgeGroup):")
    print("     - Categorical predictor: Adolescent vs Adult")
    print("     - Tests whether there's a significant difference between age groups")
    print("     - The coefficient represents the difference in density between Adult and Adolescent")
    print("     - The p-value tests if this difference is statistically significant")
    print("\n  3. Random Effect (1|MouseID):")
    print("     - Random intercept per mouse (MouseID)")
    print("     - Accounts for the fact that multiple axons come from the same mouse")
    print("     - Each mouse gets its own baseline/intercept")
    print("     - Allows the model to account for within-mouse correlation")
    print("\nWHY USE LME INSTEAD OF STANDARD TESTS?")
    print("  - Multiple axons per mouse: Axons from the same mouse are not independent")
    print("  - Variable sample sizes: Different mice have different numbers of axons")
    print("  - Proper statistical inference: Accounts for nested structure (axons within mice)")
    print("  - More appropriate: Standard tests assume independence, which is violated here")
    print("\nINTERPRETATION:")
    print("  - The model estimates an overall AgeGroup effect while accounting for")
    print("    mouse-to-mouse variation")
    print("  - The p-value indicates if the AgeGroup difference is significant after")
    print("    controlling for the random mouse effects")
    print("  - Coefficients show the estimated difference between age groups")
    print("\n" + "="*80 + "\n")


def main():
    """Main function."""
    print("="*70)
    print("BOUTON DENSITY ANALYSIS")
    print("="*70)
    print(f"\nUsing configuration parameters:")
    print(f"  DENSITY_PER_MICRONS = {DENSITY_PER_MICRONS}")
    print(f"  MIN_DENSITY_THRESHOLD = {MIN_DENSITY_THRESHOLD}")
    
    # Print LME structure explanation
    print_lme_structure_explanation()
    
    # Load and filter using configuration parameters
    # DENSITY_PER_MICRONS is passed to calculate_density() function
    # MIN_DENSITY_THRESHOLD is used in the filtering step
    df, density_col = load_and_filter(DENSITY_PER_MICRONS, MIN_DENSITY_THRESHOLD)
    
    # Plot - DENSITY_PER_MICRONS is used for axis labels and titles
    plot_density(df, density_col, DENSITY_PER_MICRONS)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)


if __name__ == "__main__":
    main()
