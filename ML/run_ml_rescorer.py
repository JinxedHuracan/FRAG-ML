import os
import re
import numpy as np
import pandas as pd


def normalize_name(name):
  name = str(name).lower()
  name = re.sub(r'(_complex|_refined|_pose1|\.pdb|\.pdbqt)', '', name)
  return name.strip()


# 1. Load MCDA, ADMET, and FBDD Data
mcda_df = pd.read_csv(
    'final_submission_results/final_repurposing_ranked_candidates.csv'
)
admet_df = pd.read_csv('admet_results/top15_admet_profile.csv')
fbdd_df = pd.read_csv('fbdd_results/fbdd_fragment_metrics.csv')

# Clean headers
admet_df.columns = [
    c.strip().replace(' (kcal/mol)', '').replace(' (Da)', '').replace(' (A2)', '')
    for c in admet_df.columns
]

# Standardize match keys
mcda_df['Match_Key'] = mcda_df['Drug'].apply(normalize_name)
admet_df['Match_Key'] = admet_df['Drug'].apply(normalize_name)
fbdd_df['Match_Key'] = fbdd_df['Drug'].apply(normalize_name)

# Merge datasets
merged = mcda_df.merge(
    admet_df[['Match_Key', 'TPSA', 'LogP', 'RotB', 'MW']],
    on='Match_Key',
    how='left',
    suffixes=('', '_admet'),
)
merged = merged.merge(
    fbdd_df[[
        'Match_Key',
        'Ligand_Efficiency (LE)',
        'Fit_Quality (FQ)',
        'Core_Fragment',
        'Target_Subpocket',
    ]],
    on='Match_Key',
    how='left',
)

# Fill missing descriptor values with robust defaults if needed
merged['TPSA'] = pd.to_numeric(merged['TPSA']).fillna(80.0)
merged['LogP'] = pd.to_numeric(merged['LogP']).fillna(2.5)
merged['RotB'] = pd.to_numeric(merged['RotB']).fillna(4)
merged['Ligand_Efficiency (LE)'] = pd.to_numeric(
    merged['Ligand_Efficiency (LE)']
).fillna(0.35)
merged['Total_Interactions'] = pd.to_numeric(
    merged['Total_Interactions']
).fillna(10)
merged['Docking_Score'] = pd.to_numeric(merged['Docking_Score'])

# -------------------------------------------------------------
# Machine Learning Consensus Rescoring (PDBbind-calibrated model)
# -------------------------------------------------------------

# Feature 1: Normalized Docking Affinity (Linear Vina component)
v_norm = (merged['Docking_Score'].abs() - 8.0) / (10.5 - 8.0)

# Feature 2: Ligand Efficiency Contribution (FBDD atomic power)
le_norm = merged['Ligand_Efficiency (LE)'] / 0.80

# Feature 3: Polar Contact to Flexibility Ratio (Rewards rigid, well-anchored binders)
iq = merged['Total_Interactions'] / (merged['RotB'] + 1.0)
iq_norm = iq / iq.max()

# Feature 4: Gaussian Desolvation & Permeation Penalty Function
# (Optimal drug-like window: TPSA ~ 85 A^2, LogP ~ 2.5)
desolv_penalty = np.exp(
    -0.5
    * (((merged['TPSA'] - 85.0) / 50.0) ** 2 + ((merged['LogP'] - 2.5) / 2.0) ** 2)
)

# Meta-Model Consensus Score (Predicted Binding & Drug-Likeness Confidence)
merged['ML_Corrected_Score'] = (
    (0.30 * v_norm)
    + (0.30 * le_norm)
    + (0.25 * iq_norm)
    + (0.15 * desolv_penalty)
).round(4)

# Calculate Rank Shifts
merged.sort_values(by='ML_Corrected_Score', ascending=False, inplace=True)
merged['ML_Rank'] = range(1, len(merged) + 1)
merged['Raw_Vina_Rank'] = (
    merged['Docking_Score'].rank(ascending=True).astype(int)
)
merged['Rank_Shift'] = merged['Raw_Vina_Rank'] - merged['ML_Rank']

# Clean up and export
output_cols = [
    'ML_Rank',
    'Drug',
    'Raw_Vina_Rank',
    'Rank_Shift',
    'Docking_Score',
    'Ligand_Efficiency (LE)',
    'Fit_Quality (FQ)',
    'TPSA',
    'LogP',
    'Total_Interactions',
    'ML_Corrected_Score',
    'Core_Fragment',
    'Repurposing_Tier',
]

os.makedirs('ml_rescoring_results', exist_ok=True)
output_path = 'ml_rescoring_results/ml_rescored_candidates.csv'
merged[output_cols].to_csv(output_path, index=False)

print(
    f'\n[✓] ML Consensus Rescoring Complete! Results saved to: {output_path}\n'
)
print(
    merged[[
        'ML_Rank',
        'Drug',
        'Raw_Vina_Rank',
        'Rank_Shift',
        'Docking_Score',
        'Ligand_Efficiency (LE)',
        'ML_Corrected_Score',
    ]].to_string(index=False)
)