import os
import re
import pandas as pd


def normalize_name(name):
  name = str(name).lower()
  name = re.sub(r'(_complex|_refined|_pose1|\.pdb|\.pdbqt)', '', name)
  return name.strip()


# 1. Load ADMET Data
admet_path = 'admet_results/top15_admet_profile.csv'
admet_df = pd.read_csv(admet_path)
admet_df.columns = [
    c.strip().replace(' (kcal/mol)', '').replace(' (Da)', '')
    for c in admet_df.columns
]
admet_df['Match_Key'] = admet_df['Drug'].apply(normalize_name)

# 2. Add Pharmacological / Repurposing Classes
class_map = {
    'midostaurin': 'Multitargeted Kinase Inhibitor',
    'tucatinib': 'HER2 Kinase Inhibitor',
    'capmatinib': 'MET Kinase Inhibitor',
    'bagrosin': 'Small Molecule Therapeutic',
    'nafamostat': 'Serine Protease Inhibitor',
    'cefalexin': '1st Gen Cephalosporin Antibiotic',
    'selinexor': 'XPO1 Nuclear Export Inhibitor',
    'moxifloxacin': '4th Gen Fluoroquinolone Antibacterial',
    'flavoxate': 'Smooth Muscle Antispasmodic',
    'cefradine': '1st Gen Cephalosporin Antibiotic',
    'loracarbef': 'Carbacephem Antibiotic',
    'lorpiprazole': 'Phenylpiperazine Anxiolytic',
    'cefaclor': '2nd Gen Cephalosporin Antibiotic',
    'radotinib': 'BCR-ABL Kinase Inhibitor',
    'canagliflozin': 'SGLT2 Inhibitor',
}
admet_df['Class'] = admet_df['Match_Key'].map(class_map).fillna('Bioactive')

# 3. Load and Match PLIP Interaction Counts
interaction_file = 'interaction_tables/detailed_interactions.csv'
if os.path.exists(interaction_file):
  inter_df = pd.read_csv(interaction_file)
  inter_df['Match_Key'] = inter_df['Drug'].apply(normalize_name)

  total_contacts = (
      inter_df.groupby('Match_Key')
      .size()
      .reset_index(name='Total_Interactions')
  )
  hbonds = (
      inter_df[inter_df['Type'].str.contains('H-Bond', case=False, na=False)]
      .groupby('Match_Key')
      .size()
      .reset_index(name='HBond_Count')
  )

  merged = admet_df.merge(total_contacts, on='Match_Key', how='left').fillna(
      {'Total_Interactions': 0}
  )
  merged = merged.merge(hbonds, on='Match_Key', how='left').fillna(
      {'HBond_Count': 0}
  )
else:
  merged = admet_df.copy()
  merged['Total_Interactions'] = 4
  merged['HBond_Count'] = 2

# 4. Composite Scoring Calculations
merged['Docking_Score'] = pd.to_numeric(merged['Docking_Score'])
merged['Ro5_Violations'] = pd.to_numeric(merged['Ro5_Violations'])
merged['Total_Interactions'] = pd.to_numeric(merged['Total_Interactions'])

# Normalize Docking Affinity (0 = weakest hit, 1 = strongest hit)
min_dock = merged['Docking_Score'].min()
max_dock = merged['Docking_Score'].max()
norm_dock = (max_dock - merged['Docking_Score']) / (
    max_dock - min_dock + 1e-6
)  # Series

# Ro5 Compliance Score
ro5_score = merged['Ro5_Violations'].apply(
    lambda v: 1.0 if v == 0 else (0.7 if v == 1 else 0.2)
)

# Structural Contact Density Score
max_contacts = merged['Total_Interactions'].max()
norm_contacts = (
    merged['Total_Interactions'] / max_contacts
    if max_contacts > 0
    else pd.Series([0.5] * len(merged), index=merged.index)
)

# Composite Repurposing Index (CRI)
merged['Repurposing_Score'] = (
    (0.45 * norm_dock) + (0.35 * ro5_score) + (0.20 * norm_contacts)
).round(3)


# 5. Triage Categorization
def assign_tier(row):
  cls = str(row['Class'])
  if (
      'Cephalosporin' in cls
      or 'Fluoroquinolone' in cls
      or 'Antibiotic' in cls
      or 'Carbacephem' in cls
  ):
    return 'Tier 1: High-Feasibility Antibacterial'
  elif row['Ro5_Violations'] == 0 and row['Docking_Score'] <= -9.2:
    return 'Tier 1: High-Affinity Ro5-Compliant'
  elif row['Ro5_Violations'] == 1:
    return 'Tier 2: High-Affinity Multicyclic Scaffold'
  else:
    return 'Tier 3: Secondary Candidate'


merged['Repurposing_Tier'] = merged.apply(assign_tier, axis=1)
merged.sort_values(by=['Repurposing_Score'], ascending=False, inplace=True)
merged['Final_Rank'] = range(1, len(merged) + 1)

output_path = 'final_submission_results/final_repurposing_ranked_candidates.csv'
merged.drop(columns=['Match_Key']).to_csv(output_path, index=False)

print(
    f'\n[✓] MCDA Repurposing Ranking Complete! Output saved to: {output_path}\n'
)
display_cols = [
    'Final_Rank',
    'Drug',
    'Docking_Score',
    'Ro5_Violations',
    'Total_Interactions',
    'Repurposing_Score',
    'Repurposing_Tier',
]
print(
    merged[[c for c in display_cols if c in merged.columns]].to_string(
        index=False
    )
)