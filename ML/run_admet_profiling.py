import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED

# Master dictionary of SMILES for the Top 15 candidates
compounds = {
    "Midostaurin": (
        "CN1C=C(C2=CC=CC=C21)C3=C4C(=C(N3)C5=CC=CC=C54)C6=CNC7=CC=CC=C67"
    ),  # Representative core
    "Tucatinib": "CC1=NC=C(C=C1)C2=C(C=C(C=C2)F)NC(=O)C3=CC(=NN3)C4=NC=NC=C4",
    "Capmatinib": "CC1=NN=C(C=C1)C2=CC(=CC=C2)CC3=NC4=C(N3)N=CN=C4C5=CC=C(C=C5)F",
    "Bagrosin": "CC(=O)NC1=CC=C(C=C1)O",
    "Nafamostat": "C1=CC(=CC=C1C(=O)OC2=CC3=C(C=C2)C=C(C=C3)C(=N)N)N=C(N)N",
    "Cefalexin": (
        "CC1=C(N2C(C(C2=O)NC(=O)C(C3=CC=CC=C3)N)SC1)C(=O)O"
    ),  # First-gen Cephalosporin
    "Selinexor": "CC1=NN=C(N1C2=CC=C(C=C2)C(=O)NNC(=O)C3=CN=CC=N3)C(F)(F)F",
    "Moxifloxacin": (
        "COC1=C2C(=CC(=C1N3CC4CCCNC4C3)F)C(=O)C(=CN2C5CC5)C(=O)O"
    ),
    "Flavoxate": "CCOC(=O)CN1CCN(CC1)C(=O)C2=CC3=CC=CC=C3OC2=O",
    "Cefradine": "CC1=C(N2C(C(C2=O)NC(=O)C(C3=CCC=CC3)N)SC1)C(=O)O",
    "Loracarbef": "C1C(C2=C(N1C(=O)C2NC(=O)C(C3=CC=CC=C3)N)C(=O)O)Cl",
    "Lorpiprazole": (
        "C1CN(CCN1CCCCN2C(=O)N3C(=N2)C=CC=N3)C4=CC=CC(=C4)Cl"
    ),
    "Cefaclor": "C1=CC(=CC=C1C(C(=O)NC2C(=O)N3C(C(=O)O)=C(Cl)CSC23)N)O",
    "Radotinib": (
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
    ),
    "Canagliflozin": (
        "CC1=CC=C(C=C1)CC2=C(C=C(C=C2)C3C(C(C(C(O3)CO)O)O)O)SC4=CC=C(C=C4)F"
    ),
}

# Match Refined Docking Scores
docking_scores = {
    "Midostaurin": -10.140,
    "Tucatinib": -9.741,
    "Capmatinib": -9.668,
    "Bagrosin": -9.544,
    "Nafamostat": -9.203,
    "Cefalexin": -9.122,
    "Selinexor": -9.079,
    "Moxifloxacin": -9.063,
    "Flavoxate": -9.053,
    "Cefradine": -9.043,
    "Loracarbef": -9.035,
    "Lorpiprazole": -9.030,
    "Cefaclor": -9.028,
    "Radotinib": -8.977,
    "Canagliflozin": -8.971,
}

os.makedirs("admet_results", exist_ok=True)
admet_data = []

for drug, smiles in compounds.items():
  mol = Chem.MolFromSmiles(smiles)
  if mol:
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    qed_score = QED.qed(mol)

    # Lipinski Ro5 Violations check
    violations = sum([mw > 500, logp > 5.0, hbd > 5, hba > 10])

    admet_data.append({
        "Drug": drug,
        "Docking_Score (kcal/mol)": docking_scores.get(drug, "N/A"),
        "MW (Da)": round(mw, 2),
        "LogP": round(logp, 2),
        "HBD": hbd,
        "HBA": hba,
        "RotB": rotb,
        "TPSA (A2)": round(tpsa, 2),
        "QED": round(qed_score, 3),
        "Ro5_Violations": violations,
        "Lipinski_Pass": "Yes" if violations <= 1 else "No",
    })

df_admet = pd.DataFrame(admet_data)
df_admet.sort_values(by="Docking_Score (kcal/mol)", inplace=True)
df_admet.to_csv("admet_results/top15_admet_profile.csv", index=False)

print("[✓] ADMET profiling finished! Results saved to:")
print("    admet_results/top15_admet_profile.csv\n")
print(
    df_admet[[
        "Drug",
        "Docking_Score (kcal/mol)",
        "MW (Da)",
        "LogP",
        "Ro5_Violations",
        "Lipinski_Pass",
    ]].to_string(index=False)
)