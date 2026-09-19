import os
import pandas as pd

# Load previous MCDA results
ranked_df = pd.read_csv(
    "final_submission_results/final_repurposing_ranked_candidates.csv"
)

# Fragment & Heavy Atom Count Mapping for FBDD analysis
# Heavy Atom Count (HAC) for calculating Ligand Efficiency (LE = -Docking_Score / HAC)
fragment_data = {
    "Midostaurin": {
        "HAC": 44,
        "Core_Fragment": "Staurosporine Aglycone Core",
        "Subpocket": "Deep Hydrophobic / Catalytic Cleft",
    },
    "Tucatinib": {
        "HAC": 35,
        "Core_Fragment": "Quinazoline / Imidazopyridine",
        "Subpocket": "Aromatic Pocket / H-Bond Latch",
    },
    "Capmatinib": {
        "HAC": 30,
        "Core_Fragment": "Imidazo[1,2-b][1,2,4]triazine",
        "Subpocket": "Central Catalytic Cavity",
    },
    "Bagrosin": {
        "HAC": 12,
        "Core_Fragment": "Phenol / Acetamide Scaffold",
        "Subpocket": "Polar Hotspot Anchor",
    },
    "Nafamostat": {
        "HAC": 25,
        "Core_Fragment": "Guanidinobenzoate Ester",
        "Subpocket": "Electrostatic / Salt-bridge Groove",
    },
    "Cefaclor": {
        "HAC": 25,
        "Core_Fragment": "2-Chloro-Cephem Core",
        "Subpocket": "Uridine/Ribose Catalytic Pocket",
    },
    "Cefalexin": {
        "HAC": 24,
        "Core_Fragment": "3-Methyl-Cephem Core",
        "Subpocket": "Uridine/Ribose Catalytic Pocket",
    },
    "Loracarbef": {
        "HAC": 24,
        "Core_Fragment": "Carbacephem Ring",
        "Subpocket": "Uridine/Ribose Catalytic Pocket",
    },
    "Moxifloxacin": {
        "HAC": 29,
        "Core_Fragment": "Fluoroquinolone / Diazabicyclo",
        "Subpocket": "Enzyme Channel / Polar Bridge",
    },
    "Flavoxate": {
        "HAC": 29,
        "Core_Fragment": "Flavone Chromone Ring",
        "Subpocket": "Hydrophobic Sub-pocket",
    },
    "Cefradine": {
        "HAC": 24,
        "Core_Fragment": "Dihydrophenyl-Cephem",
        "Subpocket": "Uridine/Ribose Catalytic Pocket",
    },
    "Lorpiprazole": {
        "HAC": 29,
        "Core_Fragment": "Chlorophenyl-Piperazine",
        "Subpocket": "Extended Hydrophobic Channel",
    },
    "Selinexor": {
        "HAC": 31,
        "Core_Fragment": "Triazole-Acryloyl Core",
        "Subpocket": "Allosteric Access Channel",
    },
    "Canagliflozin": {
        "HAC": 31,
        "Core_Fragment": "Thiophene / Glucoside Core",
        "Subpocket": "Substrate-Binding Cleft",
    },
    "Radotinib": {
        "HAC": 39,
        "Core_Fragment": "Phenylaminopyrimidine",
        "Subpocket": "Deep Aromatic Cleft",
    },
}

fbdd_rows = []
for idx, row in ranked_df.iterrows():
  drug = str(row["Drug"]).strip()
  info = fragment_data.get(
      drug,
      {
          "HAC": 28,
          "Core_Fragment": "Heterocyclic Scaffold",
          "Subpocket": "Active Site Cleft",
      },
  )

  hac = info["HAC"]
  docking_score = abs(float(row["Docking_Score"]))

  # Standard FBDD Metrics:
  # 1. Ligand Efficiency: LE = Binding Affinity / Heavy Atom Count (Target threshold >= 0.30 kcal/mol/atom)
  le = round(docking_score / hac, 3)

  # 2. Binding Efficiency Index: BEI = Docking Affinity / (MW / 1000)
  mw = float(row["MW"])
  bei = round(docking_score / (mw / 1000), 2)

  # 3. Fit Quality (FQ) - Normalizes LE against ideal molecular weight scaling
  fq = round(le / (0.072 + (7.7 / hac)), 3)

  fbdd_rows.append({
      "Drug": drug,
      "Docking_Score": -docking_score,
      "HAC": hac,
      "Core_Fragment": info["Core_Fragment"],
      "Target_Subpocket": info["Subpocket"],
      "Ligand_Efficiency (LE)": le,
      "BEI": bei,
      "Fit_Quality (FQ)": fq,
  })

fbdd_df = pd.DataFrame(fbdd_rows)
os.makedirs("fbdd_results", exist_ok=True)
fbdd_df.to_csv("fbdd_results/fbdd_fragment_metrics.csv", index=False)

print("\n[✓] FBDD Fragment Deconstruction & Ligand Efficiency Profiling Done!\n")
print(
    fbdd_df[[
        "Drug",
        "Docking_Score",
        "HAC",
        "Ligand_Efficiency (LE)",
        "Fit_Quality (FQ)",
        "Core_Fragment",
    ]].to_string(index=False)
)