import os
import re
import subprocess
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LIGANDS_DIR = os.path.join(BASE_DIR, "ligands")
HIGH_RES_DIR = os.path.join(BASE_DIR, "redock_top15_results")
RECEPTOR = os.path.join(BASE_DIR, "ml2177c_receptor.pdbqt")
VINA_EXE = os.path.join(BASE_DIR, "vina_1.2.7_win.exe")
CSV_INPUT = os.path.join(BASE_DIR, "top_15_candidates.csv")

os.makedirs(HIGH_RES_DIR, exist_ok=True)

CENTER = [-8.896, 0.493, 4.775]
SIZE = [20.0, 20.0, 20.0]

# 1. Load the top 15 compounds from the CSV
df = pd.read_csv(CSV_INPUT)
top15 = df.head(15)

print("=" * 60)
print("STARTING HIGH-PRECISION REDOCKING FOR TOP 15 CANDIDATES")
print("Settings: Exhaustiveness = 32 | Modes = 9 | CPU = All Cores")
print("=" * 60)

final_records = []

for idx, row in top15.iterrows():
    drug_name = str(row["Drug_Name"])
    ligand_file = os.path.join(LIGANDS_DIR, f"{drug_name}.pdbqt")
    out_file = os.path.join(HIGH_RES_DIR, f"{drug_name}_refined.pdbqt")
    log_file = os.path.join(HIGH_RES_DIR, f"{drug_name}_refined.log")

    if not os.path.exists(ligand_file):
        print(f"[-] Missing input file for: {drug_name}")
        continue

    print(f"\n[{len(final_records) + 1}/15] High-precision docking: {drug_name}...")

    cmd = [
        VINA_EXE,
        "--receptor", RECEPTOR,
        "--ligand", ligand_file,
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(SIZE[0]),
        "--size_y", str(SIZE[1]),
        "--size_z", str(SIZE[2]),
        "--exhaustiveness", "32",
        "--num_modes", "9",
        "--energy_range", "3",
        "--out", out_file
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    
    # Save log output
    with open(log_file, "w") as lf:
        lf.write(proc.stdout)

    # Extract mode 1 energy
    match = re.search(r"^\s*1\s+([-\d\.]+)", proc.stdout, re.MULTILINE)
    best_score = float(match.group(1)) if match else None

    if best_score is not None:
        print(f"    --> Refined Score: {best_score} kcal/mol")
        final_records.append({
            "Rank": len(final_records) + 1,
            "Drug_Name": drug_name,
            "Screening_Score_kcal_mol": row["Vina_Score_kcal_mol"],
            "Refined_Score_kcal_mol": best_score,
            "Refined_Pose_File": out_file
        })

# Save refined comparison table
refined_df = pd.DataFrame(final_records).sort_values(by="Refined_Score_kcal_mol")
refined_df["Final_Rank"] = range(1, len(refined_df) + 1)
output_summary = os.path.join(BASE_DIR, "top15_refined_docking.csv")
refined_df.to_csv(output_summary, index=False)

print("\n" + "=" * 60)
print(f"REDOCKING COMPLETE! Results saved to: {output_summary}")
print("=" * 60)
print(refined_df[["Final_Rank", "Drug_Name", "Refined_Score_kcal_mol", "Screening_Score_kcal_mol"]].to_string(index=False))