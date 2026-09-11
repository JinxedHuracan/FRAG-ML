'''
This code runs each and every compound "one at a time". While Vina uses all CPU threads internally for that one ligand, 
there is a small subprocess startup overhead per compound.
'''
import os
import re
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LIGANDS_DIR = os.path.join(BASE_DIR, "ligands")
RESULTS_DIR = os.path.join(BASE_DIR, "docking_results")
os.makedirs(LIGANDS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

RECEPTOR = os.path.join(BASE_DIR, "ml2177c_receptor.pdbqt")
VINA_EXE = os.path.join(BASE_DIR, "vina_1.2.7_win.exe")

# Frozen validated pocket box
CENTER = [-8.896, 0.493, 4.775]
SIZE = [20.0, 20.0, 20.0]

tsv_path = os.path.join(BASE_DIR, "structures.smiles.tsv")
if not os.path.exists(tsv_path):
    print("Error: structures.smiles.tsv not found in project folder!")
    exit(1)

print("Loading DrugCentral Approved Drug Library...")
df_drugs = pd.read_csv(tsv_path, sep="\t")

# Extract SMILES and Name columns
smiles_col = [c for c in df_drugs.columns if "SMILES" in c.upper()][0]
name_col = [c for c in df_drugs.columns if "NAME" in c.upper() or "DRUG" in c.upper() or "INN" in c.upper()][0]

prepared_ligands = []
print(f"Total entries in library: {len(df_drugs)}. Preparing 3D PDBQT structures...")

# Prepare first batch of small molecules
for idx, row in df_drugs.iterrows():
    raw_name = str(row[name_col]).strip()
    smi = str(row[smiles_col]).strip()
    
    if not smi or smi == "nan":
        continue
        
    clean_name = re.sub(r'[^\w\-_\.]', '_', raw_name)
    pdbqt_path = os.path.join(LIGANDS_DIR, f"{clean_name}.pdbqt")
    
    if not os.path.exists(pdbqt_path):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol_h = Chem.AddHs(mol)
            res = AllChem.EmbedMolecule(mol_h, randomSeed=42)
            if res != 0:
                continue
            AllChem.MMFFOptimizeMolecule(mol_h)
            
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol_h)
            for setup in mol_setups:
                pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(setup)
                if is_ok:
                    with open(pdbqt_path, "w", encoding="utf-8") as f:
                        f.write(pdbqt_string)
                    prepared_ligands.append((clean_name, pdbqt_path))
                    break
        except Exception:
            continue
    else:
        prepared_ligands.append((clean_name, pdbqt_path))

print(f"\n{len(prepared_ligands)} drugs successfully converted to 3D PDBQT.")
print("Starting Virtual Screen with AutoDock Vina...\n")

# --- Run Screening ---
results = []
for i, (name, ligand_file) in enumerate(prepared_ligands):
    out_file = os.path.join(RESULTS_DIR, f"{name}_docked.pdbqt")
    
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
        "--exhaustiveness", "8",
        "--num_modes", "1",
        "--out", out_file
    ]
    
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        match = re.search(r"^\s*1\s+([-\d\.]+)", proc.stdout, re.MULTILINE)
        if match:
            score = float(match.group(1))
            results.append({"Drug_Name": name, "Vina_Score_kcal_mol": score, "Output_File": out_file})
            print(f"[{i+1}/{len(prepared_ligands)}] {name}: {score} kcal/mol")
    except Exception:
        continue

# --- Save & Summary ---
df_res = pd.DataFrame(results)
if not df_res.empty:
    df_res = df_res.sort_values(by="Vina_Score_kcal_mol", ascending=True)
    summary_path = os.path.join(BASE_DIR, "virtual_screening_results.csv")
    df_res.to_csv(summary_path, index=False)
    print("\n" + "="*50)
    print("SCREENING COMPLETE")
    print(f"Full results saved to: {summary_path}")
    print("="*50)
    print("\nTOP CANDIDATES:")
    print(df_res.head(10).to_string(index=False))
