import os
import re
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LIGANDS_DIR = os.path.join(BASE_DIR, "ligands")
RESULTS_DIR = os.path.join(BASE_DIR, "docking_results")
RECEPTOR = os.path.join(BASE_DIR, "ml2177c_receptor.pdbqt")
VINA_EXE = os.path.join(BASE_DIR, "vina_1.2.7_win.exe")

CENTER = [-8.896, 0.493, 4.775]
SIZE = [20.0, 20.0, 20.0]

def dock_molecule(ligand_name):
    ligand_path = os.path.join(LIGANDS_DIR, ligand_name)
    base_name = ligand_name.replace(".pdbqt", "")
    out_file = os.path.join(RESULTS_DIR, f"{base_name}_docked.pdbqt")
    
    # Skip if already docked previously
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        try:
            with open(out_file, "r") as f:
                content = f.read()
                m = re.search(r"REMARK VINA RESULT:\s+([-\d\.]+)", content)
                if m:
                    return base_name, float(m.group(1)), out_file
        except Exception:
            pass

    cmd = [
        VINA_EXE,
        "--receptor", RECEPTOR,
        "--ligand", ligand_path,
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(SIZE[0]),
        "--size_y", str(SIZE[1]),
        "--size_z", str(SIZE[2]),
        "--cpu", "1",
        "--exhaustiveness", "8",
        "--num_modes", "1",
        "--out", out_file
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
        match = re.search(r"^\s*1\s+([-\d\.]+)", proc.stdout, re.MULTILINE)
        if match:
            return base_name, float(match.group(1)), out_file
    except Exception:
        pass
    return base_name, None, out_file

if __name__ == "__main__":
    ligand_files = [f for f in os.listdir(LIGANDS_DIR) if f.endswith(".pdbqt")]
    total = len(ligand_files)
    print(f"Total ligands to screen: {total}")
    
    # Use 8 parallel workers for your i7
    workers = min(os.cpu_count() or 8, 8)
    print(f"Running in parallel using {workers} CPU workers...")
    
    results = []
    completed = 0
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(dock_molecule, f): f for f in ligand_files}
        for future in as_completed(futures):
            completed += 1
            name, score, out_file = future.result()
            if score is not None:
                results.append({"Drug_Name": name, "Vina_Score_kcal_mol": score, "Output_File": out_file})
                print(f"[{completed}/{total}] {name}: {score} kcal/mol")
    
    # Save ranked CSV
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(by="Vina_Score_kcal_mol", ascending=True)
        summary_path = os.path.join(BASE_DIR, "virtual_screening_results.csv")
        df.to_csv(summary_path, index=False)
        print("\n" + "="*50)
        print("VIRTUAL SCREENING COMPLETED")
        print(f"Results saved to: {summary_path}")
        print("="*50)
        print("\nTOP 15 CANDIDATES:")
        print(df.head(15).to_string(index=False))