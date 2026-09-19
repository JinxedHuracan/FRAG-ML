import glob
import os
import subprocess

# Auto-locate ml2177c_receptor.pdb in current directory or 01_structure/
if os.path.exists("ml2177c_receptor.pdb"):
  RECEPTOR_PDB = "ml2177c_receptor.pdb"
elif os.path.exists("01_structure/ml2177c_receptor.pdb"):
  RECEPTOR_PDB = "01_structure/ml2177c_receptor.pdb"
else:
  candidates = glob.glob("**/ml2177c_receptor.pdb", recursive=True)
  if candidates:
    RECEPTOR_PDB = candidates[0]
  else:
    raise FileNotFoundError("Could not find ml2177c_receptor.pdb")

INPUT_DIR = "redock_top15_results"
OUTPUT_DIR = "complexes_pdb"

print(f"[*] Found Receptor: {RECEPTOR_PDB}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. Read clean receptor ATOM/HETATM lines (excluding any previous ligands)
with open(RECEPTOR_PDB, "r") as f:
  receptor_lines = [
      line
      for line in f
      if line.startswith(("ATOM", "HETATM"))
      and not any(x in line for x in ["UNL", "LIG"])
  ]

# 2. Locate all refined docked pdbqt files
pdbqt_files = glob.glob(
    os.path.join(INPUT_DIR, "*_refined.pdbqt")
) or glob.glob(os.path.join(INPUT_DIR, "**", "*_refined.pdbqt"), recursive=True)

print(f"[*] Found {len(pdbqt_files)} docked files to assemble...")

for file_path in pdbqt_files:
  base_name = os.path.basename(file_path).replace("_refined.pdbqt", "")
  temp_ligand_pdbqt = os.path.join(OUTPUT_DIR, f"{base_name}_pose1.pdbqt")
  ligand_pdb = os.path.join(OUTPUT_DIR, f"{base_name}_ligand.pdb")
  complex_pdb = os.path.join(OUTPUT_DIR, f"{base_name}_complex.pdb")

  # Extract MODEL 1 (lowest energy docked conformation)
  pose1_lines = []
  recording = False
  with open(file_path, "r") as f:
    for line in f:
      if line.startswith("MODEL 1"):
        recording = True
        continue
      elif line.startswith("ENDMDL"):
        break
      if recording:
        pose1_lines.append(line)

  with open(temp_ligand_pdbqt, "w") as f:
    f.writelines(pose1_lines)

  # Convert Pose 1 to PDB using Open Babel CLI
  subprocess.run(
      ["obabel", "-ipdbqt", temp_ligand_pdbqt, "-opdb", "-O", ligand_pdb],
      check=True,
  )

  # Re-tag ligand residue name to 'LIG A 999' for standard PLIP recognition
  with open(ligand_pdb, "r") as f:
    lig_lines = [
        line[:17] + "LIG A 999" + line[26:]
        for line in f
        if line.startswith(("ATOM", "HETATM"))
    ]

  # Write combined complex file
  with open(complex_pdb, "w") as f:
    f.writelines(receptor_lines)
    f.write("TER\n")
    f.writelines(lig_lines)
    f.write("END\n")

  # Clean up temporary pdbqt file
  if os.path.exists(temp_ligand_pdbqt):
    os.remove(temp_ligand_pdbqt)

  print(f"[✓] Assembled: {complex_pdb}")

print(f"\n[+] Complete! All complexes generated in `{OUTPUT_DIR}/`")