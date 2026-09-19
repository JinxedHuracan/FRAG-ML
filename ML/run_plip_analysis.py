import glob
import os
import pandas as pd
from plip.structure.preparation import PDBComplex

complex_files = sorted(glob.glob("complexes_pdb/*_complex.pdb"))
interaction_data = []

if not complex_files:
  raise FileNotFoundError(
      "No complex files found in complexes_pdb/. Run build_complexes_mac.py"
      " first."
  )

os.makedirs("interaction_tables", exist_ok=True)
print(f"[*] Found {len(complex_files)} complexes to profile...")

for complex_file in complex_files:
  drug_name = (
      os.path.basename(complex_file)
      .replace("_complex.pdb", "")
      .replace("_refined", "")
  )
  print(f" -> Profiling: {drug_name}")

  try:
    mol = PDBComplex()
    mol.load_pdb(complex_file)
    mol.analyze()

    for site_key, site in mol.interaction_sets.items():
      # 1. Hydrophobic Contacts
      for c in site.hydrophobic_contacts:
        interaction_data.append({
            "Drug": drug_name,
            "Residue": f"{c.restype}{c.resnr}",
            "Type": "Hydrophobic",
            "Distance_A": round(c.distance, 2),
            "Details": f"Dist: {round(c.distance, 2)} A",
        })

      # 2. Hydrogen Bonds (Protein donor & Ligand donor)
      for hb in site.hbonds_pdon:
        interaction_data.append({
            "Drug": drug_name,
            "Residue": f"{hb.restype}{hb.resnr}",
            "Type": "H-Bond (Prot Donor)",
            "Distance_A": round(hb.distance_ad, 2),
            "Details": f"{hb.dtype} -> {hb.atype} (Angle: {round(hb.angle, 1)})",
        })
      for hb in site.hbonds_ldon:
        interaction_data.append({
            "Drug": drug_name,
            "Residue": f"{hb.restype}{hb.resnr}",
            "Type": "H-Bond (Lig Donor)",
            "Distance_A": round(hb.distance_ad, 2),
            "Details": f"{hb.dtype} -> {hb.atype} (Angle: {round(hb.angle, 1)})",
        })

      # 3. Pi-Stacking
      for pi in site.pistacking:
        interaction_data.append({
            "Drug": drug_name,
            "Residue": f"{pi.restype}{pi.resnr}",
            "Type": f"Pi-Stacking ({pi.type})",
            "Distance_A": round(pi.distance, 2),
            "Details": f"Angle: {round(pi.angle, 1)} deg",
        })

      # 4. Salt Bridges
      for sb in site.saltbridge_pneg + site.saltbridge_lneg:
        interaction_data.append({
            "Drug": drug_name,
            "Residue": f"{sb.restype}{sb.resnr}",
            "Type": "Salt Bridge",
            "Distance_A": round(sb.distance, 2),
            "Details": f"Dist: {round(sb.distance, 2)} A",
        })

      # 5. Halogen Bonds
      for xb in site.halogen_bonds:
        interaction_data.append({
            "Drug": drug_name,
            "Residue": f"{xb.restype}{xb.resnr}",
            "Type": "Halogen Bond",
            "Distance_A": round(xb.distance, 2),
            "Details": f"{xb.don_res} -> {xb.acc_res}",
        })

  except Exception as e:
    print(f" [!] Error processing {drug_name}: {e}")

# Save compiled tables
df = pd.DataFrame(interaction_data)

if not df.empty:
  df.to_csv("interaction_tables/detailed_interactions.csv", index=False)

  freq = pd.crosstab(df["Residue"], df["Type"], margins=True)
  freq.sort_values(by="All", ascending=False).to_csv(
      "interaction_tables/target_hotspot_frequency.csv"
  )

  matrix = pd.crosstab(df["Drug"], df["Residue"])
  matrix.to_csv("interaction_tables/drug_residue_binding_matrix.csv")

  print("\n[✓] Done! Files saved in interaction_tables/:")
  print("  • detailed_interactions.csv")
  print("  • target_hotspot_frequency.csv")
  print("  • drug_residue_binding_matrix.csv")
else:
  print("[!] No interactions detected.")