import os
import py3Dmol

complex_targets = {
    "1_High_Binding_Midostaurin": "complexes_pdb/Midostaurin_complex.pdb",
    "2_Mid_Binding_Moxifloxacin": "complexes_pdb/Moxifloxacin_complex.pdb",
    "3_Low_Binding_Canagliflozin": "complexes_pdb/Canagliflozin_complex.pdb",
}

os.makedirs("interactive_html_views", exist_ok=True)

for name, pdb_file in complex_targets.items():
  if not os.path.exists(pdb_file):
    print(f"[!] Warning: {pdb_file} not found.")
    continue

  with open(pdb_file, "r") as f:
    pdb_str = f.read()

  viewer = py3Dmol.view(width=900, height=600)
  viewer.addModel(pdb_str, "pdb")

  # Protein Cartoon (Blue)
  viewer.setStyle(
      {"not": {"resn": ["LIG", "UNL"]}},
      {"cartoon": {"color": "#3b82f6", "opacity": 0.85}},
  )

  # Binding Site Residues within 4.5 A (Cyan sticks)
  viewer.addStyle(
      {"byres": True, "within": {"distance": 4.5, "sel": {"resn": "LIG"}}},
      {"stick": {"colorscheme": "cyanCarbon", "radius": 0.16}},
  )

  # Ligand Molecule (Magenta sticks)
  viewer.addStyle(
      {"resn": "LIG"},
      {"stick": {"colorscheme": "magentaCarbon", "radius": 0.32}},
  )

  # Pocket Surface
  viewer.addSurface(
      py3Dmol.VDW,
      {"opacity": 0.18, "color": "#f1f5f9"},
      {"byres": True, "within": {"distance": 5.0, "sel": {"resn": "LIG"}}},
  )

  viewer.zoomTo({"resn": "LIG"})
  viewer.setBackgroundColor("#0b0f19")

  out_path = f"interactive_html_views/{name}.html"
  with open(out_path, "w") as out:
    out.write(viewer._make_html())

  print(f"[✓] Saved interactive viewer: {out_path}")