import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load the drug-residue binding matrix
matrix_path = "interaction_tables/drug_residue_binding_matrix.csv"
df_matrix = pd.read_csv(matrix_path, index_col=0)

# Clean up index names (remove _refined / file extensions)
df_matrix.index = [idx.replace("_refined", "") for idx in df_matrix.index]

# Filter out residues with 0 interactions across all compounds
df_matrix = df_matrix.loc[:, (df_matrix != 0).any(axis=0)]

plt.figure(figsize=(14, 8))
sns.set_theme(style="white")

# Plot heatmap
ax = sns.heatmap(
    df_matrix,
    cmap="YlGnBu",
    annot=True,
    fmt="d",
    linewidths=0.5,
    cbar_kws={"label": "Number of Non-Covalent Interactions"},
)

plt.title(
    "ML2177c Target Binding Pocket: Drug-Residue Interaction Fingerprint",
    fontsize=14,
    fontweight="bold",
    pad=15,
)
plt.xlabel("Binding Pocket Residues", fontsize=12, fontweight="bold")
plt.ylabel("Top 15 Candidate Drugs", fontsize=12, fontweight="bold")
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(rotation=0, fontsize=10)
plt.tight_layout()

output_png = "interaction_tables/binding_interaction_heatmap.png"
plt.savefig(output_png, dpi=300)
print(f"[✓] Heatmap saved to {output_png}")