import glob
import os
import matplotlib.pyplot as plt
import pandas as pd
import py3Dmol
import seaborn as sns
import streamlit as st
import streamlit.components.v1 as components

# Page configuration
st.set_page_config(
    page_title="SIH 2026 | Leprosy Drug Repurposing Engine",
    page_icon="🔬",
    layout="wide",
)

st.title(
    "🔬 Hybrid FRAG–ML Drug Repurposing Engine for Neglected Disease - Leprosy (SIH 2026)"
)
st.markdown(
    "**Target:** *Mycobacterium leprae* ML2177c (Uridine Phosphorylase) |"
    " **Pipeline:** High-Throughput Screening → FBDD Hotspot Deconstruction →"
    " PLIP Profiling → Physics-Informed ML Rescoring"
)
st.divider()

# Sidebar Navigation
st.sidebar.header("🕹️ Pipeline Navigation")
view_mode = st.sidebar.radio(
    "Select Workflow Module to Inspect:",
    [
        "1. Master Summary & Final Triage",
        "2. 3D Pocket & Ligand Pose Viewer",
        "3. PLIP Interactions & Contact Heatmap",
        "4. FBDD & ML Rescoring Analytics",
        "5. Lipinski Ro5 Profiler",
        "6. Raw Screening & Redocking Datasets",
    ],
)


# Helper function to load and standardize CSV columns
@st.cache_data
def load_clean_csv(path):
  if os.path.exists(path):
    df = pd.read_csv(path)
    # Strip whitespace and common unit annotations
    df.columns = [
        str(c)
        .strip()
        .replace(" (kcal/mol)", "")
        .replace(" (Da)", "")
        .replace(" (A2)", "")
        .replace(" (Å²)", "")
        for c in df.columns
    ]
    return df
  return None


# Paths
ml_path = "ml_rescoring_results/ml_rescored_candidates.csv"
mcda_path = "final_submission_results/final_repurposing_ranked_candidates.csv"
admet_path = "admet_results/top15_admet_profile.csv"
fbdd_path = "fbdd_results/fbdd_fragment_metrics.csv"
inter_path = "interaction_tables/detailed_interactions.csv"
matrix_path = "interaction_tables/drug_residue_binding_matrix.csv"
redock_path = "top15_refined_docking.csv"
screen_path = "virtual_screening_results.csv"

# -------------------------------------------------------------
# 1. MASTER SUMMARY & FINAL TRIAGE
# -------------------------------------------------------------
if view_mode == "1. Master Summary & Final Triage":
  st.subheader("📊 Multi-Criteria Repurposing Prioritization")

  col1, col2, col3, col4 = st.columns(4)
  col1.metric("Screened Library", "3,982 SMILES")
  col2.metric("Control Affinity Benchmark", "-6.72 kcal/mol")
  col3.metric("Top Hit Affinity", "-10.14 kcal/mol", "Δ = -3.42")
  col4.metric("Top Clinical Pick", "Moxifloxacin", "ML Rank #4 (Climbed +4)")

  df_ml = load_clean_csv(ml_path)
  if df_ml is not None:
    st.write("### Final Prioritized Candidate Shortlist")
    st.dataframe(df_ml, use_container_width=True)
    st.download_button(
        "📥 Download Final Shortlist CSV",
        df_ml.to_csv(index=False),
        "final_repurposing_ranked_candidates.csv",
        "text/csv",
    )
  else:
    st.warning("Master ranked results not found.")

# -------------------------------------------------------------
# 2. 3D POCKET & LIGAND POSE VIEWER
# -------------------------------------------------------------
elif view_mode == "2. 3D Pocket & Ligand Pose Viewer":
  st.subheader("🧬 Interactive 3D Target-Ligand Complex Viewer")

  available_files = sorted(glob.glob("complexes_pdb/*_complex.pdb"))

  if not available_files:
    st.error("No complex PDB files found in `complexes_pdb/`.")
  else:
    options_map = {
        os.path.basename(f)
        .replace("_complex.pdb", "")
        .replace("_refined", ""): f
        for f in available_files
    }

    col_sel, col_desc = st.columns([1, 2])
    with col_sel:
      selected_drug = st.selectbox(
          "Choose Candidate Complex:", list(options_map.keys())
      )
      pdb_path = options_map[selected_drug]

      with open(pdb_path, "r") as f:
        pdb_data = f.read()

      st.info(f"**Loaded File:** `{pdb_path}`")
      st.caption(
          "• Blue Cartoon = ML2177c Backbone\n• Cyan Sticks = Active Pocket"
          " Residues (≤ 4.5 Å)\n• Magenta Sticks = Docked Repurposed Candidate"
      )

    with col_desc:
      view = py3Dmol.view(width=750, height=500)
      view.addModel(pdb_data, "pdb")
      view.setStyle(
          {"not": {"resn": ["LIG", "UNL"]}},
          {"cartoon": {"color": "#3b82f6", "opacity": 0.85}},
      )
      view.addStyle(
          {
              "byres": True,
              "within": {"distance": 4.5, "sel": {"resn": "LIG"}},
          },
          {"stick": {"colorscheme": "cyanCarbon", "radius": 0.15}},
      )
      view.addStyle(
          {"resn": "LIG"},
          {"stick": {"colorscheme": "magentaCarbon", "radius": 0.30}},
      )
      view.addSurface(
          py3Dmol.VDW,
          {"opacity": 0.18, "color": "#cbd5e1"},
          {"byres": True, "within": {"distance": 5.0, "sel": {"resn": "LIG"}}},
      )
      view.zoomTo({"resn": "LIG"})
      view.setBackgroundColor("#090d16")

      html_str = f"<!-- {selected_drug} -->\n" + view._make_html()
      components.html(html_str, width=750, height=500)


# -------------------------------------------------------------
# 3. PLIP INTERACTIONS & CONTACT HEATMAP
# -------------------------------------------------------------
elif view_mode == "3. PLIP Interactions & Contact Heatmap":
  st.subheader("🔥 Target Interaction Profiling & Heatmap")

  df_mat = load_clean_csv(matrix_path)
  if df_mat is not None:
    df_mat.set_index(df_mat.columns[0], inplace=True)
    df_mat.index = [str(i).replace("_refined", "") for i in df_mat.index]
    df_mat = df_mat.loc[:, (df_mat != 0).any(axis=0)]

    fig5, ax5 = plt.subplots(figsize=(12, 5.5))
    sns.heatmap(
        df_mat,
        cmap="YlGnBu",
        annot=True,
        fmt="d",
        linewidths=0.5,
        ax=ax5,
        cbar_kws={"label": "Non-Covalent Contacts"},
    )
    plt.title("Active Pocket Residue Fingerprint Across Top 15 Candidates")
    plt.xlabel("ML2177c Target Residues")
    plt.ylabel("Candidates")
    plt.xticks(rotation=45, ha="right")
    st.pyplot(fig5)

  df_inter = load_clean_csv(inter_path)
  if df_inter is not None:
    st.write("### Detailed Non-Covalent Bond Log (PLIP Output)")
    st.dataframe(df_inter, use_container_width=True)


# -------------------------------------------------------------
# 4. FBDD & ML RESCORING ANALYTICS
# -------------------------------------------------------------
elif view_mode == "4. FBDD & ML Rescoring Analytics":
  st.subheader("⚖️ Fragment Deconstruction & Machine Learning Rescoring")
  st.markdown(
      "**Core Innovation:** Classical docking over-rewards heavy molecular"
      " weight. Our **Physics-Informed ML Rescorer** penalizes molecular bulk"
      " and rewards high **Ligand Efficiency (LE)** and **Polar Anchoring**."
  )

  df_ml = load_clean_csv(ml_path)
  df_fbdd = load_clean_csv(fbdd_path)

  if df_ml is not None:
    col_a, col_b = st.columns(2)
    with col_a:
      st.write("#### Ligand Efficiency ($LE \\ge 0.30$)")
      fig1, ax1 = plt.subplots(figsize=(8, 4.5))
      sns.barplot(
          data=df_ml,
          x="Drug",
          y="Ligand_Efficiency (LE)",
          palette="Blues_r",
          ax=ax1,
      )
      plt.axhline(0.30, color="red", linestyle="--", label="Target LE Threshold")
      plt.xticks(rotation=45, ha="right")
      plt.ylabel("LE (kcal/mol/atom)")
      plt.legend()
      st.pyplot(fig1)

    with col_b:
      st.write("#### Rank Shift: Raw Vina Rank vs. ML Rescored Rank")
      fig2, ax2 = plt.subplots(figsize=(8, 4.5))
      colors = ["#22c55e" if x >= 0 else "#ef4444" for x in df_ml["Rank_Shift"]]
      sns.barplot(data=df_ml, x="Drug", y="Rank_Shift", palette=colors, ax=ax2)
      plt.axhline(0, color="black", linestyle="--")
      plt.xticks(rotation=45, ha="right")
      plt.ylabel("Rank Shift (Positive = Promoted by ML)")
      st.pyplot(fig2)

    st.write("### FBDD Fragment Deconstruction Table")
    if df_fbdd is not None:
      st.dataframe(df_fbdd, use_container_width=True)

# -------------------------------------------------------------
# 5.LIPINSKI RO5 PROFILER
# -------------------------------------------------------------
elif view_mode == "5. Lipinski Ro5 Profiler":
  st.subheader("💊 Physicochemical")

  df_admet = load_clean_csv(admet_path)
  if df_admet is not None:
    st.write("### Lipinski Rule of 5 Compliance & Descriptors")
    st.dataframe(df_admet, use_container_width=True)

    # Ensure numeric conversion
    for col in ["MW", "LogP", "TPSA", "Ro5_Violations"]:
      if col in df_admet.columns:
        df_admet[col] = pd.to_numeric(df_admet[col], errors="coerce")

    col1, col2 = st.columns(2)
    with col1:
      fig3, ax3 = plt.subplots(figsize=(7, 4))
      sns.scatterplot(
          data=df_admet,
          x="MW",
          y="LogP",
          hue="Ro5_Violations",
          palette="Set2",
          s=120,
          ax=ax3,
      )
      plt.axvline(500, color="red", linestyle="--", label="MW Limit (500 Da)")
      plt.axhline(5, color="orange", linestyle="--", label="LogP Limit (5.0)")
      plt.title("Lipinski Chemical Space (MW vs LogP)")
      plt.legend()
      st.pyplot(fig3)

    with col2:
      fig4, ax4 = plt.subplots(figsize=(7, 4))
      sns.barplot(data=df_admet, x="Drug", y="TPSA", palette="viridis", ax=ax4)
      plt.axhline(
          140, color="red", linestyle="--", label="Permeability Cap (140 Å²)"
      )
      plt.xticks(rotation=45, ha="right")
      plt.ylabel("TPSA (Å²)")
      plt.title("Topological Polar Surface Area")
      plt.legend()
      st.pyplot(fig4)
  else:
    st.warning(" data file not found.")


# -------------------------------------------------------------
# 6. RAW SCREENING & REDOCKING DATASETS
# -------------------------------------------------------------
elif view_mode == "6. Raw Screening & Redocking Datasets":
  st.subheader("🗄️ Full Screening Logs & Raw Computational Data")

  tab1, tab2 = st.tabs(
      ["High-Exhaustiveness Redocking", "Initial HTVS Library (~4,000 SMILES)"]
  )

  with tab1:
    df_redock = load_clean_csv(redock_path)
    if df_redock is not None:
      st.dataframe(df_redock, use_container_width=True)

  with tab2:
    df_screen = load_clean_csv(screen_path)
    if df_screen is not None:
      st.dataframe(df_screen, use_container_width=True)
      st.caption(f"Showing all {len(df_screen):,} screened molecules.")