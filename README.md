# FRAG-ML

### Fragment-Rescored Affinity Grid with Machine Learning

**FRAG-ML** is an in-silico drug repurposing pipeline designed to identify promising therapeutic candidates against **ML2177c (Uridine Nucleoside Phosphorylase)** from *Mycobacterium leprae*, the causative agent of leprosy.

The project combines **high-throughput molecular docking, Fragment-Based Drug Discovery (FBDD), Protein-Ligand Interaction Profiler (PLIP) analysis, Lipinski Rule-of-5 profiling, and physics-informed machine-learning rescoring** to reduce the heavy-atom bias present in conventional docking-based screening.

> **Core idea:** A molecule should not rank highly simply because it is large and makes many contacts. FRAG-ML evaluates how efficiently and realistically a compound binds.

---

# Authors

**FRAG-ML — Fragment-Rescored Affinity Grid with Machine Learning**

Developed as a computational drug-repurposing project focused on *Mycobacterium leprae* ML2177c.

Contributions by 
Shreyas Nigam (shreyasnigam01).
Kunsh Billa (JinxedHuracan) .


---

## Why FRAG-ML?

Conventional virtual screening often ranks compounds primarily according to their docking score.

This creates a major limitation: **larger molecules can receive favorable docking scores simply because they have more heavy atoms interacting with the target.**

FRAG-ML addresses this by incorporating:

* Docking affinity
* Ligand Efficiency (LE)
* Protein-ligand interaction quality
* Molecular flexibility
* Lipinski drug-likeness
* Desolvation/permeation penalties
* Fragment-based chemical efficiency

This produces a more balanced ranking of potential candidates.

---

# Target

**Protein:** ML2177c
**Organism:** *Mycobacterium leprae*
**Target:** Uridine Nucleoside Phosphorylase (UPase)

The receptor structure was obtained from structural databases including the **RCSB Protein Data Bank and AlphaFold Protein Structure Database**. The receptor was prepared by removing non-essential crystallographic waters, adding missing hydrogens, assigning partial charges, and converting the structure into a docking-compatible `.pdbqt` format.

---

# Pipeline

```text
Target Structure + Drug Library
             │
             ▼
┌─────────────────────────────┐
│ Phase 1: Target Preparation │
│        & Grid Benchmark     │
└──────────────┬──────────────┘
               ▼
┌─────────────────────────────┐
│ Phase 2: High-Throughput    │
│        Virtual Screening    │
└──────────────┬──────────────┘
               ▼
┌─────────────────────────────┐
│ Phase 3: High-Exhaustiveness│
│        Redocking            │
└──────────────┬──────────────┘
               ▼
┌─────────────────────────────┐
│ Phase 4: PLIP Interaction   │
│        Profiling            │
└──────────────┬──────────────┘
               ▼
┌─────────────────────────────┐
│ Phase 5: Lipinski rule      |
|          of 5 verification  │
│                             │
└──────────────┬──────────────┘
               ▼
┌─────────────────────────────┐
│ Phase 6: FBDD + ML Rescoring│
└──────────────┬──────────────┘
               ▼
       Final Candidate Ranking
               │
               ▼
      Streamlit Visualization
```

The documented project data flow follows this six-stage structure, from receptor preparation through final ML rescoring and the Streamlit interface.

---

# Methodology

## 1. Target Preparation & Grid Benchmark

A docking grid was centered on the catalytic active site of ML2177c.

A known reference ligand was redocked to establish a benchmark affinity of approximately:

**−6.72 kcal/mol**

Compounds performing better than this benchmark were considered potential binders.

---

## 2. High-Throughput Virtual Screening

Approximately **2,883 FDA-approved and bioactive compounds** were screened against ML2177c using **AutoDock Vina** with multicore processing.

The initial screen identified the **Top 15 candidates**, with docking scores ranging from approximately:

**−8.97 to −10.14 kcal/mol**

The compound library consisted of molecules sourced from open-access repositories including DrugBank, the FDA Orange Book, and ZINC20-approved libraries.

---

## 3. High-Exhaustiveness Redocking

The Top 15 compounds from the initial screen were redocked using higher exhaustiveness values of **32–64**.

This step was used to improve confidence in the predicted binding poses and energy estimates.

---

## 4. Protein-Ligand Interaction Profiling

Protein-ligand complexes were assembled into `.pdb` structures and analyzed using **PLIP**.

The interaction analysis examined:

* Hydrogen bonds
* π-stacking interactions
* Salt bridges
* Residue-level interaction patterns

The resulting interaction data was also used to construct active-site interaction heatmaps.

---

## 5. Lipinski Rule-of-5 verification

Binding affinity alone does not establish drug-likeness.

The shortlisted compounds were therefore evaluated using physicochemical descriptors including:

* Molecular Weight (MW)
* LogP
* Topological Polar Surface Area (TPSA)
* Hydrogen Bond Donors (HBD)
* Hydrogen Bond Acceptors (HBA)
* Rotatable Bonds

Lipinski Rule-of-5 violation flags were calculated using RDKit.

---

## 6. Fragment-Based Drug Discovery & ML Rescoring

FRAG-ML calculates **Ligand Efficiency (LE)** to account for the number of heavy atoms contributing to binding:

```text
Ligand Efficiency = Binding Energy / Heavy Atom Count
```

This helps distinguish genuinely efficient binders from compounds whose docking scores are primarily driven by molecular size.

The final composite rescoring framework considers:

| Component                             | Weight |
| ------------------------------------- | -----: |
| Docking Energy                        |    30% |
| Ligand Efficiency                     |    30% |
| Polar Contact Quality vs. Flexibility |    25% |
| Desolvation & Permeation Penalty      |    15% |

## The final dataset combines docking, FBDD, interaction-density and ADMET-derived features into a **ML_Corrected_Score** and **Rank_Shift**.

# Key Results

One of the major outcomes of FRAG-ML was that the ranking of several compounds changed substantially after rescoring.

| Compound         | Raw Vina Rank | ML-Corrected Rank | Rank Shift |
| ---------------- | ------------: | ----------------: | ---------: |
| **Midostaurin**  |   #1 (−10.14) |                #5 |         −4 |
| **Moxifloxacin** |    #8 (−9.06) |                #4 |         +4 |
| **Bagrosin**     |    #4 (−9.54) |                #1 |         +3 |

### Midostaurin

Midostaurin initially ranked #1 based on raw docking affinity.

However, it contains **44 heavy atoms**, and its ligand efficiency was only **0.23**. FRAG-ML therefore reduced its ranking from #1 to #5, illustrating the effect of correcting for molecular bulk.

### Moxifloxacin

Moxifloxacin moved from **#8 to #4** after rescoring.

Its documented advantages include:

* Ligand Efficiency of **0.313**
* No Lipinski violations in the analyzed profile
* Compact molecular structure
* Favorable membrane-permeation characteristics
* Existing human safety/clinical usage history

### Bagrosin

Bagrosin moved from **#4 to #1**.

With only **12 heavy atoms** and a reported ligand efficiency of **0.795**, it represents a particularly efficient chemical starting point for fragment-based lead development.

---

# Candidate Interpretation

The project identifies two different types of promising candidates:

### Direct Repurposing Candidate

**Moxifloxacin**

The project identifies Moxifloxacin as the strongest candidate for direct drug repurposing based on its combination of binding, efficiency, drug-likeness and existing clinical use.

### Novel Lead Scaffolds

**Bagrosin** and **Tucatinib**

These compounds are proposed as high-efficiency chemical scaffolds that could serve as starting points for future medicinal-chemistry optimization.

> **Important:** These are computationally prioritized candidates, not experimentally validated anti-leprosy drugs. Experimental binding, activity, toxicity and pharmacokinetic validation would still be required.

---

# Project Structure

The documented pipeline produces the following major files and directories:

```text
FRAG-ML/
│
├── 03_docking_grid/
│   └── receptor.pdbqt
│
├── virtual_screening_results.csv
│
├── top15_refined_docking.csv
├── redock_top15_results/
│
├── complexes_pdb/
│   ├── Midostaurin_complex.pdb
│   ├── Moxifloxacin_complex.pdb
│   └── ...
│
├── interaction_tables/
│   ├── detailed_interactions.csv
│   └── drug_residue_binding_matrix.csv
│
├── admet_results/
│   └── top15_admet_profile.csv
│
├── fbdd_results/
│   └── fbdd_fragment_metrics.csv
│
├── ml_rescoring_results/
│   └── ml_rescored_candidates.csv
│
└── app.py
```

These outputs correspond to the documented six-phase data flow and include the intermediate docking, interaction, ADMET, FBDD and final rescoring datasets.

---

# Interactive Dashboard

The project includes a **Streamlit-based interactive dashboard** through `app.py`.

The dashboard is designed to provide:

* Interactive 3D protein-ligand visualization
* PLIP interaction heatmaps
* Candidate comparison
* ML rank-shift visualization
* Exploration of the final candidate rankings

This allows the complete computational workflow to be presented interactively rather than as static docking tables.

---

# Technologies Used

| Category              | Tools                               |
| --------------------- | ----------------------------------- |
| Molecular Docking     | AutoDock Vina                       |
| Molecular Preparation | Open Babel / MGLTools               |
| Molecular Descriptors | RDKit                               |
| Interaction Analysis  | PLIP                                |
| Fragment Analysis     | FBDD / Ligand Efficiency            |
| Data Processing       | Python                              |
| Visualization         | Streamlit                           |
| Structural Data       | RCSB PDB / AlphaFold DB             |
| Compound Libraries    | DrugBank / FDA Orange Book / ZINC20 |

---

# Limitations

FRAG-ML is an **in-silico prioritization framework**. Its predictions should not be interpreted as experimental proof of efficacy.

Important limitations include:

* Molecular docking provides an approximation of binding affinity.
* Docking scores do not directly represent experimentally measured binding free energies.
* Lipinski and calculated ADMET properties are computational estimates.
* ML rescoring is dependent on the selected features, weighting scheme and scoring methodology.
* Candidate ranking requires experimental validation through biochemical and biological assays.

The purpose of the pipeline is therefore to **reduce the search space and prioritize candidates for further investigation**, rather than replace experimental drug development.

---

# Future Work

Potential extensions include:

* Molecular Dynamics (MD) simulations of top-ranked complexes
* MM/GBSA or related binding-energy estimation
* Experimental enzyme inhibition assays
* Cytotoxicity and selectivity testing
* Expanded compound-library screening
* Improved machine-learning models trained on experimentally validated data
* Fragment optimization of high-Ligand-Efficiency scaffolds
* ADMET prediction using additional validated models


---

## Disclaimer

This repository contains computational research and drug-repurposing predictions. The identified compounds are **not presented as clinically validated treatments for leprosy**. Experimental validation is required before any therapeutic conclusions can be drawn.

