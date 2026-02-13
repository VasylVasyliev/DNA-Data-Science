# 🧬 DNA Data Science & Molecular Dynamics: Bioinformatics Portfolio

This repository represents a comprehensive journey through bioinformatics—from algorithmic Rosalind solutions to AI-powered viral analysis and high-performance Molecular Dynamics (MD).

---

## 🚀 FEATURED 2026: Nano-Interface & Protein Stability
**Goal:** Comparative study of heavy metal (Au, Pt, Ag) interaction stability within cellular protein targets (HeLa/EGFR models).

### 1. High-Precision Silver (Ag) Docking at L858R Mutation
Targeted modeling of a 100-atom Silver cluster interacting with the **EGFR L858R** mutation (PDB: 4I22), a key driver in HeLa cell proliferation.

<p align="center">
  <img src="./reports/docking_energy_plot.png" width="600px">
  <br>
  <i>Figure 1: Targeted docking analysis of metal nanoparticles into the ATP-binding pocket of mutant EGFR.</i>
</p>

* **Target:** EGFR kinase domain (L858R mutation site).
* **Mechanism:** Physical blockade of the active site to inhibit ATP binding and downstream signaling.
* **Status:** Structural alignment completed; coordinates centered at $(9.634, -7.326, 22.414)$.

### 2. Comparative MD Stability (Au vs Pt vs Ag)
Statistical validation of binding energies at physiological temperature (37°C) using the Amber force field.

| Metal Target | Mean Energy ($\langle E \rangle$) | Std Deviation ($\sigma$) | Stability Status |
| :--- | :--- | :--- | :--- |
| **Platinum (Pt)** | **-0.838** | **±0.05** | ✅ Highly Stable |
| **Silver (Ag100)** | **-0.720** | **±0.08** | 🟢 Stable / Good Fit |
| **Gold (Au)** | -0.512 | ±0.15 | ⚠️ Moderate Flux |

* **Key Finding:** Silver (Ag100) demonstrates superior stability over Gold, effectively filling the L858R mutant pocket volume due to optimized atom count.

---

## 🤖 AI-Powered: Viral Intelligence Station
A production-grade tool combining real-time genomic data retrieval with Local AI analysis.

* **Neural Core:** Integrated with **Ollama (Llama 3 / Gemma 2)** for scientific reporting.
* **Main UI:** `app.py` (Streamlit-based Cyber-Station).
* **Tech Stack:** Python 3.11, Streamlit, Plotly, BioPython.

---

## 🧬 Genomic Algorithms & Rosalind
* **ORF Discovery:** Scanning PhiX174 for Open Reading Frames (176 candidates).
* **Rosalind Challenges:** 14+ tasks covering core molecular biology algorithms.

---

## 📂 Project Structure (v1.0)
* `app.py` — **Core Interface v1.0** (Locked Feb 13, 2026).
* `data/EGFR_Ag_Study.pse` — PyMOL session file for silver docking analysis.
* `reports/docking_energy_plot.png` — Final interaction visualization.
* `scripts/python/` — Specialized MD and protein analysis scripts.
* `results/plots/` — Statistical output and XRD comparisons.

---
📊 **Check out my live experiments on [Kaggle](https://www.kaggle.com/vasylvasylievvasyl)**

---
*Last synced: 2026-02-13 (v1.0 Release)*

![Traffic Chart](results/plots/platinum_xrd.png)