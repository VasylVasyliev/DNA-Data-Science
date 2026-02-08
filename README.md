# 🧬 DNA Data Science & Molecular Dynamics: Bioinformatics Portfolio

This repository represents a comprehensive journey through bioinformatics—from algorithmic Rosalind solutions to AI-powered viral analysis and high-performance Molecular Dynamics (MD).

---

## 🚀 FEATURED 2026: Nano-Interface & Protein Stability
**Goal:** Comparative study of heavy metal (Au vs Pt) interaction stability within cellular protein targets (HeLa/EGFR models).

* **Analysis:** Evaluated dynamic stability at physiological temperature (37°C).
* **Tech Stack:** OpenMM 8.4, Python 3.11, Amber14-all, Matplotlib.

<p align="center">
  <img src="./md_stability_animated.gif?raw=true" width="800px">
  <br>
  <i>Figure 1: Real-time dynamic stability analysis of Gold (Au) and Platinum (Pt) on the EGFR surface.</i>
</p>

### 📊 Statistical Validation
Based on the MD trajectory analysis, Platinum (Pt) demonstrates a significantly more favorable binding profile:

| Metal Target | Mean Energy ($\langle E \rangle$) | Std Deviation ($\sigma$) | Stability Status |
| :--- | :--- | :--- | :--- |
| **Platinum (Pt)** | **-0.838** | **±0.05** | ✅ Highly Stable |
| **Gold (Au)** | -0.512 | ±0.15 | ⚠️ Moderate Flux |

* **Key Finding:** Platinum (Pt) shows lower displacement amplitudes and deeper energy minima, suggesting higher affinity for the EGFR binding site compared to Gold (Au).

---

## 🤖 AI-Powered: Viral Intelligence Station
A production-grade tool combining real-time genomic data retrieval with Local AI analysis.

* **Neural Core:** Integrated with **Ollama (Llama 3)** for scientific reporting.
* **Tech Stack:** Python 3.11, Streamlit, Plotly, LangChain.

---

## 🧬 Genomic Algorithms & Rosalind
* **ORF Discovery:** Scanning PhiX174 for Open Reading Frames (176 candidates).
* **Rosalind Challenges:** 14+ tasks covering core molecular biology algorithms.

---

## 📂 Project Structure
* `main_md.py` — Core Molecular Dynamics simulation script.
* `md_stability_animated.gif` — Animated stability analysis.
* `app.py` — Viral Intelligence dashboard.
* `Bioinformatics_Portfolio/` — Rosalind and legacy genomic scripts.

---
📊 **Check out my live experiments on [Kaggle](https://www.kaggle.com/vasylvasylievvasyl)**

---
*Last synced: 2026-02-08*
