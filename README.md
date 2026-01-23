# 🧬 DNA Data Science: Bioinformatics & Sequence Analysis

Repository of **14+ Rosalind solutions** and high-performance genomic data analysis tools. This portfolio demonstrates the transition from fundamental biological algorithms to optimized production-grade sequence analysis and AI-driven insights.

---

## 🚀 FEATURED: Viral Intelligence Station (AI-Powered)
My latest production-grade tool that combines genomic data retrieval with Local AI analysis.

* **Cyberpunk UI:** Interactive dashboard for viral telemetry and sequence exploration.
* **Neural Core:** Integrated with **Ollama** (Llama 3) for automated scientific reporting in Russian.
* **Live Data:** Fetches real-time virus sequences from NCBI Entrez database.
* **Tech Stack:** Python 3.11, Streamlit, Plotly, LangChain.

---

## 🛠 Tech Stack & Methodology
* **Languages:** Python 3.11/3.14 (Latest Features), R (ggplot2).
* **Key Libraries:** Biopython, Ollama, Plotly, NetworkX, Matplotlib.
* **Methodology:** Algorithm optimization (**$O(N)$ complexity**), neighborhood search, and memory caching.

---

## 🧬 Featured Analysis: ORF Discovery & PhiX174 Case Study
Identification of **Open Reading Frames (ORFs)** is a critical step in gene prediction. This tool scans both strands of DNA for start and stop codons to find potential protein-coding sequences.

* **PhiX174 Case Study:** Analyzed the first-ever sequenced DNA using a custom Python detector.
* **Biological Insight:** The analysis identified **176 ORF candidates**. Visualization shows a clear separation between short random sequences and real functional genes (generally >100 amino acids).

---

## 🔍 Ultra-Fast Genomic Marker Discovery
A specialized tool designed to identify unique regulatory motifs or genetic markers across multiple genomes.
* **Fuzzy Pattern Matching:** Finds markers even with mutations ($max\_err=1$) using a Neighborhood Search algorithm.
* **Performance:** Reduced processing time for the PhiX174 genome from 10s to **<0.1s**.
* **Evolution:** Progressed from Brute Force ($10.0s$) to Hash Indexing ($2.0s$) and finally Neighborhood Search (**Instant**).

---

## ✅ Rosalind Challenges (14+ Tasks)
*Located in the `/Bioinformatics_Portfolio` directory.*

### 🧬 Molecular Biology & Translation
* **INI, DNA, RNA, REVC:** Core nucleotide manipulation and transcription/translation.
* **PROT:** Protein translation simulation.
* **CONS:** Profile matrices and consensus sequence discovery.
* **ORF:** Open Reading Frame detection across 6 frames.

### 📊 Genetics, Evolution & Probability
* **GC:** Analysis of DNA thermal stability via Guanine-Cytosine content.
* **HAMM:** Evolutionary distance measurement via point mutations.
* **IPRB, IEV:** Mendelian inheritance and genotype distribution modeling.
* **FIB, FIBD:** Population dynamics using recursive relations and mortality models.

---

## 📂 Project Structure
* `app.py`: Main Viral Intelligence dashboard.
* `Bioinformatics_Portfolio/`: Legacy Windows projects, Rosalind scripts, and R visualizations.
* `analyzer.py`: Core genomic logic for the AI Station.
* `bio_agent.py`: LangChain-powered AI analyst.