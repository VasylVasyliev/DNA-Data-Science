import streamlit as st
import pandas as pd
import plotly.express as px
import ollama
import os
import numpy as np
from Bio import Entrez, SeqIO
from datetime import datetime

# --- 1. CONFIGURATION ---
st.set_page_config(page_title="Bio-Cyber DNA Station", layout="wide", page_icon="🧬")
Entrez.email = "cyber_bio@example.com"

# --- 2. STYLING ---
st.markdown("""
    <style>
    .main { background-color: #0a0a0a; color: #39ff14; }
    .dna-thinking { font-size: 30px; display: inline-block; animation: dna-spin 3s linear infinite; }
    @keyframes dna-spin { 0% { transform: rotateY(0deg); } 100% { transform: rotateY(360deg); } }
    .block-container { padding-top: 1rem; padding-bottom: 1rem; }
    </style>
    """, unsafe_allow_html=True)

# --- 3. FUNCTIONS ---
def fetch_genome(accession_id):
    try:
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records: return None
            rec = records[0]
            seq_str = str(rec.seq).upper()
            L = len(seq_str)
            return {
                "ID": rec.id, 
                "Length": L, 
                "Count A": seq_str.count('A'),
                "Count T": seq_str.count('T'),
                "Count G": seq_str.count('G'),
                "Count C": seq_str.count('C'),
                "GC%": round((seq_str.count('G') + seq_str.count('C')) / L * 100, 2), 
                "Protein Count (ATG)": seq_str.count("ATG")
            }
    except Exception: return None

def smooth_data(data, window=5):
    return pd.Series(data).rolling(window=window, min_periods=1).mean()

# --- 4. SIDEBAR (DYNAMICAL MODELS) ---
with st.sidebar:
    st.markdown('<h2 style="color: #39ff14;">⚙️ SYSTEM CONTROL</h2>', unsafe_allow_html=True)
    try:
        models_info = ollama.list()
        # EN: Get names from your local list / RU: Список имен из вашего локального Ollama
        model_names = [m['name'] for m in models_info['models']]
        if not model_names: model_names = ["llama3:latest"]
    except:
        model_names = ["llama3:latest", "gemma2:9b", "llama3.2:1b"]
    
    selected_model = st.selectbox("🤖 Select AI Model:", model_names, index=0)
    st.divider()
    st.info(f"Model Active: {selected_model}")
    st.caption("OS: macOS | Python 3.11")

# --- 5. MAIN INTERFACE ---
st.title("🧬 Bio-Intelligence Cyber-Station")
t1, t2 = st.tabs(["🦠 Viral Intelligence", "🧪 Protein Dynamics (MD)"])

# --- TAB 1: VIROLOGY ---
with t1:
    col_top_in, col_top_chart = st.columns([1, 1.2])
    
    with col_top_in:
        st.subheader("📥 Genomic Input")
        inp = st.text_area("Accession IDs:", "NC_045512, NC_001422", height=100)
        if st.button("🚀 EXECUTE ANALYSIS"):
            res = [fetch_genome(i.strip()) for i in inp.split(",") if i.strip()]
            st.session_state.v_list = [r for r in res if r]

    with col_top_chart:
        if 'v_list' in st.session_state:
            st.subheader("📊 Visualization")
            df_plot = pd.DataFrame(st.session_state.v_list)
            fig = px.bar(df_plot, x="ID", y="Length", color="ID", template="plotly_dark", height=300)
            st.plotly_chart(fig, use_container_width=True)

    if 'v_list' in st.session_state:
        st.divider()
        st.subheader("📋 Detailed Genomic Results")
        df_final = pd.DataFrame(st.session_state.v_list)
        df_final.index = range(1, len(df_final)+1)
        
        # EN: Narrowing ATG column / RU: Сужение колонки ATG
        st.dataframe(
            df_final, 
            use_container_width=True,
            column_config={
                "Protein Count (ATG)": st.column_config.NumberColumn("ATG", width="small"),
                "ID": st.column_config.TextColumn("Accession ID", width="medium")
            }
        )

# --- TAB 2: MOLECULAR DYNAMICS ---
with t2:
    st.header("🧪 EGFR Binding Stability (MD)")
    cv, cs = st.columns([1.5, 1])
    with cv:
        st.subheader("📈 Smoothed Trajectory (Au, Pt, Ag)")
        t = np.linspace(0, 100, 100)
        # Raw data with noise
        r_gold = -0.5 + 0.1 * np.sin(t/5) + 0.08 * np.random.randn(100)
        r_platinum = -0.8 + 0.05 * np.sin(t/10) + 0.04 * np.random.randn(100)
        r_silver = -0.7 + 0.07 * np.sin(t/7) + 0.05 * np.random.randn(100)
        
        fig_md = px.line(template="plotly_dark")
        fig_md.add_scatter(x=t, y=smooth_data(r_gold), name="Gold (Au)", line=dict(color="#FFD700", width=2))
        fig_md.add_scatter(x=t, y=smooth_data(r_platinum), name="Platinum (Pt)", line=dict(color="#000000", width=2))
        fig_md.add_scatter(x=t, y=smooth_data(r_silver), name="Silver (Ag100)", line=dict(color="#808080", width=3))
        
        fig_md.update_layout(xaxis_title="Time (ps)", yaxis_title="Energy (kcal/mol)")
        st.plotly_chart(fig_md, use_container_width=True)
        
        st.subheader("🧬 Molecular Docking Analysis")
        img_p = "reports/docking_energy_plot.png"
        if os.path.exists(img_p):
            st.image(img_p, use_container_width=True)
        else:
            st.warning("Image 'reports/docking_energy_plot.png' not found.")

    with cs:
        st.subheader("📊 MD Stats")
        md_data = {
            "Target": ["Platinum (Pt)", "Silver (Ag100)", "Gold (Au)"], 
            "Energy": [-0.838, -0.712, -0.512],
            "Status": ["✅ Stable", "🟢 Stable", "⚠️ Moderate"]
        }
        st.table(pd.DataFrame(md_data, index=range(1, 4)))
    
    if st.button("🤖 GENERATE AI REPORT"):
        with st.chat_message("assistant"):
            st.markdown("<div class='dna-thinking'>🧬</div> *Analyzing...*", unsafe_allow_html=True)
            try:
                prompt = f"Data: {md_data}. Explain metal stability for EGFR project."
                resp = ollama.chat(model=selected_model, messages=[{'role': 'user', 'content': prompt}])
                st.markdown(resp['message']['content'])
            except Exception as e:
                st.error(f"AI Connection Error: Ensure {selected_model} is running.")

st.sidebar.caption(f"DNA-Data-Science | {datetime.now().strftime('%Y-%m-%d')}")