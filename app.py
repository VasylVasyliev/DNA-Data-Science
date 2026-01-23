import streamlit as st
import pandas as pd
import plotly.express as px
import subprocess
from Bio import Entrez
from Bio.SeqUtils import gc_fraction
from langchain_ollama import OllamaLLM

# --- 1. CONFIGURATION ---
Entrez.email = "vasj5722814@gmail.com" \
""

# --- 2. DYNAMIC OLLAMA MODEL DISCOVERY ---
def get_installed_ollama_models():
    """Detects models installed on the user's local machine."""
    try:
        result = subprocess.run(['ollama', 'list'], capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')[1:]
        models = [line.split()[0] for line in lines if line]
        return models if models else ["llama3.2:1b"]
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ["llama3.2:1b", "llama3"]

# --- 3. BIOINFORMATICS LOGIC ---
@st.cache_data
def fetch_virus_data(input_strings):
    results = []
    unique_items = list(set([i.strip() for i in input_strings if i.strip()]))
    for item in unique_items:
        try:
            target_id = None
            # Search logic for Name or Accession ID
            if not (item.startswith("NC_") or item.startswith("NM_")):
                search_term = f"{item}[Organism] AND complete genome AND srcdb_refseq[PROP]"
                with Entrez.esearch(db="nucleotide", term=search_term, retmax=1) as h:
                    res = Entrez.read(h)
                    if res["IdList"]: target_id = res["IdList"][0]
            else:
                target_id = item

            if target_id:
                with Entrez.efetch(db="nucleotide", id=target_id, rettype="gb", retmode="xml") as h:
                    records = Entrez.read(h)
                    if records:
                        rec = records[0]
                        seq = rec.get('GBSeq_sequence', '').upper()
                        results.append({
                            "Virus ID": target_id,
                            "Molecule": rec.get('GBSeq_moltype', 'N/A').upper(),
                            "Length (bp)": int(rec.get('GBSeq_length', 0)),
                            "GC %": round(gc_fraction(seq) * 100, 2) if seq else 0,
                            "Description": rec.get('GBSeq_definition', 'Unknown'),
                            "Sequence": seq
                        })
        except Exception as e:
            st.error(f"System Error: {item} -> {e}")
    return pd.DataFrame(results)

# --- 4. CYBERPANK INTERFACE ---
st.set_page_config(page_title="Viral Intelligence Station", layout="wide", page_icon="🧬")

# Title and Visual Identity
st.title("🧬 Viral Intelligence Station")

# Main Instruction Block (Black text on Neon)
with st.expander("📖 USER MANUAL / SYSTEM PROTOCOL", expanded=True):
    st.markdown("""
    <style>
        .cyber-box {
            background: linear-gradient(90deg, #00ffcc 0%, #0077ff 100%);
            color: #000000;
            font-family: 'Segoe UI', sans-serif;
            font-weight: bold;
            padding: 20px;
            border-radius: 5px;
            border: 2px solid #ff00ff;
            line-height: 1.5;
        }
    </style>
    <div class="cyber-box">
        1. AI INITIALIZATION: Ensure <strong>Ollama</strong> is running. Select your model from the "NEURAL CORE" menu.<br>
        2. DATA ACQUISITION: Select viruses from presets or enter custom IDs. Click "EXECUTE ANALYTICS".<br>
        3. SEQUENCE ANALYSIS: Review GC-stability and Genome Capacity charts using dark telemetry templates.<br>
        4. INTELLIGENCE REPORT: Click "GENERATE REPORT" to activate the AI analysis and translation module.<br>
        5. DATA EXPORT: Save your research batch using the "DOWNLOAD DATASET" button in .FASTA format.
    </div>
    """, unsafe_allow_html=True)

# Sidebar Configuration
st.sidebar.header("⚙️ SYSTEM CONFIG")
installed_models = get_installed_ollama_models()
final_model = st.sidebar.selectbox("🤖 NEURAL CORE:", options=installed_models)

virus_presets = {
    "SARS-CoV-2": "NC_045512", 
    "Ebola virus": "NC_002549", 
    "Zika virus": "NC_012532", 
    "HIV-1": "NC_001802", 
    "Influenza A": "NC_007374"
}
selected = st.sidebar.multiselect("DATA SAMPLES:", options=list(virus_presets.keys()), default=["SARS-CoV-2"])
custom = st.sidebar.text_input("CUSTOM TARGETS (ID/Name):")

query_list = [virus_presets[name] for name in selected]
if custom: query_list.extend(custom.split(","))

run_btn = st.sidebar.button("🚀 EXECUTE ANALYTICS")

# --- DATA VISUALIZATION (CYBER NEON) ---
if run_btn or 'data' in st.session_state:
    if run_btn:
        with st.spinner("📡 SCANNING GLOBAL GENE BANK..."):
            st.session_state.data = fetch_virus_data(query_list)
    
    df = st.session_state.data
    if not df.empty:
        st.subheader("📋 RAW DATA LOG")
        st.dataframe(df.drop(columns=['Sequence']), use_container_width=True, hide_index=True)

        st.divider()
        col1, col2 = st.columns(2)
        
        # Neon palette for charts
        cyber_colors = ['#00ffcc', '#ff00ff', '#0077ff', '#ffcc00']

        with col1:
            st.subheader("🧬 GC-STABILITY SCAN")
            fig_gc = px.bar(df, x='Virus ID', y='GC %', color='Molecule', 
                             text='GC %', color_discrete_sequence=cyber_colors,
                             template="plotly_dark")
            st.plotly_chart(fig_gc, use_container_width=True)

        with col2:
            st.subheader("📏 GENOME CAPACITY (BP)")
            fig_len = px.bar(df, x='Virus ID', y='Length (bp)', color='Molecule', 
                              text='Length (bp)', color_discrete_sequence=cyber_colors,
                              template="plotly_dark")
            st.plotly_chart(fig_len, use_container_width=True)

        # AI Analysis Module
        st.divider()
        if st.button("🧠 GENERATE INTELLIGENCE REPORT"):
            try:
                llm = OllamaLLM(model=final_model)
                summary_data = df[['Virus ID', 'Molecule', 'Length (bp)', 'GC %', 'Description']].to_string(index=False)
                prompt = f"""
                Act as a Cyberpunk Bioinformatician. 
                Analyze these viruses and provide a professional report in RUSSIAN. 
                Translate the descriptions into Russian.
                Data: {summary_data}
                """
                with st.spinner(f"CORE {final_model} PROCESSING..."):
                    st.info(llm.invoke(prompt))
            except Exception as e:
                st.error(f"NEURAL CORE CONNECTION FAILURE: {e}")

        # Dataset Export
        fasta_data = "\n".join([f">{r['Virus ID']} {r['Description']}\n{r['Sequence']}" for _, r in df.iterrows()])
        st.download_button("💾 DOWNLOAD DATASET (.FASTA)", data=fasta_data, file_name="biocode_export.fasta")
else:
    st.info("System Standby. Please select data samples and click 'EXECUTE ANALYTICS'.")