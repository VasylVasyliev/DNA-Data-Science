import streamlit as st
import pandas as pd
import plotly.express as px
import ollama
from Bio import Entrez, SeqIO

# --- 1. КОНФИГУРАЦИЯ ---
st.set_page_config(page_title="Bio-Cyber DNA Station", layout="wide", page_icon="🧬")
Entrez.email = "cyber_bio@example.com"

# --- 2. СТИЛИЗАЦИЯ (Cyber-Bio Theme) ---
st.markdown("""
    <style>
    .main { background-color: #0a0a0a; color: #39ff14; }
    .black-header { 
        color: #000000 !important; 
        font-weight: bold;
        margin-top: 15px;
        margin-bottom: 15px;
    }
    @keyframes dna-spin {
        0% { transform: rotateY(0deg); }
        100% { transform: rotateY(360deg); }
    }
    .dna-loader {
        font-size: 60px;
        display: inline-block;
        animation: dna-spin 2s linear infinite;
        transform-style: preserve-3d;
    }
    .dna-thinking {
        font-size: 30px;
        display: inline-block;
        animation: dna-spin 3s linear infinite;
    }
    div[data-testid="stNumberInput"] button { 
        background-color: #1a1a1a; 
        color: #39ff14; 
    }
    </style>
    """, unsafe_allow_html=True)

# --- 3. ФУНКЦИИ ---
def fetch_genome(accession_id):
    try:
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records: return None
            rec = records[0]
            L = len(rec.seq)
            counts = {b: rec.seq.count(b) for b in "ATGC"}
            return {
                "ID": rec.id, "Length": L, 
                "GC%": round((counts['G'] + counts['C']) / L * 100, 2),
                "A%": round(counts['A'] / L * 100, 2), 
                "T%": round(counts['T'] / L * 100, 2),
                "G%": round(counts['G'] / L * 100, 2), 
                "C%": round(counts['C'] / L * 100, 2),
                "Sequence": str(rec.seq) 
            }
    except: return None

def calculate_hamming(seq1, seq2):
    if len(seq1) != len(seq2): return None
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

# --- 4. ПАМЯТЬ СЕССИИ ---
if 'virus_list' not in st.session_state: st.session_state.virus_list = []
if 'comparison_result' not in st.session_state: st.session_state.comparison_result = None

# --- 5. БОКОВАЯ ПАНЕЛЬ (КОНТРОЛЬ МОДЕЛЕЙ) ---
with st.sidebar:
    st.markdown('<h2 class="black-header">⚙️ SYSTEM CONTROL</h2>', unsafe_allow_html=True)
    try:
        models_info = ollama.list()
        model_names = [m['name'] for m in models_info['models']]
        if not model_names:
            model_names = ["llama3:latest", "llama3.2:1b"]
    except: 
        model_names = ["llama3:latest", "llama3.2:1b"]
    
    selected_model = st.selectbox("🤖 Select AI Model:", model_names, index=0)
    st.divider()
    st.success("Platform: MacBook Air 2026")

# --- 6. ОСНОВНОЙ ИНТЕРФЕЙС (TABS) ---
st.title("🧬 Bio-Intelligence Cyber-Station")

tab_viral, tab_md = st.tabs(["🦠 Viral Intelligence", "🧪 Protein Dynamics (MD)"])

# --- ВКЛАДКА 1: ВИРУСОЛОГИЯ ---
with tab_viral:
    col_in, col_chart = st.columns([1, 2])

    with col_in:
        st.subheader("📥 Data Input")
        ids_input = st.text_area("Accession IDs:", "NC_045512, NC_003391", height=100)
        if st.button("🚀 EXECUTE"):
            raw_ids = [i.strip() for i in ids_input.replace(',', ' ').split() if i.strip()]
            st.session_state.virus_list = []
            st.session_state.comparison_result = None 
            with st.status("Sequencing Data...") as status:
                st.markdown("<div style='text-align: center;'><div class='dna-loader'>🧬</div></div>", unsafe_allow_html=True)
                for rid in raw_ids:
                    data = fetch_genome(rid)
                    if data: st.session_state.virus_list.append(data)
                status.update(label="Complete!", state="complete")
            st.rerun()

    with col_chart:
        st.subheader("📊 Comparative View")
        if st.session_state.virus_list:
            df = pd.DataFrame(st.session_state.virus_list)
            fig = px.bar(df, x="ID", y="Length", color="ID", log_y=True, template="plotly_dark", height=250)
            st.plotly_chart(fig, width="stretch") # Обновлено для 2026
            st.dataframe(df.drop(columns=["Sequence"]), width="stretch", hide_index=True)

    st.divider()
    st.markdown('<h2 class="black-header">🔍 Mutation Tracker</h2>', unsafe_allow_html=True)

    if len(st.session_state.virus_list) >= 2:
        virus_options = [v['ID'] for v in st.session_state.virus_list]
        c1, c2, c3, c4 = st.columns([2, 2, 2, 2])
        with c1: ref_id = st.selectbox("Reference:", virus_options, index=0)
        with c2: target_id = st.selectbox("Target:", virus_options, index=1)
        with c3: start_p = st.number_input("Start Position:", min_value=0, value=0, step=1, format="%d")
        with c4: end_p = st.number_input("End Position:", min_value=1, value=100, step=1, format="%d")
        
        if st.button("🧪 COMPARE SEQUENCES"):
            ref_seq = next(v['Sequence'] for v in st.session_state.virus_list if v['ID'] == ref_id)
            tar_seq = next(v['Sequence'] for v in st.session_state.virus_list if v['ID'] == target_id)
            s1, s2 = ref_seq[start_p:end_p], tar_seq[start_p:end_p]
            dist = calculate_hamming(s1, s2)
            if dist is not None:
                st.session_state.comparison_result = {
                    "dist": dist,
                    "diff": "".join([b if a == b else f":red[{b}]" for a, b in zip(s1, s2)]),
                    "range": f"{start_p}:{end_p}",
                    "identity": f"{(len(s1)-dist)/len(s1):.1%}"
                }

        if st.session_state.comparison_result:
            res = st.session_state.comparison_result
            st.info(f"Range {res['range']} | Mutations: {res['dist']} | Identity: {res['identity']}")
            st.markdown(f"**Sequence Map:**\n**{res['diff']}**")

    st.subheader(f"💬 AI Genomic Insights")
    if prompt := st.chat_input("Спроси ИИ об анализе геномов..."):
        with st.chat_message("user"): st.markdown(prompt)
        with st.chat_message("assistant"):
            thinking = st.empty()
            thinking.markdown("<div class='dna-thinking'>🧬</div> *Analyzing...*", unsafe_allow_html=True)
            try:
                clean_ctx = [{k: v for k, v in d.items() if k != 'Sequence'} for d in st.session_state.virus_list]
                resp = ollama.chat(model=selected_model, messages=[
                    {'role': 'system', 'content': 'You are a Bioinformatician. Use Russian.'},
                    {'role': 'user', 'content': f"Data: {clean_ctx}. Analysis: {st.session_state.comparison_result}. Question: {prompt}"}
                ])
                thinking.empty()
                st.markdown(resp['message']['content'])
            except Exception as e:
                thinking.empty()
                st.error(f"AI Error: {e}")

# --- ВКЛАДКА 2: МОЛЕКУЛЯРНАЯ ДИНАМИКА ---
with tab_md:
    st.markdown('<h2 class="black-header">🧪 EGFR Binding Stability (Molecular Dynamics)</h2>', unsafe_allow_html=True)
    
    col_gif, col_stats = st.columns([1, 1])
    
    with col_gif:
        st.subheader("Trajectory Visualization")
        try:
            st.image("./md_stability_animated.gif", width="stretch") # Обновлено для 2026
            st.caption("Simulation: Interaction of EGFR with Metallic Nano-Targets (37°C)")
        except:
            st.warning("GIF-файл не найден. Убедитесь, что он в корневой папке.")

    with col_stats:
        st.subheader("Statistical Validation")
        md_data = {
            "Metal Target": ["Platinum (Pt)", "Gold (Au)"],
            "Mean Energy (⟨E⟩)": [-0.838, -0.512],
            "Std Deviation (σ)": [0.05, 0.15],
            "Stability Status": ["✅ Highly Stable", "⚠️ Moderate Flux"]
        }
        st.table(pd.DataFrame(md_data))
        
        st.info("""
        **Scientific Conclusion:** The Platinum (Pt) complex shows a deeper energy minimum and lower fluctuations. 
        This suggests a more rigid and stable docking compared to Gold (Au).
        """)

    st.divider()
    st.subheader("🤖 AI Structural Analysis")
    if st.button("🧬 GENERATE MOLECULAR REPORT"):
        with st.chat_message("assistant"):
            thinking_md = st.empty()
            thinking_md.markdown("<div class='dna-thinking'>🧬</div> *Analyzing MD Trajectories...*", unsafe_allow_html=True)
            try:
                prompt_md = f"Данные симуляции EGFR: {md_data}. Объясни научно, почему платина лучше золота для терапии."
                resp_md = ollama.chat(model=selected_model, messages=[
                    {'role': 'system', 'content': 'You are a Structural Biologist. Use Russian.'},
                    {'role': 'user', 'content': prompt_md}
                ])
                thinking_md.empty()
                st.markdown(resp_md['message']['content'])
            except Exception as e:
                thinking_md.empty()
                st.error("Убедитесь, что Ollama запущена локально.")