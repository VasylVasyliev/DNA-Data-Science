import streamlit as st
import pandas as pd
import plotly.express as px
import ollama
from Bio import Entrez, SeqIO

# --- 1. КОНФИГУРАЦИЯ ---
st.set_page_config(page_title="Bio-Cyber DNA Station", layout="wide", page_icon="🧬")
Entrez.email = "cyber_bio@example.com"

# --- 2. СТИЛИЗАЦИЯ ---
st.markdown("""
    <style>
    .main { background-color: #0a0a0a; color: #39ff14; }
    
    /* Черные заголовки */
    .black-header { 
        color: #000000 !important; 
        font-weight: bold;
        margin-top: 15px;
        margin-bottom: 15px;
    }
    
    /* Анимация вращения ДНК вокруг своей оси (Y-axis) */
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

    /* Стилизация кнопок +/- */
    div[data-testid="stNumberInput"] button { 
        background-color: #1a1a1a; 
        color: #39ff14; 
    }
    </style>
    """, unsafe_allow_html=True)

# --- 3. ФУНКЦИИ ---
def fetch_genome(accession_id):
    """Загрузка данных и расчет частоты нуклеотидов"""
    try:
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records: return None
            rec = records[0]
            L = len(rec.seq)
            # Полный расчет состава
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

# --- 5. БОКОВАЯ ПАНЕЛЬ (SYSTEM CONTROL) ---
with st.sidebar:
    st.markdown('<h2 class="black-header">⚙️ SYSTEM CONTROL</h2>', unsafe_allow_html=True)
    try:
        models_info = ollama.list()
        model_names = [m['name'] for m in models_info['models']]
    except: model_names = ["llama3.2:1b"]
    
    selected_model = st.selectbox("🤖 Select AI Model:", model_names, index=0)
    st.divider()
    st.success("Platform: Silicon Air 2026")

# --- 6. ОСНОВНОЙ ИНТЕРФЕЙС ---
st.title("🧬 Bio-Intelligence Cyber-Station")

col_in, col_chart = st.columns([1, 2])

with col_in:
    st.subheader("📥 Data Input")
    ids_input = st.text_area("Accession IDs (comma separated):", "NC_045512, NC_003391", height=100)
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
        st.plotly_chart(fig, use_container_width=True)
        # Таблица с A%, T%, G%, C%
        st.dataframe(df.drop(columns=["Sequence"]), use_container_width=True, hide_index=True)

st.divider()

# --- 7. MUTATION TRACKER ---
st.markdown('<h2 class="black-header">🔍 Mutation Tracker</h2>', unsafe_allow_html=True)

if len(st.session_state.virus_list) >= 2:
    virus_options = [v['ID'] for v in st.session_state.virus_list]
    
    # Расширенные колонки (2, 2, 2, 2) чтобы кнопки +/- точно влезли
    c1, c2, c3, c4 = st.columns([2, 2, 2, 2])
    
    with c1: ref_id = st.selectbox("Reference:", virus_options, index=0)
    with c2: target_id = st.selectbox("Target:", virus_options, index=1)
    
    # Параметр format="%d" убирает десятичные точки, step=1 добавляет кнопки
    with c3: start_p = st.number_input("Start Position:", min_value=0, value=0, step=1, format="%d")
    with c4: end_p = st.number_input("End Position:", min_value=1, value=100, step=1, format="%d")
    
    if st.button("🧪 COMPARE"):
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
        else: st.error("Ошибка: Отрезки разной длины!")

    if st.session_state.comparison_result:
        res = st.session_state.comparison_result
        st.info(f"Range {res['range']} | Mutations: {res['dist']} | Identity: {res['identity']}")
        st.markdown(f"**Sequence Map:**\n**{res['diff']}**")

st.divider()

# --- 8. AI CONSULTANT ---
st.subheader(f"💬 AI Genomic Insights ({selected_model})")
if prompt := st.chat_input("Спроси ИИ об анализе..."):
    with st.chat_message("user"): st.markdown(prompt)
    with st.chat_message("assistant"):
        thinking = st.empty()
        thinking.markdown("<div class='dna-thinking'>🧬</div> *Analyzing...*", unsafe_allow_html=True)
        try:
            clean_ctx = [{k: v for k, v in d.items() if k != 'Sequence'} for d in st.session_state.virus_list]
            comp_info = str(st.session_state.comparison_result)
            resp = ollama.chat(model=selected_model, messages=[
                {'role': 'system', 'content': 'You are a Bioinformatician. Use Russian.'},
                {'role': 'user', 'content': f"Data: {clean_ctx}. Analysis: {comp_info}. Question: {prompt}"}
            ])
            thinking.empty()
            st.markdown(resp['message']['content'])
        except Exception as e:
            thinking.empty()
            st.error(f"AI Error: {e}")