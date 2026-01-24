import streamlit as st
import pandas as pd
import plotly.express as px
import ollama
from Bio import Entrez, SeqIO
import time

# 1. Настройка страницы
st.set_page_config(page_title="Bio-Cyber DNA Station", layout="wide")
Entrez.email = "cyber_bio@example.com"

# Кастомный CSS для стиля, таблицы и анимации
st.markdown("""
    <style>
    .main { background-color: #0a0a0a; color: #39ff14; }
    div[data-testid="stDataFrame"] td { white-space: nowrap !important; }
    h2 { color: #00f3ff !important; text-shadow: 0 0 10px #00f3ff; }
    
    @keyframes dna-rotate {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
    }
    .dna-loader {
        font-size: 60px;
        display: inline-block;
        animation: dna-rotate 2s linear infinite;
    }
    /* Медленное вращение для чата */
    .dna-thinking {
        font-size: 30px;
        display: inline-block;
        animation: dna-rotate 4s linear infinite;
    }
    </style>
    """, unsafe_allow_html=True)

def fetch_genome(accession_id):
    try:
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records: return None
            record = records[0]
            seq = record.seq
            L = len(seq)
            a, t, g, c = (seq.count('A')/L*100), (seq.count('T')/L*100), (seq.count('G')/L*100), (seq.count('C')/L*100)
            return {
                "ID": record.id, "Length": L, "GC%": round(g+c, 1),
                "A%": round(a, 1), "T%": round(t, 1), "G%": round(g, 1), "C%": round(c, 1)
            }
    except Exception:
        return None

if 'virus_list' not in st.session_state:
    st.session_state.virus_list = []

st.title("🧬 Bio-Intelligence Cyber-Station")

col_ctrl, col_main = st.columns([1, 2.5])

with col_ctrl:
    st.header("⚙️ CONTROL")
    try:
        models = [m['name'] for m in ollama.list()['models']]
    except: models = ["llama3:latest", "llama3.2:1b"]
    model_name = st.selectbox("🤖 NEURAL AGENT:", models)
    
    st.divider()
    st.subheader("📥 DATA INPUT")
    new_ids_input = st.text_area("ID (через пробел/запятую):", placeholder="NC_045512, NC_003391", height=100)
    
    c1, c2 = st.columns(2)
    with c1:
        if st.button("➕ EXECUTE"):
            raw_ids = list(set(new_ids_input.replace(',', ' ').split())) # Убираем дубликаты из ввода
            if raw_ids:
                with st.empty():
                    st.markdown("<div style='text-align: center;'><div class='dna-loader'>🧬</div><br><b style='color:#39ff14'>Sequencing Data...</b></div>", unsafe_allow_html=True)
                    
                    # Очищаем текущий список перед новым расчетом, чтобы избежать дублирования строк
                    st.session_state.virus_list = [] 
                    
                    for acc in raw_ids:
                        data = fetch_genome(acc)
                        if data:
                            st.session_state.virus_list.append(data)
                    time.sleep(0.5)
                st.rerun()
    with c2:
        if st.button("🗑 RESET"):
            st.session_state.virus_list = []
            st.rerun()

with col_main:
    st.header("📊 COMPARATIVE ANALYTICS")
    if st.session_state.virus_list:
        df = pd.DataFrame(st.session_state.virus_list)
        neon_pal = ["#00f3ff", "#ff00ff", "#39ff14", "#ffff00"]
        
        # 1. График длины (Логарифмический)
        fig_len = px.bar(df, x="ID", y="Length", color="ID", log_y=True,
                         title="Genome Size (Log Scale)",
                         color_discrete_sequence=neon_pal, template="plotly_dark", height=300)
        
        fig_len.update_layout(
            plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
            yaxis=dict(showgrid=True, gridcolor='rgba(255,255,255,0.05)', minor_showgrid=False),
            xaxis=dict(showgrid=False)
        )
        st.plotly_chart(fig_len, use_container_width=True)
        
        # 2. График состава
        df_m = df.melt(id_vars=["ID"], value_vars=["A%", "T%", "G%", "C%"])
        fig_comp = px.bar(df_m, x="ID", y="value", color="variable", barmode="group",
                          title="Nucleotide Composition (%)",
                          color_discrete_sequence=neon_pal, template="plotly_dark", height=300)
        
        fig_comp.update_layout(
            plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
            yaxis=dict(showgrid=True, gridcolor='rgba(255,255,255,0.05)', minor_showgrid=False),
            xaxis=dict(showgrid=False)
        )
        st.plotly_chart(fig_comp, use_container_width=True)
        
        # 3. Таблица
        st.subheader("📑 Detailed Sequence Data")
        st.dataframe(
            df, 
            hide_index=True, 
            use_container_width=True,
            column_config={
                "Length": st.column_config.NumberColumn("Size(bp)", format="%d"),
                "GC%": st.column_config.NumberColumn("GC%", format="%.1f%%"),
                "A%": "A%", "T%": "T%", "G%": "G%", "C%": "C%"
            }
        )
    else:
        st.info("Input IDs and press EXECUTE.")

st.divider()

# --- ЧАТ С МЕДЛЕННОЙ ДНК ---
st.subheader(f"💬 CYBER-CONSULTANT ({model_name})")
if prompt := st.chat_input("Ask about genomes..."):
    with st.chat_message("user"): st.markdown(prompt)
    with st.chat_message("assistant"):
        thinking_status = st.empty()
        thinking_status.markdown("<div class='dna-thinking'>🧬</div> *Analyzing sequences...*", unsafe_allow_html=True)
        
        try:
            instruction = "You are a professional Bioinformatician. Compare viruses biologically based on the data. Avoid code."
            data_ctx = f"Data: {str(st.session_state.virus_list)}. "
            resp = ollama.chat(model=model_name, messages=[
                {'role': 'system', 'content': instruction},
                {'role': 'user', 'content': data_ctx + prompt}
            ])
            thinking_status.empty()
            st.markdown(resp['message']['content'])
        except Exception as e:
            thinking_status.empty()
            st.error(f"AI Error: {e}")