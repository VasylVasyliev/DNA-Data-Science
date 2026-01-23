import streamlit as st
import pandas as pd
from Bio import Entrez
import os

# --- 1. НАСТРОЙКИ ---
# Укажи свой настоящий email, чтобы NCBI не блокировал запросы
Entrez.email = "your_email@example.com" 
st.set_page_config(page_title="Bio-Analyzer Pro", layout="wide", page_icon="🧬")

# --- 2. ЯДРО АНАЛИЗА ---
def analyze_viral_dna(accession_list):
    results = []
    for acc in accession_list:
        acc = acc.strip()
        if not acc: continue
        try:
            # Запрос к базе данных нуклеотидов
            with Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="xml") as handle:
                records = Entrez.read(handle)
            
            if records:
                features = records[0].get('GBSeq_feature-table', [])
                # Считаем количество белок-кодирующих последовательностей (CDS)
                protein_count = sum(1 for f in features if f.get('GBFeature_key') == 'CDS')
                
                results.append({
                    "Virus": acc,
                    "Length": int(records[0].get('GBSeq_length', 0)),
                    "Protein_Count": int(protein_count)
                })
        except Exception as e:
            # Если один ID не сработал, программа выведет ошибку и пойдет дальше
            st.error(f"Ошибка с ID {acc}: {e}")
            
    return pd.DataFrame(results)

# --- 3. ИНТЕРФЕЙС ---
st.title("🧬 Viral Genome Analysis Platform")
st.markdown("---")

st.sidebar.header("🕹️ Control Panel")
accession_input = st.sidebar.text_area(
    "Enter NCBI Accession IDs:", 
    value="NC_045512\nNC_012532\nNC_001422\nNC_001416"
)

if st.sidebar.button("🚀 Start Full Analysis"):
    ids = [i.strip() for i in accession_input.split('\n') if i.strip()]
    
    with st.status("🧬 Analyzing genomes...", expanded=True) as status:
        df = analyze_viral_dna(ids)
        status.update(label="✅ Analysis Complete!", state="complete", expanded=False)

    if not df.empty:
        # УБИРАЕМ ИНДЕКС (0, 1, 2...): делаем Virus главным столбцом
        df_display = df.set_index('Virus')
        
        st.subheader("📊 Protein Count Comparison")
        # Строим график по очищенным данным
        st.bar_chart(df_display[['Protein_Count']])
        
        st.subheader("📋 Genomic Summary Table")
        # Выводим таблицу без первого столбца с нулями
        st.dataframe(df_display, use_container_width=True)

        # Подготовка данных для скачивания
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("📂 Download Results CSV", csv, "virus_data.csv", "text/csv")
    else:
        st.error("No data found. Please check your IDs.")
else:
    st.info("👈 Enter NCBI Accession IDs and click 'Start'.")