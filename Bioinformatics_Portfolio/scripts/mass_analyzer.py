import os
import pandas as pd
from Bio import Entrez, SeqIO

# Твои данные для NCBI
Entrez.email = "your_email@example.com" 

def count_proteins(sequence):
    """Считает количество потенциальных белков (от ATG до Стоп-кодона)"""
    count = 0
    # Ищем старт-кодон ATG
    for i in range(len(sequence) - 3):
        if sequence[i:i+3] == "ATG":
            # Если нашли старт, ищем ближайший стоп-кодон в той же рамке
            for j in range(i + 3, len(sequence) - 3, 3):
                codon = sequence[j:j+3]
                if codon in ["TAA", "TAG", "TGA"]:
                    if (j - i) > 100: # Берем только белки длиннее 33 аминокислот
                        count += 1
                    break
    return count

def mass_process(input_list_path):
    summary_data = []
    with open(input_list_path, 'r') as f:
        ids = [line.strip() for line in f if line.strip()]

    for accession in ids:
        print(f"🔬 Анализируем {accession}...")
        try:
            with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
                record = SeqIO.read(handle, "fasta")
            
            seq = str(record.seq)
            gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
            
            # Считаем белки
            proteins = count_proteins(seq)
            
            summary_data.append({
                "Accession": accession,
                "Organism": record.description.split(',')[0],
                "Length_bp": len(seq),
                "GC_Content_%": round(gc_content, 2),
                "Proteins_Found": proteins # ВОТ ТА САМАЯ КОЛОНКА
            })
        except Exception as e:
            print(f"⚠️ Ошибка с {accession}: {e}")

    df = pd.DataFrame(summary_data)
    # Сохраняем в результаты
    os.makedirs("results", exist_ok=True)
    df.to_csv("results/virus_report.csv", index=False)
    
    print("\n✅ Новый отчет готов:")
    print(df.to_string(index=False))

if __name__ == "__main__":
    mass_process("virus_list.txt")