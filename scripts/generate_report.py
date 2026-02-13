import pandas as pd
import re

# EN: Calculate GC content / RU: Расчет GC-состава
def calculate_gc_content(sequence):
    if not sequence: return 0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

# EN: Estimate proteins using start codon search / RU: Оценка белков через поиск старт-кодонов
def estimate_protein_count(sequence):
    return len(re.findall("ATG", sequence.upper()))

# EN: Final table generator with four pillars: ID, Length, GC, Proteins
# RU: Генератор финальной таблицы с четырьмя столпами: ID, Длина, GC, Белки
def generate_biostat_table(accession_list, data_map):
    results_data = []
    columns = ['ID', 'Length', 'GC%', 'Protein Count']
    for acc in accession_list:
        seq = data_map.get(acc, "")
        results_data.append({
            'ID': acc,
            'Length': len(seq),
            'GC%': f"{calculate_gc_content(seq):.2f}%",
            'Protein Count': estimate_protein_count(seq)
        })
    return pd.DataFrame(results_data, columns=columns)