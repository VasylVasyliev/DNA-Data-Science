import os
from Bio import Entrez, SeqIO
from utils import get_data_path

# Всегда указывай свой Email, иначе NCBI может заблокировать доступ
Entrez.email = "vasj5722814@gmail.com" 

def download_genome(accession_id):
    """Downloads a genome from NCBI by Accession ID."""
    print(f"🚀 Searching for {accession_id} in NCBI databases...")
    
    try:
        # 1. Запрос к базе данных Nucleotide (nuccore)
        with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")
            
        # 2. Формируем путь для сохранения
        filename = f"{accession_id}.fasta"
        # Сохраняем в папку data, используя наш utils.py
        output_path = os.path.join(os.path.dirname(get_data_path("phix174.txt")), filename)
        
        # 3. Записываем файл
        SeqIO.write(record, output_path, "fasta")
        print(f"✅ Success! Genome saved to: {output_path}")
        print(f"🧬 Organism: {record.description}")
        print(f"📏 Length: {len(record.seq)} bp")
        
    except Exception as e:
        print(f"❌ Error downloading from NCBI: {e}")

if __name__ == "__main__":
    # Давай скачаем геном вируса Зика (Zika virus) для примера
    # Его Accession ID: NC_012532
    target_id = "NC_012532"
    download_genome(target_id)