import matplotlib.pyplot as plt
import os
# Импортируем наши инструменты из utils.py
from utils import read_fasta, get_data_path

def generate_dna_distribution(filename="phix174.txt"):
    """Counts nucleotides and saves a bar chart to the results folder."""
    try:
        # 1. Читаем последовательность
        record = read_fasta(filename)
        seq = record.seq
        
        # 2. Считаем количество нуклеотидов
        counts = {
            'Adenine (A)': seq.count('A'),
            'Cytosine (C)': seq.count('C'),
            'Guanine (G)': seq.count('G'),
            'Thymine (T)': seq.count('T')
        }
        
        # 3. Создаем столбчатую диаграмму
        # Используем мягкие цвета (палитра 'pastel')
        colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']
        plt.figure(figsize=(10, 6))
        bars = plt.bar(counts.keys(), counts.values(), color=colors, edgecolor='black', alpha=0.8)
        
        # Добавляем числа над каждым столбцом
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval + 10, int(yval), 
                     ha='center', va='bottom', fontsize=11, fontweight='bold')

        plt.title(f"Nucleotide Composition: {record.description[:40]}...", fontsize=14, fontweight='bold')
        plt.ylabel("Base Count", fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.3)
        
        # 4. Определяем путь сохранения в папку results
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        output_path = os.path.join(base_dir, "results", "dna_distribution.png")
        
        # Сохраняем файл
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✅ Success! Chart saved as: {output_path}")

    except Exception as e:
        print(f"❌ Error generating chart: {e}")

if __name__ == "__main__":
    generate_dna_distribution()