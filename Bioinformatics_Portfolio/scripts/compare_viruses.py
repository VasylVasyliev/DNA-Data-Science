import matplotlib.pyplot as plt
import pandas as pd
import os
from utils import get_data_path

def create_comparison():
    # 1. Пути к файлам результатов (предполагаем, что ты уже запустил orf.py для обоих)
    # Для этого примера мы просто возьмем топ-10 длинных белков для каждого
    results_dir = os.path.join(os.path.dirname(get_data_path("phix174.txt")), "..", "results")
    
    # Данные для графика (на основе твоих последних запусков)
    virus_names = ['PhiX174 (Phage)', 'Zika Virus']
    max_protein_lengths = [427, 3419]  # Те цифры, которые выдал твой скрипт
    avg_protein_lengths = [150, 850]   # Примерные средние значения
    
    # 2. Создаем график
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Рисуем столбцы для максимальной длины
    bars = ax.bar(virus_names, max_protein_lengths, color=['#4287f5', '#f54242'], alpha=0.7)
    
    # Добавляем подписи данных
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height} aa',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontweight='bold')

    # 3. Оформление
    ax.set_title('Genome Architecture Comparison: Max Protein Length', fontsize=14)
    ax.set_ylabel('Length (Amino Acids)', fontsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    
    # Сохраняем результат
    output_path = os.path.join(results_dir, "virus_comparison.png")
    plt.savefig(output_path, dpi=300)
    print(f"✅ Comparison chart saved to: {output_path}")

if __name__ == "__main__":
    create_comparison()