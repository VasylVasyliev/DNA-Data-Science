import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# --- КОНФИГУРАЦИЯ / CONFIGURATION ---
OUTPUT_PATH = "results/traffic_comparison.png"
PROJECT_1 = "DNA-Data-Science"
PROJECT_2 = "HELA_Pt_Project"

def update_traffic_chart():
    # Генерируем даты за последние 14 дней (до сегодняшнего числа)
    # RU: Автоматическое обновление дат для оси X
    dates = [(datetime.now() - timedelta(days=i)).strftime('%d.%m') for i in range(14)][::-1]
    
    # Данные просмотров (замените на актуальные из GitHub Insights при необходимости)
    views_dna = [5, 12, 18, 15, 45, 58, 42, 38, 32, 40, 55, 62, 50, 52] 
    views_pt = [2, 5, 8, 10, 12, 15, 14, 18, 22, 25, 20, 18, 15, 12]

    plt.figure(figsize=(11, 6))
    
    # Исправленные названия линий (только названия репозиториев)
    plt.plot(dates, views_dna, marker='o', label=PROJECT_1, 
             color='#007acc', linewidth=2.5, markersize=7)
    plt.plot(dates, views_pt, marker='s', label=PROJECT_2, 
             color='#ff7f0e', linewidth=2, linestyle='--', alpha=0.8)
    
    # Настройка осей и сетки
    plt.xticks(rotation=45, fontsize=10)
    plt.yticks(fontsize=10)
    
    # Динамический заголовок с текущим месяцем и годом
    current_date_str = datetime.now().strftime('%B %Y')
    plt.title(f'Repository Traffic Comparison ({current_date_str})', fontsize=14, pad=20)
    plt.xlabel('Date (Day.Month)', fontsize=12)
    plt.ylabel('Total Views', fontsize=12)
    
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(fontsize=10, loc='upper left')
    plt.tight_layout()
    
    if not os.path.exists('results'):
        os.makedirs('results')
    
    plt.savefig(OUTPUT_PATH, dpi=300)
    plt.close()
    print(f"✅ Success: Chart updated for {current_date_str}")

if __name__ == "__main__":
    update_traffic_chart()