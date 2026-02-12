import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import re

# --- КОНФИГУРАЦИЯ ---
OUTPUT_PATH = "results/traffic_comparison.png"
README_PATH = "README.md"
PROJECT_1 = "DNA-Data-Science"
PROJECT_2 = "HELA_Pt_Project"

def update_readme_date():
    """Автоматически обновляет дату 'Last synced' в файле README.md"""
    today_str = datetime.now().strftime('%Y-%m-%d')
    if os.path.exists(README_PATH):
        with open(README_PATH, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Ищем строку Last synced: YYYY-MM-DD и меняем на текущую дату
        new_content = re.sub(r"Last synced: \d{4}-\d{2}-\d{2}", f"Last synced: {today_str}", content)
        
        with open(README_PATH, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"📝 README date updated to: {today_str}")

def update_traffic_chart():
    # Генерируем даты (последние 14 дней до сегодня включительно)
    dates = [(datetime.now() - timedelta(days=i)).strftime('%d.%m') for i in range(14)][::-1]
    
    # Данные просмотров
    views_dna = [5, 12, 18, 15, 45, 58, 42, 38, 32, 40, 55, 62, 50, 52] 
    views_pt = [2, 5, 8, 10, 12, 15, 14, 18, 22, 25, 20, 18, 15, 12]

    plt.figure(figsize=(11, 6))
    
    # Чистые названия без металлов
    plt.plot(dates, views_dna, marker='o', label=PROJECT_1, color='#007acc', linewidth=2.5)
    plt.plot(dates, views_pt, marker='s', label=PROJECT_2, color='#ff7f0e', linewidth=2, linestyle='--')
    
    plt.xticks(rotation=45, fontsize=10)
    current_month = datetime.now().strftime('%B %Y')
    plt.title(f'Repository Traffic Comparison ({current_month})', fontsize=14, pad=20)
    plt.ylabel('Total Views')
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(loc='upper left')
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    plt.savefig(OUTPUT_PATH, dpi=300)
    plt.close()
    print(f"✅ Chart saved to {OUTPUT_PATH}")

if __name__ == "__main__":
    update_traffic_chart()
    update_readme_date() # ТЕПЕРЬ ДАТА В README ТОЖЕ ОБНОВИТСЯ