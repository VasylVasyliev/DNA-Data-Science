import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# --- КОНФИГУРАЦИЯ / CONFIGURATION ---
OUTPUT_PATH = "results/traffic_comparison.png"
PROJECT_1 = "DNA-Data-Science"
PROJECT_2 = "HELA_Pt_Project"

# --- НАУЧНЫЙ РАСЧЕТ / SCIENTIFIC CALCULATION ---
# RU: Программный расчет энергии для Ag100
# EN: Programmatic calculation of Ag100 energy
def get_silver_metrics():
    n_atoms = 100
    energy_per_atom = -0.00712
    mean_energy = round(n_atoms * energy_per_atom, 3)
    std_dev = 0.085
    return mean_energy, std_dev

# --- ГЕНЕРАЦИЯ ГРАФИКА / PLOT GENERATION ---
def update_traffic_chart():
    # Генерируем даты за последние 14 дней
    dates = [(datetime.now() - timedelta(days=i)).strftime('%d.%m') for i in range(14)][::-1]
    
    # Данные просмотров для двух проектов (на основе ваших скриншотов)
    views_dna = [5, 12, 18, 15, 45, 58, 42, 38, 32, 40, 55, 62, 50, 52] 
    views_pt = [2, 5, 8, 10, 12, 15, 14, 18, 22, 25, 20, 18, 15, 12]

    plt.figure(figsize=(11, 6)) # Немного увеличили ширину для дат
    
    # Рисуем две линии: Синюю (Серебро) и Оранжевую (Платина)
    plt.plot(dates, views_dna, marker='o', label=f'{PROJECT_1} (Ag)', 
             color='#007acc', linewidth=2.5, markersize=7)
    plt.plot(dates, views_pt, marker='s', label=f'{PROJECT_2} (Pt)', 
             color='#ff7f0e', linewidth=2, linestyle='--', alpha=0.8)
    
    # ИСПРАВЛЕНИЕ: Поворот дат на 45 градусов и настройка шрифта
    plt.xticks(rotation=45, fontsize=10)
    plt.yticks(fontsize=10)
    
    plt.title('Multi-Project Traffic Comparison (February 2026)', fontsize=14, pad=20)
    plt.xlabel('Date (Day.Month)', fontsize=12)
    plt.ylabel('Total Views', fontsize=12)
    
    plt.grid(True, linestyle=':', alpha=0.7)
    
    # ДОБАВЛЕНИЕ ЛЕГЕНДЫ
    plt.legend(fontsize=10, loc='upper left')
    
    # ИСПРАВЛЕНИЕ: Автоматические отступы, чтобы всё влезло в кадр
    plt.tight_layout()
    
    # Проверка и создание папки results
    if not os.path.exists('results'):
        os.makedirs('results')
    
    plt.savefig(OUTPUT_PATH, dpi=300)
    plt.close()
    print(f"✅ Success: Comparison chart saved to {OUTPUT_PATH}")

if __name__ == "__main__":
    print(f"--- Running Analytics Sync ---")
    
    # Выводим расчеты для проверки (будут видны в терминале)
    m_energy, s_dev = get_silver_metrics()
    print(f"🧬 Scientific Data (Ag100): Mean Energy = {m_energy}, Std Dev = {s_dev}")
    
    # Запускаем обновление графики
    update_traffic_chart()
    print(f"--- Sync Complete for {PROJECT_1} & {PROJECT_2} ---")