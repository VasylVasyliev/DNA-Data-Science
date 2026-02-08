import os
import requests
import datetime
import re
import matplotlib.pyplot as plt

# --- Конфигурация ---
# English: Repository list and credentials
# Russian: Список репозиториев и данные доступа
REPOS = ["DNA-Data-Science", "HELA_Pt_Project"]
USERNAME = "VasylVasyliev"
TOKEN = "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"
OUTPUT_PATH = "results/traffic_comparison.png"

def update_dashboard():
    # English: Get today's date for synchronization
    # Russian: Получаем сегодняшнюю дату для синхронизации
    today = datetime.date.today().strftime("%Y-%m-%d")
    print(f"🚀 Starting synchronization: {today}")
    
    counts = []
    headers = {"Authorization": f"token {TOKEN}"}

    # 1. English: Fetch data from GitHub API
    # 1. Russian: Получаем данные из GitHub API
    for repo in REPOS:
        url = f"https://api.github.com/repos/{USERNAME}/{repo}/traffic/views"
        try:
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                data = response.json()
                counts.append(data.get('count', 0))
            else:
                print(f"⚠️ Error fetching {repo}: {response.status_code}")
                counts.append(0)
        except Exception as e:
            print(f"❌ Connection error for {repo}: {e}")
            counts.append(0)

    # 2. English: Create chart with data labels
    # 2. Russian: Создаем график с числовыми метками
    plt.figure(figsize=(10, 6))
    bars = plt.bar(REPOS, counts, color=['#4A90E2', '#50E3C2'])
    
    # English: Add exact numbers above each bar
    # Russian: Добавляем точные числа над каждым столбиком
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                 f'{int(height)}', ha='center', va='bottom', 
                 fontsize=12, fontweight='bold', color='#333333')
    
    plt.title(f'Research Activity (Last 14 Days) - {today}', fontsize=14)
    plt.ylabel('Total Views', fontsize=12)
    plt.ylim(0, max(counts) * 1.2 if counts and max(counts) > 0 else 10)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # English: Save the visualization
    # Russian: Сохраняем визуализацию
    os.makedirs("results", exist_ok=True)
    plt.savefig(OUTPUT_PATH)
    plt.close()
    print(f"✅ Chart saved to {OUTPUT_PATH}")

    # 3. English: Update the date in README.md
    # 3. Russian: Обновляем дату в файле README.md
    if os.path.exists("README.md"):
        with open("README.md", "r", encoding="utf-8") as f:
            content = f.read()
        
        # English: Regex to find and replace the date
        # Russian: Регулярное выражение для поиска и замены даты
        new_content = re.sub(r"Last synced: \d{4}-\d{2}-\d{2}", f"Last synced: {today}", content)
        
        with open("README.md", "w", encoding="utf-8") as f:
            f.write(new_content)
        print(f"✅ README date updated to {today}")

if __name__ == "__main__":
    update_dashboard()