import os
import requests
import datetime
import re
import matplotlib.pyplot as plt
from dotenv import load_dotenv

# English: Load environment variables from .env file
# Russian: Загружаем переменные окружения из файла .env
load_dotenv()

# English: Get token from environment, avoid hardcoding for security
# Russian: Получаем токен из окружения, избегаем жесткой вставки для безопасности
TOKEN = os.getenv("GITHUB_TOKEN")
USERNAME = "VasylVasyliev"
REPOS = ["DNA-Data-Science", "HELA_Pt_Project"]
OUTPUT_PATH = "results/traffic_comparison.png"

# English: Debug line to check if the 'black box' works
# Russian: Отладочная строка, чтобы проверить, работает ли «черный ящик»
print(f"--- DEBUG: Token loaded: {'YES (Success)' if TOKEN else 'NO (Check .env file)'} ---")

def get_traffic_data(repo):
    url = f"https://api.github.com/repos/{USERNAME}/{repo}/traffic/views"
    headers = {
        "Authorization": f"token {TOKEN}",
        "Accept": "application/vnd.github.v3+json"
    }
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()['views']
    else:
        print(f"⚠️ Error fetching {repo}: {response.status_code}")
        return []

def update_readme():
    today = datetime.date.today().strftime("%Y-%m-%d")
    readme_path = "README.md"
    if os.path.exists(readme_path):
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()
        
        # English: Update the date in the README using regex
        # Russian: Обновляем дату в README с помощью регулярного выражения
        new_content = re.sub(r"Last Update: \d{4}-\d{2}-\d{2}", f"Last Update: {today}", content)
        
        with open(readme_path, "w", encoding="utf-8") as f:
            f.write(new_content)
        print(f"✅ README date updated to {today}")

def plot_traffic(data_dict):
    plt.figure(figsize=(10, 5))
    for repo, views in data_dict.items():
        dates = [v['timestamp'][:10] for v in views]
        counts = [v['count'] for v in views]
        plt.plot(dates, counts, marker='o', label=repo)
    
    plt.title("GitHub Traffic Comparison (2026)")
    plt.xlabel("Date")
    plt.ylabel("Views")
    plt.legend()
    plt.grid(True)
    os.makedirs("results", exist_ok=True)
    plt.savefig(OUTPUT_PATH)
    print(f"✅ Chart saved to {OUTPUT_PATH}")

if __name__ == "__main__":
    print(f"🚀 Starting synchronization: {datetime.date.today()}")
    
    if not TOKEN:
        print("❌ CRITICAL ERROR: GITHUB_TOKEN not found in .env!")
    else:
        all_data = {}
        for repo in REPOS:
            all_data[repo] = get_traffic_data(repo)
        
        plot_traffic(all_data)
        update_readme()