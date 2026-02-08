import os, requests, datetime, re
import matplotlib.pyplot as plt

# Конфигурация
REPOS = ["DNA-Data-Science", "HELA_Pt_Project"]
USERNAME = "VasylVasyliev"
TOKEN = "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"
OUTPUT_PATH = "results/traffic_comparison.png"

def update_data():
    today = datetime.date.today().strftime("%Y-%m-%d")
    print(f"🚀 Запуск синхронизации: {today}")
    
    # 1. Получаем данные и рисуем график
    counts = []
    for repo in REPOS:
        url = f"https://api.github.com/repos/{USERNAME}/{repo}/traffic/views"
        r = requests.get(url, headers={"Authorization": f"token {TOKEN}"})
        counts.append(r.json().get('count', 0) if r.status_code == 200 else 0)
    
    plt.figure(figsize=(10, 6))
    plt.bar(REPOS, counts, color=['#4A90E2', '#50E3C2'])
    os.makedirs("results", exist_ok=True)
    plt.savefig(OUTPUT_PATH)
    print("✅ График обновлен.")

    # 2. Обновляем README
    if os.path.exists("README.md"):
        with open("README.md", "r", encoding="utf-8") as f:
            content = f.read()
        # Ищем 'Last synced' и меняем дату после него
        new_content = re.sub(r"Last synced: \d{4}-\d{2}-\d{2}", f"Last synced: {today}", content)
        with open("README.md", "w", encoding="utf-8") as f:
            f.write(new_content)
        print(f"✅ Дата в README обновлена на {today}")

if __name__ == "__main__":
    update_data()
