import numpy as np
import matplotlib.pyplot as plt
import os

# EN: Use absolute path to avoid confusion between directories
# RU: Используем абсолютный путь, чтобы избежать путаницы с директориями
PROJECT_ROOT = "/Users/vasylvasyliev/Documents/Project_HeLa_Research"
SAVE_PATH = os.path.join(PROJECT_ROOT, "reports/docking_energy_plot.png")

# EN: Force delete old file if exists
# RU: Принудительное удаление старого файла, если он существует
if os.path.exists(SAVE_PATH):
    os.remove(SAVE_PATH)
    print(f"🗑 EN: Old plot removed. / RU: Старый график удален.")

# EN: Physical constants including Silver (Ag)
# RU: Физические константы, включая Серебро (Ag)
PARAMS = {
    "AU": {"eps": 5.29, "sig": 2.951}, 
    "PT": {"eps": 8.00, "sig": 2.754}, 
    "AG": {"eps": 4.40, "sig": 2.940}, 
    "S":  {"eps": 0.25, "sig": 3.55}  
}

def lj_potential(r, m_type):
    eps = np.sqrt(PARAMS[m_type]["eps"] * PARAMS["S"]["eps"])
    sig = (PARAMS[m_type]["sig"] + PARAMS["S"]["sig"]) / 2
    return 4 * eps * ((sig/r)**12 - (sig/r)**6)

r_axis = np.linspace(2.8, 8.0, 200)
plt.figure(figsize=(10, 6))
plt.style.use('dark_background')

# EN: Plotting Platinum, Gold, and the new Silver curve
# RU: Построение графиков для Платины, Золота и новой кривой Серебра
plt.plot(r_axis, [lj_potential(r, "PT") for r in r_axis], label='Platinum (Pt)', color='#708090')
plt.plot(r_axis, [lj_potential(r, "AG") for r in r_axis], label='Silver (Ag) Potential', color='#C0C0C0', lw=3, ls='--')
plt.plot(r_axis, [lj_potential(r, "AU") for r in r_axis], label='Gold (Au)', color='#FFD700')

plt.axhline(0, color='white', alpha=0.3)
plt.title("EGFR Energy Landscape: Silver (Ag) Included")
plt.ylim(-2.0, 1.0)
plt.legend()

# EN: Ensure directory exists and save
# RU: Убеждаемся, что директория существует, и сохраняем
os.makedirs(os.path.dirname(SAVE_PATH), exist_ok=True)
plt.savefig(SAVE_PATH)
print(f"✅ EN: New plot saved to: {SAVE_PATH} / RU: Новый график сохранен в: {SAVE_PATH}")