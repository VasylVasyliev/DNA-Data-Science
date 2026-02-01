import numpy as np
import matplotlib.pyplot as plt
import os

# Используем те же константы для честности
PARAMS = {
    "AU": {"eps": 5.29, "sig": 2.951}, 
    "PT": {"eps": 8.00, "sig": 2.754}, 
    "S":  {"eps": 0.25, "sig": 3.55}  
}

def lj_potential(r, m_type):
    eps = np.sqrt(PARAMS[m_type]["eps"] * PARAMS["S"]["eps"])
    sig = (PARAMS[m_type]["sig"] + PARAMS["S"]["sig"]) / 2
    return 4 * eps * ((sig/r)**12 - (sig/r)**6)

# Создаем диапазон расстояний от 2.5 до 8.0 Ангстрем
r_axis = np.linspace(2.8, 8.0, 100)
energy_au = [lj_potential(r, "AU") for r in r_axis]
energy_pt = [lj_potential(r, "PT") for r in r_axis]

plt.figure(figsize=(10, 6))
plt.plot(r_axis, energy_au, label='Gold (Au) Potential', color='#FFD700', lw=2)
plt.plot(r_axis, energy_pt, label='Platinum (Pt) Potential', color='#708090', lw=2)

# Линия нуля
plt.axhline(0, color='black', linestyle='--', alpha=0.3)

# Оформление
plt.title("Энергетический ландшафт взаимодействия Металл-EGFR (HeLa Project 2026)")
plt.xlabel("Расстояние между центрами атомов (Å)")
plt.ylabel("Энергия Леннард-Джонса (kcal/mol)")
plt.ylim(-2.0, 1.0) # Масштабируем, чтобы видеть "дно" ямы
plt.legend()
plt.grid(True, alpha=0.2)

# Сохранение
save_path = "Project_HeLa_Research/reports/docking_energy_plot.png"
os.makedirs(os.path.dirname(save_path), exist_ok=True)
plt.savefig(save_path)
print(f"📊 График сохранен в: {save_path}")