import numpy as np
import matplotlib.pyplot as plt
import os

# Путь для сохранения графиков
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, "../../results/plots"))
os.makedirs(SAVE_PATH, exist_ok=True)

def generate_xrd_pattern(name, peaks, color):
    x = np.linspace(20, 90, 1000)
    y = np.zeros_like(x)
    
    # Моделируем пики дифракции (функция Лоренца)
    for peak in peaks:
        y += 100 * (0.5**2 / ((x - peak)**2 + 0.5**2))
    
    plt.figure(figsize=(10, 5))
    plt.plot(x, y, color=color, lw=2, label=f'{name} Nanoparticles (Model 2026)')
    plt.fill_between(x, y, color=color, alpha=0.2)
    plt.title(f"Simulated XRD Pattern: {name}", fontsize=14)
    plt.xlabel("2-Theta (degrees)", fontsize=12)
    plt.ylabel("Intensity (a.u.)", fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    file_name = os.path.join(SAVE_PATH, f"{name.lower()}_xrd.png")
    plt.savefig(file_name)
    print(f"[Успех] График XRD для {name} сохранен в: {file_name}")
    plt.close()

if __name__ == "__main__":
    # Пики для Au (золото) и Pt (платина) - стандартные плоскости (111), (200), (220)
    generate_xrd_pattern("Gold", [38.2, 44.4, 64.6, 77.5], "gold")
    generate_xrd_pattern("Platinum", [39.8, 46.2, 67.5, 81.3], "slategray")