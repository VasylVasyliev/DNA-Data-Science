import numpy as np
import matplotlib.pyplot as plt
import os

# Физические константы
TEMP = 310  # Температура тела (K)
KB = 0.001987 # Константа Больцмана (kcal/mol*K)
TIMESTEPS = 1000
DT = 0.1

def run_stochastic_sim(energy_well):
    # Коэффициент "вязкости" (упрощенно)
    friction = 0.1
    pos = 0.0
    vel = 0.0
    trajectory = []
    
    # Жесткость связи (зависит от глубины ямы)
    k_binding = abs(energy_well) * 5.0 
    
    for _ in range(TIMESTEPS):
        # Случайный удар от молекул воды (Тепловой шум)
        thermal_force = np.random.normal(0, np.sqrt(2 * friction * KB * TEMP / DT))
        
        # Сила притяжения белка (возвращающая сила)
        binding_force = -k_binding * pos
        
        # Уравнение движения (Langevin-like)
        accel = (thermal_force + binding_force - friction * vel)
        vel += accel * DT
        pos += vel * DT
        trajectory.append(pos)
        
    return trajectory

# Твои данные из "мокрого" теста
traj_au = run_stochastic_sim(-0.493)
traj_pt = run_stochastic_sim(-0.838)

plt.figure(figsize=(12, 6))
plt.plot(traj_au, label=f'Gold (Au) - Энергия: -0.493', color='#FFD700', alpha=0.6)
plt.plot(traj_pt, label=f'Platinum (Pt) - Энергия: -0.838', color='#708090', alpha=0.9)

plt.axhline(0, color='black', lw=1, ls='--')
plt.title("Динамическая устойчивость на поверхности EGFR (Симуляция при 37°C)")
plt.xlabel("Время (условные фемтосекунды)")
plt.ylabel("Смещение от центра сайта (Å)")
plt.legend()
plt.grid(True, alpha=0.2)

# Сохранение результата
report_path = "Project_HeLa_Research/reports/md_stability.png"
plt.savefig(report_path)
print(f"✅ Динамический график сохранен: {report_path}")
plt.show()