import numpy as np
import matplotlib.pyplot as plt
import os

# Константы жесткости (упрощенно)
K_PROTEIN = 2.0  # Насколько "трудно" согнуть аминокислоту
TIMESTEPS = 1000
DT = 0.05

def simulate_damage(energy_well):
    # Сила связи (чем глубже яма, тем сильнее тяга)
    k_bind = abs(energy_well) * 10.0
    
    pos_m = 1.0  # Стартовая позиция металла
    pos_p = 0.0  # Стартовая позиция белка (норма)
    
    m_traj = []
    p_traj = []
    
    for _ in range(TIMESTEPS):
        # Расстояние между металлом и белком
        dist = pos_m - pos_p
        
        # Сила связи (тянет их друг к другу)
        f_bind = -k_bind * (dist - 0.5) # 0.5 - идеальный зазор
        
        # Сила упругости белка (сопротивляется деформации)
        f_prot = -K_PROTEIN * pos_p
        
        # Обновляем позиции (упрощенная динамика без инерции для ясности)
        pos_m += f_bind * DT
        pos_p += (f_prot - f_bind) * DT # Белок тянется за металлом
        
        m_traj.append(pos_m)
        p_traj.append(pos_p)
        
    return m_traj, p_traj

# Данные из твоего "мокрого" теста
m_au, p_au = simulate_damage(-0.493)
m_pt, p_pt = simulate_damage(-0.838)

plt.figure(figsize=(10, 6))
plt.plot(p_au, label='Деформация от Золота (Au)', color='#FFD700', lw=2)
plt.plot(p_pt, label='Деформация от Платины (Pt)', color='#708090', lw=3)

plt.axhline(0, color='black', linestyle='--')
plt.title("Деформация структуры рецептора EGFR (Induced Fit)")
plt.xlabel("Время взаимодействия")
plt.ylabel("Смещение аминокислот от нормы (Å)")
plt.legend()
plt.grid(True, alpha=0.2)

save_path = "Project_HeLa_Research/reports/protein_deformation.png"
os.makedirs(os.path.dirname(save_path), exist_ok=True)
plt.savefig(save_path)
print(f"✅ Анализ деформации сохранен: {save_path}")
plt.show()