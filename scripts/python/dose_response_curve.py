import numpy as np
import matplotlib.pyplot as plt

def get_survival_at_dose(dose_multiplier, base_tes=66.37):
    cells = 1000
    receptors = 100
    # Эффективность растет с дозой (упрощенная модель насыщения)
    eff = (base_tes * dose_multiplier) / 100
    if eff > 0.99: eff = 0.99 # Потолок эффективности
    
    surviving_rec = np.random.binomial(receptors, 1 - eff, cells)
    return np.sum(surviving_rec > 30) # Выжившие

# Проверяем разные дозы (от 0.5x до 3.0x от текущей)
doses = np.arange(0.5, 3.1, 0.25)
survival_stats = [get_survival_at_dose(d) for d in doses]

plt.figure(figsize=(10, 6))
plt.plot(doses, survival_stats, 'o-', color='#32CD32', lw=2)
plt.axhline(500, color='red', linestyle='--', label='Точка LD50 (50% гибели)')
plt.title("Зависимость выживаемости HeLa от концентрации Smart Pt")
plt.xlabel("Множитель дозы (от базового TES)")
plt.ylabel("Количество живых клеток (из 1000)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()