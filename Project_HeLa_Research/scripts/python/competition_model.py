import numpy as np
import matplotlib.pyplot as plt

# Параметры сродства (Affinity)
# Чем выше число, тем охотнее металл связывается
EGFR_AFFINITY = {"AU": 0.493, "PT": 0.838}
GSH_AFFINITY =  {"AU": 0.300, "PT": 0.950} # Платина ОЧЕНЬ любит ловушки (GSH)

def simulate_cellular_race(metal):
    particles = 1000 # Запускаем 1000 "доз" металла
    trapped = 0
    reached_target = 0
    
    # Вероятность связывания пропорциональна сродству
    prob_trap = GSH_AFFINITY[metal] / (GSH_AFFINITY[metal] + 0.5)
    prob_target = EGFR_AFFINITY[metal] / (EGFR_AFFINITY[metal] + 0.5)

    for _ in range(particles):
        # 1. Попытка пройти через ловушки
        if np.random.random() < prob_trap:
            trapped += 1
        # 2. Если прошел, попытка связаться с целью
        elif np.random.random() < prob_target:
            reached_target += 1
            
    return trapped, reached_target

# Запуск симуляции
au_trapped, au_target = simulate_cellular_race("AU")
pt_trapped, pt_target = simulate_cellular_race("PT")

# Визуализация
labels = ['Gold (Au)', 'Platinum (Pt)']
trapped_data = [au_trapped, pt_trapped]
target_data = [au_target, pt_target]

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(x - width/2, trapped_data, width, label='Поймано ловушками (GSH)', color='#d9534f')
ax.bar(x + width/2, target_data, width, label='Достигло EGFR (Удар по раку)', color='#5cb85c')

ax.set_ylabel('Количество частиц (из 1000)')
ax.set_title('Эффективность доставки металла к мишени в клетке HeLa')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

plt.grid(axis='y', alpha=0.3)
plt.savefig("Project_HeLa_Research/reports/competition_results.png")
print("✅ Анализ конкуренции завершен и сохранен.")
plt.show()