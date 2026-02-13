import numpy as np
import matplotlib.pyplot as plt

def simulate_cell_death(tes_score):
    # Начальная популяция HeLa
    cells = 1000
    # Рецепторов на одну клетку
    receptors_per_cell = 100
    # Эффективность блокировки (базируется на нашем TES)
    block_efficiency = tes_score / 100
    
    # Симуляция: сколько рецепторов выжило у каждой клетки
    # Используем биномиальное распределение
    surviving_receptors = np.random.binomial(receptors_per_cell, 1 - block_efficiency, cells)
    
    # Клетка выживает, если у нее осталось более 30 работающих рецепторов
    # (Критический порог выживания раковой клетки)
    alive_cells = np.sum(surviving_receptors > 30)
    dead_cells = cells - alive_cells
    
    return alive_cells, dead_cells

# Сравниваем Золото и Умную Платину
au_alive, au_dead = simulate_cell_death(19.22)
pt_smart_alive, pt_smart_dead = simulate_cell_death(66.37)

print(f"--- Результаты воздействия на популяцию (1000 клеток HeLa) ---")
print(f"Gold (Au): Выжило {au_alive} | Погибло {au_dead}")
print(f"Smart Pt: Выжило {pt_smart_alive} | Погибло {pt_smart_dead}")