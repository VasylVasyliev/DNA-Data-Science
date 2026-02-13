import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio
import os

# Имитация данных (так как у нас нет доступа к твоим логам прямо сейчас)
# В реальном проекте здесь будет загрузка из твоего .xvg или .log файла
frames = np.linspace(0, 100, 50)
energy_au = -0.5 + 0.1 * np.random.randn(50)  # Данные для Gold
energy_pt = -0.838 + 0.05 * np.random.randn(50) # Данные для Platinum

filenames = []

print("Создаем кадры анимации...")
for i in range(1, len(frames)):
    plt.figure(figsize=(10, 6))
    plt.plot(frames[:i], energy_au[:i], color='gold', label='Gold (Au)', linewidth=2)
    plt.plot(frames[:i], energy_pt[:i], color='silver', label='Platinum (Pt)', linewidth=2)
    
    plt.ylim(-1.2, 0)
    plt.xlim(0, 100)
    plt.title('EGFR Binding Stability Dynamics (37°C)', fontsize=14)
    plt.xlabel('Simulation Time (ps)', fontsize=12)
    plt.ylabel('Interaction Energy', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Сохраняем каждый кадр
    filename = f'frame_{i}.png'
    plt.savefig(filename)
    filenames.append(filename)
    plt.close()

# Собираем в GIF
print("Склеиваем в GIF...")
with imageio.get_writer('md_stability_animated.gif', mode='I', duration=0.1) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
        os.remove(filename) # Удаляем временные картинки

print("Готово! Файл md_stability_animated.gif создан.")