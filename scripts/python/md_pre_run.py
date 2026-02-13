# Файл для проверки параметров перед отправкой на GPU
import os

def check_md_params(metal_name, softness):
    print(f"--- [LOCAL CHECK] Подготовка данных для GPU ---")
    # Рассчитываем параметры для силового поля MD
    # Энергия Ван-дер-Ваальса для MD на основе DFT данных
    vdw_radius = 1.75 if metal_name == "Pt" else 1.66
    
    print(f"Металл: {metal_name}")
    print(f"Квантовая мягкость (DFT): {softness}")
    print(f"Радиус Ван-дер-Ваальса: {vdw_radius} Å")
    
    # Путь для сохранения будущей траектории
    traj_path = f"Project_HeLa_Research/reports/{metal_name}_trajectory.dcd"
    print(f"Будущий файл траектории: {traj_path}")

# Запуск проверки для Платины
check_md_params("Pt", 1.232)