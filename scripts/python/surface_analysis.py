import os
from Bio.PDB import PDBParser

# Настройка путей
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Путь к белку EGFR (1ivo)
PDB_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, "../../data/protein_pdb/1ivo.pdb"))
# Путь для сохранения отчета
REPORT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../results"))

def run_surface_scan():
    if not os.path.exists(PDB_PATH):
        print(f"Ошибка: Файл {PDB_PATH} не найден!")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("EGFR", PDB_PATH)
    model = structure[0]

    exposed_cys = 0
    buried_cys = 0

    print(f"\n--- Анализ доступности поверхности EGFR (1ivo) ---")
    
    # В исследовании 2026 мы используем упрощенный геометрический фильтр 
    # для предварительной оценки перед Docking
    for residue in model.get_residues():
        if residue.get_resname() == "CYS":
            if "SG" in residue:
                # Считаем расстояние от центра масс или используем Z-фильтр
                coord = residue["SG"].get_coord()
                # Простая логика: атомы на периферии (далеко от центра 0,0,0) считаем открытыми
                dist_from_center = (coord[0]**2 + coord[1]**2 + coord[2]**2)**0.5
                
                if dist_from_center > 25: # Порог доступности для наночастицы 2нм
                    exposed_cys += 1
                else:
                    buried_cys += 1

    print(f"Всего сайтов CYS (сера): {exposed_cys + buried_cys}")
    print(f"Открыты для удара (Exposed): {exposed_cys} ✅")
    print(f"Спрятаны внутри (Buried): {buried_cys} ❌")

    # Сохраняем в папку results
    os.makedirs(REPORT_DIR, exist_ok=True)
    with open(os.path.join(REPORT_DIR, "surface_report.txt"), "w") as f:
        f.write(f"EGFR Surface Analysis:\nExposed: {exposed_cys}\nBuried: {buried_cys}\n")

if __name__ == "__main__":
    run_surface_scan()