import numpy as np
from Bio.PDB import PDBParser
import os

# Пути (используем проверенные абсолютные пути)
BASE_DIR = "/Users/vasylvasyliev/Documents/bio_ai-linux/Project_HeLa_Research"
PROTEIN_PATH = os.path.join(BASE_DIR, "data/protein_pdb/1ivo.pdb")
PARTICLE_PATH = os.path.join(BASE_DIR, "data/nanoparticle/Au_2018_model.pdb")

def run_docking_score(protein_pdb, particle_pdb):
    parser = PDBParser(QUIET=True)
    prot_struct = parser.get_structure("PROT", protein_pdb)
    part_struct = parser.get_structure("PART", particle_pdb)

    # Собираем координаты всех атомов серы (CYS) в белке
    target_atoms = []
    for model in prot_struct:
        for residue in model.get_residues():
            if residue.get_resname() == "CYS" and "SG" in residue:
                target_atoms.append(residue["SG"].get_coord())
    
    # Координаты атомов наночастицы
    particle_coords = [atom.get_coord() for atom in part_struct.get_atoms()]
    
    # Простой алгоритм: ищем среднее расстояние
    # Чем оно меньше, тем лучше "стыковка"
    min_avg_dist = float('inf')
    
    # Для теста просто посчитаем близость к самому доступному сайту
    for target in target_atoms:
        distances = [np.linalg.norm(target - p) for p in particle_coords]
        avg_dist = np.mean(distances)
        if avg_dist < min_avg_dist:
            min_avg_dist = avg_dist

    return min_avg_dist

if __name__ == "__main__":
    print("--- Запуск Docking анализа 2026 ---")
    
    # Считаем для Золота
    score_au = run_docking_score(PROTEIN_PATH, os.path.join(BASE_DIR, "data/nanoparticle/Au_2018_model.pdb"))
    # Считаем для Платины
    score_pt = run_docking_score(PROTEIN_PATH, os.path.join(BASE_DIR, "data/nanoparticle/Pt_2018_model.pdb"))
    
    print(f"Docking Score (Au-EGFR): {score_au:.2f} Å")
    print(f"Docking Score (Pt-EGFR): {score_pt:.2f} Å")
    
    if score_pt < score_au:
        print("\nРезультат: Платина имеет более высокий потенциал токсичности (лучшая стыковка).")
    else:
        print("\nРезультат: Золото связывается сильнее.")