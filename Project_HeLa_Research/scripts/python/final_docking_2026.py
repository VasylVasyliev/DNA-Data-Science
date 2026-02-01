import numpy as np
from Bio.PDB import PDBParser
import os

# --- КОНФИГУРАЦИЯ ---
BASE_DIR = "/Users/vasylvasyliev/Documents/bio_ai-linux/Project_HeLa_Research"
PROTEIN_PATH = os.path.join(BASE_DIR, "data/protein_pdb/1ivo.pdb")
DATA_DIR = os.path.join(BASE_DIR, "data/nanoparticle")

# ФИЗИЧЕСКИЕ ПАРАМЕТРЫ (Lennard-Jones)
# eps: глубина потенциальной ямы (kcal/mol), sig: равновесное расстояние (A)
PARAMS = {
    "AU": {"eps": 5.29, "sig": 2.951}, 
    "PT": {"eps": 8.00, "sig": 2.754}, 
    "S":  {"eps": 0.25, "sig": 3.55}  
}

def calc_lj_energy(r, eps_m, sig_m, eps_s, sig_s):
    """Классическая формула потенциала Леннард-Джонса 12-6"""
    eps = np.sqrt(eps_m * eps_s) 
    sig = (sig_m + sig_s) / 2    
    
    # Чтобы избежать бесконечной энергии при очень малых r
    if r < 0.1: return 0
    
    term6 = (sig / r)**6
    term12 = term6**2
    return 4 * eps * (term12 - term6)

def run_physics_docking(protein_pdb, particle_pdb):
    parser = PDBParser(QUIET=True)
    
    try:
        prot_struct = parser.get_structure("P", protein_pdb)
        part_struct = parser.get_structure("M", particle_pdb)
    except Exception as e:
        return None

    # Собираем мишени (Сера)
    s_coords = [a.get_coord() for a in prot_struct.get_atoms() if "SG" in a.get_name()]
    if not s_coords: return None
    target_site = s_coords[0]

    # Определяем параметры металла
    element = "AU" if "Au" in particle_pdb else "PT"
    p = PARAMS[element]
    ps = PARAMS["S"]

    # Исходные координаты частицы
    m_coords_orig = np.array([a.get_coord() for a in part_struct.get_atoms()])
    m_center = np.mean(m_coords_orig, axis=0)

    best_energy = float('inf')

    # ПОИСК ИСТИНЫ: Сканируем расстояние от 2.0 до 6.0 Ангстрем
    # Это имитирует мягкое приближение частицы к белку
    for dist in np.linspace(2.0, 6.0, 41):
        current_total_energy = 0
        # Смещаем частицу на дистанцию 'dist' от первого атома серы
        m_coords = m_coords_orig - m_center + target_site + np.array([dist, 0, 0])
        
        for m in m_coords:
            for s in s_coords:
                r = np.linalg.norm(m - s)
                if r < 10.0: # Считаем только ближние взаимодействия
                    current_total_energy += calc_lj_energy(r, p["eps"], p["sig"], ps["eps"], ps["sig"])
        
        if current_total_energy < best_energy:
            best_energy = current_total_energy
            
    return best_energy

if __name__ == "__main__":
    print("\n" + "="*60)
    print("  SCIENTIFIC BINDING ENERGY ANALYSIS (kcal/mol)  ")
    print("="*60)

    particles = {
        "Gold (Au)": "Au_2018_model.pdb",
        "Platinum (Pt)": "Pt_2018_model.pdb"
    }

    results = {}

    for name, p_file in particles.items():
        path = os.path.join(DATA_DIR, p_file)
        if os.path.exists(path):
            energy = run_physics_docking(PROTEIN_PATH, path)
            results[name] = energy
            status = "ПРИТЯЖЕНИЕ" if energy < 0 else "ОТТАЛКИВАНИЕ"
            print(f">>> {name:15}: {energy:10.2f} kcal/mol [{status}]")

    print("\n" + "="*60)
    if len(results) == 2:
        winner = min(results, key=results.get)
        print(f"ВЕРДИКТ: {winner} обладает более сильным сродством к EGFR.")
        print("Это физически обоснованный показатель токсичности.")
    print("="*60 + "\n")