import numpy as np
from Bio.PDB import PDBParser
import os

# Параметры среды (Water at 310K - Body Temp)
DIELECTRIC_WATER = 78.4 
SURFACE_TENSION = 0.025 # kcal/mol/A^2 (энергия на вытеснение воды)

PARAMS = {
    "AU": {"eps": 5.29, "sig": 2.951, "radius": 1.44}, 
    "PT": {"eps": 8.00, "sig": 2.754, "radius": 1.35}, 
    "S":  {"eps": 0.25, "sig": 3.55,  "radius": 1.80}  
}

def calc_advanced_energy(r, m_type):
    p = PARAMS[m_type]
    ps = PARAMS["S"]
    
    # 1. Потенциал Леннард-Джонса (Ван-дер-Ваальс)
    eps_mix = np.sqrt(p["eps"] * ps["eps"])
    sig_mix = (p["sig"] + ps["sig"]) / 2
    vdw = 4 * eps_mix * ((sig_mix/r)**12 - (sig_mix/r)**6)
    
    # 2. Поправка на воду (Solvation Penalty)
    # Чем ближе частица, тем больше воды она вытесняет (упрощенная модель SASA)
    # Мы считаем объем перекрытия сфер
    if r < (p["radius"] + ps["radius"] + 2.8): # 2.8A - диаметр молекулы воды
        solvation_cost = SURFACE_TENSION * (4 * np.pi * p["radius"]**2)
    else:
        solvation_cost = 0
        
    return vdw + solvation_cost

# --- Основной цикл расчета ---
def run_solvent_docking():
    # Мы имитируем сближение в воде
    r_range = np.linspace(2.5, 7.0, 50)
    
    print(f"{'Dist (A)':<10} | {'Au (Water) ':<15} | {'Pt (Water)':<15}")
    print("-" * 45)
    
    results = {"AU": [], "PT": []}
    
    for r in r_range:
        e_au = calc_advanced_energy(r, "AU")
        e_pt = calc_advanced_energy(r, "PT")
        results["AU"].append(e_au)
        results["PT"].append(e_pt)
        if r % 1.0 < 0.1: # Печатаем через каждый 1 Ангстрем
             print(f"{r:<10.1f} | {e_au:<15.4f} | {e_pt:<15.4f}")

    # Поиск минимума в воде
    min_au = min(results["AU"])
    min_pt = min(results["PT"])
    
    print("-" * 45)
    print(f"ФИНАЛЬНЫЙ БАЛАНС В ВОДЕ (kcal/mol):")
    print(f"Gold: {min_au:.3f} | Platinum: {min_pt:.3f}")
    
    if min_pt < min_au:
        diff = abs(min_pt - min_au)
        print(f"\nВЕРДИКТ: Платина сохраняет преимущество в воде на {diff:.3f} kcal/mol.")
    else:
        print("\nВЕРДИКТ: В воде золото сравнялось с платиной!")

if __name__ == "__main__":
    run_solvent_docking()