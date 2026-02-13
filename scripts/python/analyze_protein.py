import os

# Пути к данным
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROTEIN_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../../data/protein_pdb"))

# Список наших белков
protein_files = {
    "1ao6.pdb": "Albumin (Transport)",
    "1ivo.pdb": "EGFR (Cancer Receptor)",
    "1qx3.pdb": "Caspase-3 (Apoptosis)"
}

def analyze_all_proteins():
    print(f"{'Protein':<10} | {'Type':<20} | {'CYS (S-sites)':<12} | {'HIS (N-sites)':<12}")
    print("-" * 65)

    for file_name, p_type in protein_files.items():
        file_path = os.path.join(PROTEIN_DIR, file_name)
        
        if not os.path.exists(file_path):
            print(f"File {file_name} not found!")
            continue

        cys_s = 0  # Точки для Золота (Сера)
        his_n = 0  # Точки для Платины (Азот)

        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    
                    # Считаем серу в Цистеине (идеал для Au)
                    if res_name == "CYS" and "SG" in atom_name:
                        cys_s += 1
                    
                    # Считаем азот в Гистидине (важно для Pt в 2026 году)
                    if res_name == "HIS" and ("ND1" in atom_name or "NE2" in atom_name):
                        his_n += 1
        
        print(f"{file_name:<10} | {p_type:<20} | {cys_s:<12} | {his_n:<12}")

if __name__ == "__main__":
    analyze_all_proteins()