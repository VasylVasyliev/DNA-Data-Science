import os

# Тот же абсолютный путь
DATA_DIR = "/Users/vasylvasyliev/Documents/bio_ai-linux/Project_HeLa_Research/data/nanoparticle"

def convert():
    files = {"Au_2018_model": "AU", "Pt_2018_model": "PT"}
    for name, element in files.items():
        xyz_path = os.path.join(DATA_DIR, f"{name}.xyz")
        pdb_path = os.path.join(DATA_DIR, f"{name}.pdb")

        if os.path.exists(xyz_path):
            with open(xyz_path, 'r') as f:
                lines = f.readlines()[2:]
            with open(pdb_path, 'w') as f:
                for i, line in enumerate(lines):
                    p = line.split()
                    f.write(f"ATOM  {i+1:5d}  {element:<3} NPS A   1    {float(p[1]):8.3f}{float(p[2]):8.3f}{float(p[3]):8.3f}  1.00  0.00          {element:>2}\n")
                f.write("TER\nEND\n")
            print(f"🚀 Сконвертировано в PDB: {pdb_path}")
        else:
            print(f"❌ Всё еще не вижу: {xyz_path}")

if __name__ == "__main__":
    convert()