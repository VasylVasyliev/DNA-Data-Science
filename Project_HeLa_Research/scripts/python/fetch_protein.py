import os
import urllib.request

# Настройка путей
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../../"))
SAVE_DIR = os.path.join(PROJECT_ROOT, "data/protein_pdb")

# Список белков исследования 2026
# 1ao6: Альбумин (Транспорт)
# 1ivo: EGFR (Рецептор рака HeLa)
# 1qx3: Caspase-3 (Маркер гибели клетки)
PROTEINS = {
    "1ao6": "Human Serum Albumin",
    "1ivo": "EGFR Receptor",
    "1qx3": "Caspase-3"
}

def sync_proteins():
    os.makedirs(SAVE_DIR, exist_ok=True)
    
    for pdb_id, name in PROTEINS.items():
        file_path = os.path.join(SAVE_DIR, f"{pdb_id}.pdb")
        
        if os.path.exists(file_path):
            print(f"[OK] {pdb_id} ({name}) уже загружен.")
        else:
            print(f"[...] Загрузка {pdb_id} ({name})...")
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            try:
                urllib.request.urlretrieve(url, file_path)
                print(f"[Успех] Сохранено в {file_path}")
            except Exception as e:
                print(f"[Ошибка] Не удалось загрузить {pdb_id}: {e}")

if __name__ == "__main__":
    sync_proteins()