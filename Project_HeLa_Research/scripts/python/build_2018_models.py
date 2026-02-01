import os
import numpy as np

BASE_DIR = "/Users/vasylvasyliev/Documents/bio_ai-linux/Project_HeLa_Research"
TARGET_DIR = os.path.join(BASE_DIR, "data/nanoparticle")

def generate_dense_cluster(filename, element, num_atoms):
    path = os.path.join(TARGET_DIR, f"{filename}.xyz")
    with open(path, "w") as f:
        f.write(f"{num_atoms}\n")
        f.write(f"{element} dense cluster\n")
        for i in range(num_atoms):
            # Создаем плотный шар: атомы близко друг к другу
            phi = np.random.uniform(0, 2*np.pi)
            costheta = np.random.uniform(-1, 1)
            u = np.random.uniform(0, 1)
            
            theta = np.arccos(costheta)
            # Радиус зависит от количества атомов (кубический корень)
            r = (num_atoms**0.33) * 2.8 * (u**0.33)
            
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            f.write(f"{element} {x:.3f} {y:.3f} {z:.3f}\n")
    print(f"✅ Плотная модель {element} ({num_atoms} атомов) создана.")

if __name__ == "__main__":
    generate_dense_cluster("Au_2018_model", "Au", 13)
    generate_dense_cluster("Pt_2018_model", "Pt", 150) # Сделаем Платину реально массивной