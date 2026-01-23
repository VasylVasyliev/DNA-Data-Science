import matplotlib.pyplot as plt

def plot_mismatches(file_path):
    with open(file_path, "r") as f:
        s = f.readline().strip()
        t = f.readline().strip()

    # Создаем список: 1 если буквы разные, 0 если одинаковые
    mismatches = [1 if a != b else 0 for a, b in zip(s, t)]
    
    plt.figure(figsize=(15, 3))
    # Рисуем "штрих-код" мутаций
    plt.vlines([i for i, val in enumerate(mismatches) if val == 1], 0, 1, colors='red', alpha=0.5)
    
    plt.title(f"Mutation Map (Hamming Distance: {sum(mismatches)})")
    plt.xlabel("Position in Sequence")
    plt.yticks([]) # Убираем лишние оси
    
    save_path = "data/plots/hamm_map.png"
    plt.savefig(save_path)
    print(f"Mutation map saved to {save_path}")

if __name__ == "__main__":
    plot_mismatches("data/hamm.txt")
    