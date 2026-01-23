import os

def analyze_portfolio():
    data_dir = "data/"
    if not os.path.exists(data_dir):
        print("Папка data не найдена!")
        return

    stats = []
    # Собираем все файлы .txt
    files = [f for f in os.listdir(data_dir) if f.endswith(".txt")]
    
    for filename in sorted(files):
        path = os.path.join(data_dir, filename)
        with open(path, "r") as f:
            content = f.read().strip()
            total_len = len(content.replace('\n', '').replace(' ', ''))
            
            # Проверяем, является ли файл последовательностью ДНК/РНК
            bases = set("GCATU")
            is_sequence = any(c in bases for c in content.upper())

            if is_sequence and total_len > 0:
                g = content.upper().count('G')
                c = content.upper().count('C')
                gc_content = ((g + c) / total_len) * 100
                gc_str = f"{gc_content:.2f}%"
            else:
                gc_str = "N/A (Numeric/Other)"

            stats.append(f"| {filename:<12} | {total_len:>8} bp | {gc_str:>18} |")

    print("\n### 📊 Data Statistics (Full List: 14 Tasks)")
    print("| File Name    | Length      | GC-Content         |")
    print("|--------------|-------------|--------------------|")
    for row in stats:
        print(row)

if __name__ == "__main__":
    analyze_portfolio()
