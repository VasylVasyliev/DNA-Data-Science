import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
from io import StringIO

# --- 1. Utility Functions ---

def get_gc_content(kmer):
    """Returns the GC percentage of a given k-mer."""
    return gc_fraction(kmer) * 100

# --- 2. Visualization ---

def save_visual_map(seq_len, results, file_name, results_folder):
    """Generates and saves a scatter plot showing marker distribution."""
    plt.figure(figsize=(12, 3))
    plt.hlines(1, 0, seq_len, colors='gray', lw=2, alpha=0.5)
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    
    # Map top 5 markers
    for i, (marker, count, gc, pos_list) in enumerate(results[:5]):
        plt.scatter(pos_list, [1] * len(pos_list), color=colors[i % len(colors)], 
                    label=f"{marker} (x{count})", s=70, edgecolors='black', zorder=3)
    
    plt.ylim(0.8, 1.2)
    plt.yticks([])
    plt.title(f"Genomic Marker Map: {file_name}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    output_path = os.path.join(results_folder, file_name.replace(".txt", ".png"))
    plt.savefig(output_path)
    plt.close()

# --- 3. Ultra-Fast Analytical Engine ---

def run_full_analysis(data_folder, results_folder, k=7, max_err=1):
    """Professional-grade speed using Neighborhood Search and Memory Caching."""
    files = [f for f in os.listdir(data_folder) if f.endswith('.txt')]
    report_path = os.path.join(results_folder, "summary_report.txt")
    
    # Pre-load all files into memory
    all_files_content = {}
    for f_name in files:
        with open(os.path.join(data_folder, f_name), 'r', errors='ignore') as f:
            all_files_content[f_name] = f.read().upper()

    with open(report_path, "w", encoding="utf-8") as report:
        report.write("ULTRA-FAST BIOINFORMATICS REPORT (Neighborhood Search)\n")
        report.write("-" * 50 + "\n\n")

        print(f"🚀 Speed-of-light analysis started for {len(files)} files...")

        for f_name, raw_content in all_files_content.items():
            try:
                record = SeqIO.read(StringIO(raw_content), "fasta-pearson")
                sequence = str(record.seq).upper()
            except:
                continue

            # STEP 1: Fast Indexing
            kmer_index = {}
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                kmer_index.setdefault(kmer, []).append(i)

            # STEP 2: Pre-calculate other data pool for uniqueness check
            other_data_pool = "".join([c for name, c in all_files_content.items() if name != f_name])

            # STEP 3: Neighborhood Search (The "Secret Sauce" for speed)
            final_results = []
            bases = ['A', 'C', 'G', 'T']
            unique_patterns = list(kmer_index.keys())
            
            for target_kmer in unique_patterns:
                # Start with exact matches
                all_matching_positions = list(kmer_index[target_kmer])
                
                # Generate all possible neighbors with 1 mismatch
                if max_err == 1:
                    for pos in range(k):
                        for char in bases:
                            if char != target_kmer[pos]:
                                neighbor = target_kmer[:pos] + char + target_kmer[pos+1:]
                                if neighbor in kmer_index:
                                    all_matching_positions.extend(kmer_index[neighbor])
                
                # Filter results
                if len(all_matching_positions) >= 2:
                    gc = get_gc_content(target_kmer)
                    if 40 <= gc <= 60 and target_kmer not in other_data_pool:
                        final_results.append((target_kmer, len(all_matching_positions), gc, sorted(all_matching_positions)))

            # Sorting and reporting
            final_results.sort(key=lambda x: x[1], reverse=True)

            if final_results:
                print(f"✅ Instant Process: {f_name}")
                save_visual_map(len(sequence), final_results, f_name, results_folder)
                
                top = final_results[0]
                report.write(f"File: {f_name}\n")
                report.write(f"  - Markers found: {len(final_results)}\n")
                report.write(f"  - Best marker: {top[0]} (found {top[1]} times)\n\n")

    print(f"\n✨ Done! Analysis finished in the blink of an eye.")

# --- 4. Main Execution ---

if __name__ == "__main__":
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_PATH = os.path.normpath(os.path.join(BASE_DIR, "..", "data"))
    RESULTS_PATH = os.path.normpath(os.path.join(BASE_DIR, "..", "results"))

    # Cleanup results folder before starting
    import shutil
    if os.path.exists(RESULTS_PATH):
        shutil.rmtree(RESULTS_PATH)
    os.makedirs(RESULTS_PATH)

    run_full_analysis(DATA_PATH, RESULTS_PATH, k=7, max_err=1)