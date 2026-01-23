import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
# Import our new utility helper
from utils import get_data_path

def find_proteins(dna_sequence):
    """
    Scans a DNA sequence for all possible Open Reading Frames (ORFs).
    """
    proteins = []
    dna_len = len(dna_sequence)
    stop_codons = ["TAA", "TAG", "TGA"]
    
    for i in range(dna_len):
        if dna_sequence[i:i+3] == "ATG":
            for j in range(i + 3, dna_len - 2, 3):
                codon = dna_sequence[j:j+3]
                if codon in stop_codons:
                    dna_fragment = dna_sequence[i:j]
                    protein = dna_fragment.translate()
                    proteins.append(str(protein))
                    break 
    return proteins

def solve_orf(filename="phix174.txt"):
    # Use our utils to get the correct file path
    try:
        file_path = get_data_path(filename)
        record = SeqIO.read(file_path, "fasta")
        dna_seq = record.seq
        print(f"--- Analyzing: {record.description} ---")
        
        results = []
        results.extend(find_proteins(dna_seq))
        results.extend(find_proteins(dna_seq.reverse_complement()))
        
        unique_results = sorted(set(results), key=len, reverse=True)
        
        print(f"Found {len(unique_results)} potential proteins.")
        
        # --- NEW: Saving data for R visualization ---
        csv_path = get_data_path("protein_lengths.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["length"]) # CSV Header
            for p in unique_results:
                writer.writerow([len(p)])
        print(f"Data for visualization saved to {csv_path}")
        # --------------------------------------------

        print("\nTop 5 longest proteins:")
        for p in unique_results[:5]:
            print(f"Len {len(p)}: {p[:30]}...") 
            
    except Exception as e:
        print(f"Analysis failed: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        target_file = sys.argv[1]
    else:
        target_file = "phix174.txt"
    
    print(f"--- Analyzing file: {target_file} ---")
    solve_orf(target_file)