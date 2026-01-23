import os

# Set up paths to the data folder
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "subs.txt")

def find_motif_locations():
    try:
        with open(input_path, "r") as f:
            lines = f.read().splitlines()
            if len(lines) < 2:
                print("Error: File must contain two lines (DNA and Motif).")
                return
            
            dna, motif = lines[0].strip(), lines[1].strip()
        
        positions = []
        # Sliding window approach: check every possible substring of motif length
        for i in range(len(dna) - len(motif) + 1):
            if dna[i : i + len(motif)] == motif:
                # Biological indexing starts at 1, so we add 1 to the index
                positions.append(str(i + 1))
                
        # Join positions with spaces as required by Rosalind
        print(" ".join(positions))

    except FileNotFoundError:
        print(f"Error: Could not find file at {input_path}")

if __name__ == "__main__":
    find_motif_locations()