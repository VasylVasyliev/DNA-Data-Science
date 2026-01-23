import os

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "gc.txt")

def calculate_gc(sequence):
    """Calculates GC content percentage of a DNA sequence."""
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return (g_count + c_count) / len(sequence) * 100

try:
    with open(input_path, "r") as file:
        # Dictionary to store {ID: sequence}
        sequences = {}
        current_id = ""

        # Parse FASTA format
        for line in file:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                # Remove '>' and set as current ID
                current_id = line[1:]
                sequences[current_id] = ""
            else:
                # Append DNA string to the current ID
                sequences[current_id] += line

    # Find the sequence with the highest GC content
    max_id = ""
    max_gc = 0.0

    for seq_id, dna in sequences.items():
        gc_val = calculate_gc(dna)
        if gc_val > max_gc:
            max_gc = gc_val
            max_id = seq_id

    # Print the result as required by Rosalind
    print(max_id)
    print(f"{max_gc:.6f}")

except FileNotFoundError:
    print(f"Error: {input_path} not found.")