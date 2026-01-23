import os

# Set up the paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "dna.txt")

try:
    with open(input_path, "r") as file:
        # Read and remove all whitespace characters
        dna = "".join(file.read().split()).upper()

    if not dna:
        print("Error: The data file is empty.")
    else:
        # Print sequence length for verification
        print(f"Sequence length: {len(dna)}") 

        # Count occurrences of each nucleotide in the specified order
        nucleotides = "ACGT"
        counts = [str(dna.count(n)) for n in nucleotides]

        # Join counts into a single space-separated string
        print(" ".join(counts))
        # Save results for R in the data folder
        # Сохраняем результаты для R в папку data
        output_path = os.path.join(script_dir, "..", "data", "nucleotide_counts.csv")
        with open(output_path, "w") as out_file:
            out_file.write("Nucleotide,Count\n")
            for n in nucleotides:
                out_file.write(f"{n},{dna.count(n)}\n")
        
        print(f"Results saved to {output_path}")
except FileNotFoundError:
    print(f"Error: File not found at {input_path}")