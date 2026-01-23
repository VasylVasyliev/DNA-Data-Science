import os
from Bio.Seq import Seq

# Path configuration
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "prot.txt")

try:
    with open(input_path, "r") as file:
        rna_string = "".join(file.read().split()).upper()

    if rna_string:
        # Translate the RNA sequence directly to protein
        # to_stop=True ensures it stops at the first stop codon
        messenger_rna = Seq(rna_string)
        protein = messenger_rna.translate(to_stop=True)
        print(protein)

except FileNotFoundError:
    print(f"Error: {input_path} not found.")
