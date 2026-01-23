import os
from Bio.Seq import Seq

# Path configuration
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "revc.txt")

try:
    with open(input_path, "r") as file:
        dna_string = "".join(file.read().split()).upper()

    if dna_string:
        # Create a Seq object and use built-in reverse_complement method
        my_dna = Seq(dna_string)
        print(my_dna.reverse_complement())

except FileNotFoundError:
    print(f"Error: {input_path} not found.")