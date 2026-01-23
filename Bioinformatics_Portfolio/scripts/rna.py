import os

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "rna.txt")

try:
    with open(input_path, "r") as file:
        dna = "".join(file.read().split()).upper()

    if not dna:
        print("Error: File is empty.")
    else:
        # Replace Thymine with Uracil
        print(dna.replace("T", "U"))

except FileNotFoundError:
    print(f"Error: {input_path} not found.")
    