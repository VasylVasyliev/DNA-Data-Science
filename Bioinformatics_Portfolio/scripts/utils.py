import os
from Bio import SeqIO

def get_data_path(filename):
    """Returns the absolute path to a file in the data folder."""
    # This finds the project root (Bioinformatics_Portfolio) from the scripts folder
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_dir, "data", filename)

def read_fasta(filename):
    """Reads a FASTA file and returns the sequence object."""
    path = get_data_path(filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    return SeqIO.read(path, "fasta")