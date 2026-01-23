from Bio import SeqIO

def solve_cons():
    # 1. Read FASTA file - SeqIO handles multi-line sequences automatically
    records = list(SeqIO.parse("data/cons.txt", "fasta"))
    sequences = [str(record.seq) for record in records]
    
    if not sequences:
        print("Error: No sequences found in data/cons.txt")
        return

    seq_length = len(sequences[0])
    
    # 2. Initialize the Profile Matrix in the exact order: A, C, G, T
    profile = {base: [0] * seq_length for base in "ACGT"}
    
    # Fill the matrix
    for seq in sequences:
        for i, base in enumerate(seq):
            profile[base][i] += 1
            
    # 3. Create the Consensus String
    consensus = []
    for i in range(seq_length):
        # Find the most frequent base at each position
        max_base = max("ACGT", key=lambda b: profile[b][i])
        consensus.append(max_base)
        
    # 4. Final Output (Rosalind format)
    print("".join(consensus))
    for base in "ACGT":
        counts = " ".join(map(str, profile[base]))
        print(f"{base}: {counts}")

if __name__ == "__main__":
    solve_cons()