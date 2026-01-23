from Bio import SeqIO
from collections import defaultdict

def solve_grph():
    """
    Identifies all edges in an overlap graph for a collection of DNA strings.
    An edge exists from string s to string t if the suffix of s (length k) 
    matches the prefix of t (length k).
    """
    # 1. Parse the FASTA file using Biopython
    # We use list() to keep all records in memory for comparison
    records = list(SeqIO.parse("data/grph.txt", "fasta"))
    k = 3
    
    # 2. Build a prefix index (Hash Map) for optimization
    # Instead of O(n^2) nested loops, we map each 3-mer prefix 
    # to the IDs of sequences that start with it.
    prefix_map = defaultdict(list)
    for rec in records:
        prefix = str(rec.seq[:k])
        prefix_map[prefix].append(rec.id)
    
    # 3. Traverse the records and find matching prefixes for each suffix
    for s_rec in records:
        suffix = str(s_rec.seq[-k:])
        
        # Look up the suffix in our pre-built prefix index
        if suffix in prefix_map:
            for t_id in prefix_map[suffix]:
                # Biology rule: A sequence cannot overlap with itself
                if s_rec.id != t_id:
                    # Output the directed edge as required by Rosalind
                    print(f"{s_rec.id} {t_id}")

if __name__ == "__main__":
    solve_grph()