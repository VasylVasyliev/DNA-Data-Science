import os

def solve_hamm():
    """
    Calculates the Hamming distance between two DNA strings.
    The Hamming distance is the number of positions at which 
    the corresponding symbols are different.
    """
    input_path = "data/hamm.txt"
    
    # 1. Check if the data file exists
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found. Please create it and add two DNA strings.")
        return

    # 2. Read the sequences
    with open(input_path, "r") as f:
        # We use .strip() to remove any trailing newlines or spaces
        s = f.readline().strip()
        t = f.readline().strip()

    # 3. Safety check: Hamming distance is only defined for strings of equal length
    if len(s) != len(t):
        print("Error: Sequences must be of equal length!")
        return

    # 4. Optimized calculation using a generator expression
    # zip(s, t) pairs up characters: (s0, t0), (s1, t1)...
    # 'a != b' returns True (which equals 1) or False (which equals 0)
    distance = sum(1 for a, b in zip(s, t) if a != b)
    
    # 5. Output the result
    print(distance)

if __name__ == "__main__":
    solve_hamm()