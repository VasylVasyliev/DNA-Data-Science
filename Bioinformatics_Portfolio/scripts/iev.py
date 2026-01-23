import os

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "iev.txt")

def calculate_expected_offspring():
    try:
        with open(input_path, "r") as file:
            line = file.read().strip()
            if not line:
                return
            # counts: numbers of couples for each of the 6 genotypes
            counts = list(map(int, line.split()))

        # Expected value multipliers for each couple (probability * 2 offspring)
        # 1: AA-AA -> 2
        # 2: AA-Aa -> 2
        # 3: AA-aa -> 2
        # 4: Aa-Aa -> 1.5
        # 5: Aa-aa -> 1
        # 6: aa-aa -> 0
        multipliers = [2, 2, 2, 1.5, 1, 0]
        
        expected_value = 0
        for i in range(len(counts)):
            expected_value += counts[i] * multipliers[i]
            
        print(expected_value)

    except FileNotFoundError:
        print(f"Error: Could not find file at {input_path}")

if __name__ == "__main__":
    calculate_expected_offspring()