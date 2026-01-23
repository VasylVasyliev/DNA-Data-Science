import os

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "iprb.txt")

try:
    with open(input_path, "r") as file:
        line = file.read().strip()
        if not line:
            print("Error: File is empty.")
        else:
            # k: AA, m: Aa, n: aa
            k, m, n = map(int, line.split())
            total = k + m + n

            # Calculate the probability of picking two organisms that can produce 'aa'
            # We use the complement probability: 1 - P(recessive offspring)
            
            # Total combinations of any two organisms
            all_combos = total * (total - 1)
            
            # Combinations that can result in 'aa' offspring:
            # 1. Both are 'aa' (n, n): all offspring are 'aa' (100%)
            prob_nn = (n * (n - 1)) / all_combos
            
            # 2. Both are 'Aa' (m, m): 1/4 of offspring are 'aa' (25%)
            prob_mm = (m * (m - 1)) / all_combos * 0.25
            
            # 3. One 'Aa' and one 'aa' (m, n) or (n, m): 1/2 of offspring are 'aa' (50%)
            prob_mn = (m * n * 2) / all_combos * 0.5
            
            # Dominant probability is 1 minus the sum of recessive probabilities
            dominant_prob = 1 - (prob_nn + prob_mm + prob_mn)
            
            print(f"{dominant_prob:.5f}")

except FileNotFoundError:
    print(f"Error: {input_path} not found.")