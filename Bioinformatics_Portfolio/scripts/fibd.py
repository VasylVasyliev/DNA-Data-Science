import os

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "fibd.txt")

def solve_mortal_rabbits():
    try:
        with open(input_path, "r") as file:
            line = file.read().strip()
            if not line:
                return
            n, m = map(int, line.split())

        # Each index i represents rabbits that are i+1 months old
        # ages[0] are newborns, ages[m-1] are those about to die
        ages = [0] * m
        ages[0] = 1 # Start with one newborn pair

        for month in range(1, n):
            # 1. New rabbits are born from everyone except newborns
            new_borns = sum(ages[1:])
            
            # 2. Shift all age groups (rabbits get older)
            # Using list slicing to shift: [newborns, age1, age2, ..., age_m-1]
            ages = [new_borns] + ages[:-1]

        print(sum(ages))

    except FileNotFoundError:
        print(f"Error: Could not find file at {input_path}")

if __name__ == "__main__":
    solve_mortal_rabbits()