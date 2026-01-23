import os

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(script_dir, "..", "data", "fib.txt")

def calculate_rabbits():
    try:
        with open(input_path, "r") as file:
            line = file.read().strip()
            if not line:
                print("Error: The file is empty.")
                return
            
            # n: number of months, k: offspring per litter
            n, k = map(int, line.split())

        # Base cases: for month 1 and month 2, we have 1 pair
        # we use variables to store the population of the previous two months
        parent, child = 1, 1
        
        # We start calculating from the 3rd month up to n
        for _ in range(n - 2):
            # The new population is:
            # All rabbits from previous month (parent) 
            # PLUS the newborns (child * k)
            child, parent = parent, parent + (child * k)
        
        print(parent)

    except FileNotFoundError:
        print(f"Error: Could not find file at {input_path}")
    except ValueError:
        print("Error: Could not parse integers from the file.")

if __name__ == "__main__":
    calculate_rabbits()
