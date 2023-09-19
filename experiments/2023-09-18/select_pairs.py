import random
import itertools

# Function to randomly select n unique pairs of files from two lists
def select_random_unique_pairs(list1, list2, n):
    unique_pairs = set()
    while len(unique_pairs) < n:
        c1 = random.choice(list1)
        c2 = random.choice(list2)
        if c1 == c2 or (random.choice(list1), random.choice(list2)) in unique_pairs:
            continue
        pair = (random.choice(list1), random.choice(list2))
        unique_pairs.add(pair)
    return list(unique_pairs)

# Function to write a list of pairs to a new file
def write_file_pairs(pairs, output_file):
    with open(output_file, 'w') as file:
        for pair in pairs:
            file.write(f"{pair[0]}\t{pair[1]}\n")

# Input file paths
input_file1 = 'HFF_files.txt'  # Replace with your first input file path
input_file2 = 'GM12878_files.txt'  # Replace with your second input file path



# Read the input file paths into lists
with open(input_file1, 'r') as file:
    file_paths1 = file.read().splitlines()

with open(input_file2, 'r') as file:
    file_paths2 = file.read().splitlines()

# Randomly select 100 unique pairs of files
random_pairs = select_random_unique_pairs(file_paths1, file_paths2, 100)

# Write the selected pairs to the output file
write_file_pairs(random_pairs, 'GM12878_HFF_selected_pairs.txt')

# Randomly select 100 unique pairs of files
random_pairs = select_random_unique_pairs(file_paths1, file_paths1, 100)

# Write the selected pairs to the output file
write_file_pairs(random_pairs, 'HFF_selected_pairs.txt')

# Randomly select 100 unique pairs of files
random_pairs = select_random_unique_pairs(file_paths2, file_paths2, 100)

# Write the selected pairs to the output file
write_file_pairs(random_pairs, "GM12878_selected_pairs.txt")


print("Random unique pair selection and writing completed.")
