import os

# Parse list of cells into two groups with roughly the same number of cells per each batch

file_path = '/Users/megan/Projects/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/file_lists/MG_10kb_file_list.txt'  # Replace with the actual file path
output_file_A = '/Users/megan/Projects/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/file_lists/MG_10kb_Batch_A_file_list.txt'
output_file_B = '/Users/megan/Projects/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/file_lists/MG_10kb_Batch_B_file_list.txt'

filenames = []
filenames_A = []
filenames_B = []

with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        filenames.append(file_name)

batches = {}
for f in filenames:
    num, age = f.split('_')[0], f.split('_')[1]
    if (num, age) not in batches:
        batches[(num, age)] = 0
    batches[(num, age)] += 1

for batch in batches:
    batches[batch] = batches[batch] // 2

with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        num, age = file_name.split('_')[0], file_name.split('_')[1]
        # Check if the current file should be added to filenames_A or filenames_B
        if batches[(num, age)] > 0:
            filenames_A.append(file_name)
            batches[(num, age)] -= 1
        else:
            if len(filenames_B)-len(filenames_A) == 1:
                filenames_A.append(file_name)
            else:
                filenames_B.append(file_name)

with open(output_file_A, "w") as file:
    for item in filenames_A:
        file.write(item + '\n')
    file.flush()

with open(output_file_B, "w") as file:
    for item in filenames_B:
        file.write(item + '\n')
    file.flush()