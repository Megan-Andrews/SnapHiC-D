import cooler


file_path = '/Users/megan/Projects/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/file_lists/MG_10kb_Batch_A_file_list.txt'
output_path = 'scHi-C_link/Lee2019/pseudo-bulk_data/MG/MG_10kb_Batch_A.cool'

filenames = []

with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        filenames.append(file_name)


cooler.merge_coolers(output_path, filenames, mergebuf=1000000)