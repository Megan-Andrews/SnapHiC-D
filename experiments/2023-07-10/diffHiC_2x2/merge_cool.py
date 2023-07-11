import cooler


file_path = '/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/file_lists/MG_10kb_Batch_A_file_list.txt'
directory = '/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_cool/'
output_path = '/project/compbio-lab/scHi-C/Lee2019/pseudo-bulk_data/MG/MG_10kb_Batch_A.cool'

filenames = []

with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        filenames.append(directory + file_name)


cooler.merge_coolers(output_path, filenames, mergebuf=1000000)