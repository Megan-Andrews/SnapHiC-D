
import os
import time
from collections import defaultdict 

def main():
    unprocessed_file_list_A = "/home/maa160/SnapHiC-D/experiments/2023-06-30/file_lists/ODC_file_listA_batches.txt"
    unprocessed_file_list_B = "/home/maa160/SnapHiC-D/experiments/2023-06-30/file_lists/ODC_file_listB_batches.txt"
    inputdir = "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_cool/"
    cell_type = "_ODC_"

    batches = defaultdict(list)

    for file_name in os.listdir(inputdir):
        if file_name.find(cell_type) != -1:
            batch = tuple(file_name.split("_")[:2])
            batches[batch].append(file_name)

    file_list_A = []
    file_list_B = []

    for batch_files in batches.values():
        batch_size = len(batch_files)
        split_index = batch_size // 2
        file_list_A.extend(batch_files[:split_index])
        file_list_B.extend(batch_files[split_index:])

    with open(unprocessed_file_list_A, "w") as file_A, open(unprocessed_file_list_B, "w") as file_B:
        for item in file_list_A:
            file_A.write(item + '\n')
        for item in file_list_B:
            file_B.write(item + '\n')
        

if __name__ == "__main__":
    main()
