
import os
import cooler
import time

def main():
    unprocessed_file_list = "/home/maa160/SnapHiC-D/file_lists/MG_150kcounts_file_list.txt"
    inputdir = "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_cool/"
    cell_type = "_MG_"
    file_list =[]
    count = 0
    for file_name in os.listdir(inputdir):
            file_path = os.path.join(inputdir, file_name)
            if os.path.isfile(file_path) and file_name.find(cell_type) != -1:
                i = cooler.Cooler(file_path).info['sum'] # the sum of the nnz matrix cells
                if i >= 150000:
                    count += 1
                    print(count, file_name)
                    file_list.append(file_name)

    with open(unprocessed_file_list, "w") as file:
        for item in file_list:
            file.write(item + '\n')
        file.flush()
        
        


if __name__ == "__main__":
    main()
