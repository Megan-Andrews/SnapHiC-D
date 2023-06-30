import os
import shutil

def move_files(source_dir, destination_dir, file_list):
    for file_name in file_list:
        source_path = os.path.join(source_dir, file_name)
        destination_path = os.path.join(destination_dir, file_name)

        if os.path.isfile(source_path):
            shutil.move(source_path, destination_path)
            print(f"Moved file: {file_name}")
        else:
            print(f"File not found: {file_name}")

source_directory = "/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/ODC_1_100kb_imputed_cool"
destination_directory = "/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/ODC_2_100kb_imputed_cool"
file_path = "/home/maa160/SnapHiC-D/experiments/2023-06-30/file_lists/ODC_file_listB_batches.txt"

with open(file_path, "r") as file:
    for line in file:
        line = line.strip() 
        line = line.replace("_10kb_", "_100kb_").replace(".cool", "_imputed.cool")
        file_list = [line]  # Create a list containing the file name
        move_files(source_directory, destination_directory, file_list)