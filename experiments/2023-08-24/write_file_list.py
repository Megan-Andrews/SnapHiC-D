import os

# Directory path to read files from
directory_path = "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/"

# Get a list of all files in the directory
file_list = os.listdir(directory_path)

# Create or open the file_list.txt for writing
with open("/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/file_list.txt", "w") as file:
    # Write each file name to the file_list.txt, one per line
    for filename in file_list:
        file.write(directory_path +filename + "\n")