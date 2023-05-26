import os

def create_file_list(dir_a, dir_b, out_dir):
    file_list = [('name','group')]

    # Process files from directory A
    for file_name in os.listdir(dir_a):
        file_path = os.path.join(dir_a, file_name)
        if os.path.isfile(file_path):
            file_list.append((file_name, 'A'))

    # Process files from directory B
    for file_name in os.listdir(dir_b):
        file_path = os.path.join(dir_b, file_name)
        if os.path.isfile(file_path):
            file_list.append((file_name, 'B'))

    # Write file_list to file
    with open(out_dir+'file_list.txt', 'w') as file:
        for item in file_list:
            file.write(','.join(item) + '\n')

    print("File list created successfully.")

# Example usage
directory_a = '/project/compbio-lab/scHi-C/SnapHiC-D_example_data/sox2_ESC/'
directory_b = '/project/compbio-lab/scHi-C/SnapHiC-D_example_data/sox2_NPC/'
output_dir = '/project/compbio-lab/scHi-C/SnapHiC-D_example_data/sox2_outdir/'
create_file_list(directory_a, directory_b, output_dir)