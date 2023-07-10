import pandas as pd
import numpy as np  

file_path = '/Users/megan/Projects/SnapHiC-D/ext/Lee_MG_samples.txt'  # Replace with the actual file path
output_file = '/Users/megan/Projects/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/file_lists/MG_10kb_file_list.txt'

filenames = []

with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        file_name = line.strip()  # Remove leading/trailing whitespace and newline characters
        filenames.append(file_name)

metadata = pd.read_csv('/Users/megan/Projects/SnapHiC-D/ext/sn-m3c-seq-metadata.txt', sep='\t')
id2ct = {metadata['id'][i]: metadata['cell_type'][i] for i in range(metadata.shape[0])}

def get_id(filename):
    base = filename.split('.')[0]
    id = np.array(base.split('_'))[[1,4,9,10,11]]
    id = '_'.join(id)
    return id

def get_raw_coolname(id,ct):
    return '{}_{}_10kb_contacts.cool'.format(id,ct)

file_ids = [get_id(f) for f in filenames]
file_cts = [id2ct[id] for id in file_ids]
coolnames = [get_raw_coolname(id,ct) for id, ct in zip(file_ids,file_cts)]

with open(output_file, "w") as file:
    for item in coolnames:
        file.write(item + '\n')
    file.flush()