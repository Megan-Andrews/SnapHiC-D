import cooler
import os
from hicrep.utils import readMcool
from hicrep import hicrepSCC
import random


# cooler.merge_coolers(
# bulk in 
#pseudobulk sizes: 1, 2, 4, 8, 16, 32, 64, 128
#experiment sizes: 128, 64, 32, 16, 8, 4, 2, 1

bulk_path = "/project/compbio-lab/Hi-C"
cells =[ "GM12878" , "H1-hESC",  "HFF",  "IMR90"  ]
cell_files = ["GM12878_files.txt", "H1hESC_files.txt", "HFF_files.txt", "IMR90_files.txt"]
# Chr 21:

experiment_sizes = [128, 112, 96, 80, 64, 48, 32, 16, 8, 4, 2]
out_dir = "/project/compbio-lab/scHi-C/Kim2020/pseudobulk"

Kim2020_cools = "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool"
rwr_cools = "/project/compbio-lab/scHi-C/Kim2020/kim2020_cool_rwr3"
scVI_cools = "/project/compbio-lab/scHi-C/Kim2020/results/2023-09-14/Kim2020_scVI_model/imputed_cools"
higashi_k0_cools = "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool_higashi_k0"
higashi_k5_cools =  "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool_higashi_k5"

Kim2020_cools_out = os.path.join(out_dir, "Kim2020_cools")
rwr_cools_out = os.path.join(out_dir, "rwr_cools")
scVI_cools_out = os.path.join(out_dir, "scVI_cools")
higashi_k0_cools_out = os.path.join(out_dir, "higashi_k0_cools")
higashi_k5_cools_out = os.path.join(out_dir, "higashi_k5_cools")
conditions = [(Kim2020_cools, Kim2020_cools_out), (rwr_cools, rwr_cools_out), (scVI_cools, scVI_cools_out), (higashi_k0_cools, higashi_k0_cools_out), (higashi_k5_cools, higashi_k5_cools_out)]

for in_dir, out_dir in conditions[4:]:
    for cell_file in cell_files:
        with open(cell_file, 'r') as file:
            file_paths = file.read().splitlines()
        #print(cooler.Cooler(os.path.join(in_dir,file_paths[0])).info)
        for e in experiment_sizes:
            if len(file_paths) >= e:               
                selected_files = random.sample(file_paths, e)
                if in_dir == rwr_cools:
                    selected_files = [os.path.join(in_dir, file.replace(".cool","_rwr.cool")) for file in selected_files]
                else:
                    selected_files = [os.path.join(in_dir, file) for file in selected_files]
                out_path = os.path.join(out_dir, f"{os.path.basename(cell_file).replace('_files.txt', '')}_{e}.cool")
                cooler.merge_coolers(out_path,selected_files,100000) #,agg={'columns': 'mean'})
                print(f"merged {e} files from {cell_file} to {out_path}")
print("finished merging")
