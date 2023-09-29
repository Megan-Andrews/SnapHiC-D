

import hicrep
import cooler 
from hicrep.utils import readMcool
from hicrep import hicrepSCC
import os 
import pandas as pd 
import numpy as np

cell1 = 'GM12878' 
cell2 = 'HFF' #from Kim2020


GM12878_HFF_pairs = "/home/maa160/SnapHiC-D/experiments/2023-09-18/GM12878_HFF_selected_pairs.txt"
GM12878_pairs = "/home/maa160/SnapHiC-D/experiments/2023-09-18/GM12878_selected_pairs.txt"
HFF_pairs =  "/home/maa160/SnapHiC-D/experiments/2023-09-18/HFF_selected_pairs.txt"
pairs_list = [GM12878_HFF_pairs, GM12878_pairs, HFF_pairs]

out_dir = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores"

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
#conditions = [(Kim2020_cools, Kim2020_cools_out),(scVI_cools, scVI_cools_out), (higashi_k0_cools, higashi_k0_cools_out), (higashi_k5_cools, higashi_k5_cools_out)]

for in_dir, o_dir in conditions:
    for file in pairs_list:
        print(in_dir, file)
        df = pd.read_csv(file, sep='\t', header=None, names=["cool1", "cool2"])

        def calc_similarity(row):
            if in_dir != rwr_cools:
                cool1 = os.path.join(in_dir, row['cool1'])
                cool2 = os.path.join(in_dir, row['cool2'])
            else:
                cool1 = os.path.join(in_dir, row['cool1'].replace(".cool", "_rwr.cool"))
                cool2 = os.path.join(in_dir, row['cool2'].replace(".cool", "_rwr.cool"))

            cool1, binSize1 = readMcool(cool1, -1)
            cool2, binSize2 = readMcool(cool2, -1)
            binSize = binSize1
            # smoothing window half-size
            h = 1
            # maximal genomic distance to include in the calculation
            dBPMax = -1
            # whether to perform down-sampling or not 
            # if set True, it will bootstrap the data set # with larger contact counts to
            # the same number of contacts as in the other data set; otherwise, the contact 
            # matrices will be normalized by the respective total number of contacts
            bDownSample = False
            # compute the SCC score
            # this will result in a SCC score for each chromosome available in the data set
            # listed in the same order as the chromosomes are listed in the input Cooler files
            scc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, np.array(['chr22',], dtype=str))
            return scc[0]
        
        df['HiCRep_SCC'] = df.apply(calc_similarity, axis=1)
        print(os.path.join(o_dir, os.path.basename(file)))
        df.to_csv(os.path.join(o_dir, os.path.basename(file)), sep='\t', index=False)





       
