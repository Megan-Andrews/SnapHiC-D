import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
Kim2020_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/Kim2020_cools"
rwr_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/rwr_cools"
scVI_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/scVI_cools"
higashi_k0_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/higashi_k0_cools"
higashi_k5_cools =  "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/higashi_k5_cools"
dirs = [Kim2020_cools, rwr_cools, scVI_cools, higashi_k0_cools, higashi_k5_cools]
#names = ["raw_SCC", "rwr_SCC","scVI_SCC","higashi_k0_SCC","higashi_k5_SCC"]

GM12878_HFF_pairs = "h0_GM12878_HFF_selected_pairs.txt"
GM12878_pairs = "h0_GM12878_selected_pairs.txt"
HFF_pairs =  "h0_HFF_selected_pairs.txt"
pairs_list = [GM12878_HFF_pairs, GM12878_pairs, HFF_pairs]
names = ["GM12878_HFF_SCC", "GM12878_SCC", "HFF_SCC"]

for i, d in enumerate(dirs):
    all_violin_data = pd.read_csv(os.path.join(d, pairs_list[0]), sep='\t')
    del all_violin_data["HiCRep_SCC"]
    for j, p in enumerate(pairs_list):
        temp_df = pd.read_csv(os.path.join(d, p), sep='\t')
        temp_df["HiCRep_SCC"] = pd.to_numeric(temp_df["HiCRep_SCC"])
        all_violin_data[names[j]] = temp_df["HiCRep_SCC"]
    print(all_violin_data)
    csv_filename = os.path.join(d, "h0_all_pairs.csv")  # Specify the desired output filename
    all_violin_data.to_csv(csv_filename, index=False)


