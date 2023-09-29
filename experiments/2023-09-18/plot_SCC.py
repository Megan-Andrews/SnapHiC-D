import pandas as pd
import os
import matplotlib.pyplot as plt

Kim2020_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/Kim2020_cools"
rwr_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/rwr_cools"
scVI_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/scVI_cools"
higashi_k0_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/higashi_k0_cools"
higashi_k5_cools =  "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/higashi_k5_cools"
dirs = [Kim2020_cools, rwr_cools, scVI_cools, higashi_k0_cools, higashi_k5_cools]

GM12878_HFF_pairs = "GM12878_HFF_selected_pairs.txt"
GM12878_pairs = "GM12878_selected_pairs.txt"
HFF_pairs =  "HFF_selected_pairs.txt"
pairs_list = [GM12878_HFF_pairs, GM12878_pairs, HFF_pairs]

plt.figure(figsize=(15, 5))  

violins = []
for j, p in enumerate(pairs_list):
    plt.subplot(131+j)
    plt.xticks([1,2,3,4,5], ["Raw","RWR", "scVI","Higashi K0","Higashi K5"])
    for i, d in enumerate(dirs):
        temp_df  = pd.read_csv(os.path.join(d,p), sep='\t')
        # temp_df.columns = ["cool1", "cool2", "HiCRep_SCC"]
        #subplot_index = len(dirs) * j + i + 1
        # Create subplot with four columns
        #plt.subplot(len(pairs_list), len(dirs), subplot_index)
        temp_df["HiCRep_SCC"] = pd.to_numeric(temp_df["HiCRep_SCC"])
        v = plt.violinplot(temp_df["HiCRep_SCC"], showmeans=True, showmedians=True)
        violins.append(v['bodies'][0])
        plt.xlabel('Groups')
        plt.ylabel('Similarity Scores (SCC)')
	# plt.legend()
        plt.title(p.replace(".txt", ""))
        plt.grid(axis='y')
    plt.legend(violins,  ["Raw","scVI","Higashi K0","Higashi K5"])
    plt.tight_layout()
plt.savefig(os.path.join("/project/compbio-lab/scHi-C/Kim2020/similarity_scores/SCC_plot_2.png"))
plt.savefig(os.path.join("/home/maa160/SnapHiC-D/experiments/2023-09-18/SCC_plot_2.png"))
