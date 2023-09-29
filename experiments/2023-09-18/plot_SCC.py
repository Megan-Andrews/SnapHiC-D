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

#plt.figure(figsize=(15, 5))  

violins = []
fig, axes = plt.subplots(1, 1, figsize=(15, 5))

# Set the common x-axis ticks and labels
x_ticks = [1, 2, 3]
x_tick_labels = ["GM12878-HFF", "GM12878-GM12878", "HFF-HFF"]

for j, p in enumerate(pairs_list):
    ax = axes  # Get the current axis
    
    for i, d in enumerate(dirs):
        temp_df = pd.read_csv(os.path.join(d, p), sep='\t')
        temp_df["HiCRep_SCC"] = pd.to_numeric(temp_df["HiCRep_SCC"])
        
        v = ax.violinplot(temp_df["HiCRep_SCC"], showmeans=True, showmedians=True)
        violins.append(v['bodies'][0])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels)
        ax.set_xlabel('Groups')
        ax.set_ylabel('Similarity Scores (SCC)')
        ax.set_title(p.replace(".txt", ""))
        ax.grid(axis='y')

# Add a legend to the first subplot (you can customize this as needed)
axes.legend(violins, ["Raw", "RWR", "scVI", "Higashi K0", "Higashi K5"], loc='upper left', bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig(os.path.join("/project/compbio-lab/scHi-C/Kim2020/similarity_scores/SCC_plot_2.png"))
plt.savefig(os.path.join("/home/maa160/SnapHiC-D/experiments/2023-09-18/SCC_plot_2.png"))

