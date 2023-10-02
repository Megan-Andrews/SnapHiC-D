import pandas as pd
import os
import matplotlib.pyplot as plt

Kim2020_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/Kim2020_cools"
rwr_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/rwr_cools"
scVI_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/scVI_cools"
higashi_k0_cools = "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/higashi_k0_cools"
higashi_k5_cools =  "/project/compbio-lab/scHi-C/Kim2020/similarity_scores/higashi_k5_cools"
dirs = [Kim2020_cools, rwr_cools, scVI_cools, higashi_k0_cools, higashi_k5_cools]

GM12878_HFF_pairs = "h0_GM12878_HFF_selected_pairs.txt"
GM12878_pairs = "h0_GM12878_selected_pairs.txt"
HFF_pairs =  "h0_HFF_selected_pairs.txt"
pairs_list = [GM12878_HFF_pairs, GM12878_pairs, HFF_pairs]

#plt.figure(figsize=(15, 5))  

violins = []
fig, axes = plt.subplots(1, 1, figsize=(8, 5))

# Set the common x-axis ticks and labels
x_ticks = [1, 2, 3]
x_tick_labels = ["GM12878-HFF", "GM12878-GM12878", "HFF-HFF"]
colors = ['blue', 'green', 'red', 'cyan', 'orange']
for j, p in enumerate(pairs_list):
    ax = axes  # Get the current axis
    all_violin_data = []
    for i, d in enumerate(dirs):
        temp_df = pd.read_csv(os.path.join(d, p), sep='\t')
        temp_df["HiCRep_SCC"] = pd.to_numeric(temp_df["HiCRep_SCC"])
        all_violin_data.append(temp_df["HiCRep_SCC"])

    pos = [j+1] * 5
    v = ax.violinplot(all_violin_data, positions=pos, showmeans=False, showmedians=False)                
        # Set colors for the violin bodies
    for violin, color in zip(v['bodies'], colors):
        violin.set_facecolor(color)
    
    # Set colors for the mean line and median line
#    for partname, color in zip(('cbars', 'cmedians'), colors):
#        part = v[partname]
#        part.set_edgecolor(color)
    violins.append(v['bodies'][0])
    violins.append(v['bodies'][1])
    violins.append(v['bodies'][2])
    violins.append(v['bodies'][3])
    violins.append(v['bodies'][4])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_xlabel('Groups')
    ax.set_ylabel('Similarity Scores (SCC)')
    ax.set_title("SCC")
    ax.grid(axis='y')

# Add a legend to the first subplot (you can customize this as needed)
axes.legend(violins, ["Raw", "RWR", "scVI", "Higashi K0", "Higashi K5"], loc='upper left', bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig(os.path.join("/project/compbio-lab/scHi-C/Kim2020/similarity_scores/SCC_plot_h0.png"))
plt.savefig(os.path.join("/home/maa160/SnapHiC-D/experiments/2023-09-18/SCC_plot_h0.png"))
plt.clf()

## Other plot

violins = []
fig, axes = plt.subplots(1, 5, figsize=(15, 5))

# Set the common x-axis ticks and labels
x_ticks = [1, 2, 3]
x_tick_labels = ["Both", "GM12878", "HFF"]
colors = ['blue', 'green', 'red', 'cyan', 'orange']
for i, d in enumerate(dirs):
    ax = axes[i]  # Get the current axis
    all_violin_data = []
    for j, p in enumerate(pairs_list):
        temp_df = pd.read_csv(os.path.join(d, p), sep='\t')
        temp_df["HiCRep_SCC"] = pd.to_numeric(temp_df["HiCRep_SCC"])
        all_violin_data.append(temp_df["HiCRep_SCC"])
    print(len(all_violin_data))
    # pos = [j+1] * 5
    v = ax.violinplot(all_violin_data, positions=[1,2,3], showmeans=False, showmedians=False)                
        # Set colors for the violin bodies
#    for violin, color in zip(v['bodies'], colors):
#        violin.set_facecolor(color)
    
    # Set colors for the mean line and median line
#    for partname, color in zip(('cbars', 'cmedians'), colors):
#        part = v[partname]
#        part.set_edgecolor(color)
#    violins.append(v['bodies'][0])
#    violins.append(v['bodies'][1])
#    violins.append(v['bodies'][2])
#    violins.append(v['bodies'][3])
#    violins.append(v['bodies'][4])
    ax.set_ylim(-1, 1)  # Set y-axis limits to [-1, 3]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_xlabel('Groups')
    ax.set_ylabel('Similarity Scores (SCC)')
    ax.set_title(os.path.basename(d))
    ax.grid(axis='y')

# Add a legend to the first subplot (you can customize this as needed)
#axes.legend(violins, ["Raw", "RWR", "scVI", "Higashi K0", "Higashi K5"], loc='upper left', bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig(os.path.join("/project/compbio-lab/scHi-C/Kim2020/similarity_scores/SCC_plot_h0_subplots.png"))
plt.savefig(os.path.join("/home/maa160/SnapHiC-D/experiments/2023-09-18/SCC_plot_h0_subplots.png"))



