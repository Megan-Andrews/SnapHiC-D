import cooler 
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import cooler

chrom_size_path = "../../ext/hg19.chrom.sizes.ordered"

def visualize_hic(coolfile, matrixfile, outputfile): 
    c = cooler.Cooler(coolfile)

    # Get the contact matrix
    chrom = "chr22"
    cool_matrix = c.matrix(balance=False).fetch(f"{chrom}")

    columns = ["binId_1", "binId_2", "counts","normalized_counts", "1st_chr", "2nd_chr"]
    matrix_df = pd.read_csv(matrixfile, delimiter='\t', header=None)
    matrix_df.columns = columns
    matrix_df = matrix_df[matrix_df["1st_chr"] == matrix_df["2nd_chr"]] # only include intra-chromosome pairs
    matrix_df = matrix_df[matrix_df["1st_chr"] == "human_chr22"] # only include intra-chromosome pairs

    def get_chrom_offsets(bins_df):
        chrom_offset = {chrom: bins_df[bins_df['chrom'] == chrom].index[0]
                        for chrom in bins_df['chrom'].cat.categories}
        return chrom_offset

    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t',index_col=0, header=None).squeeze(axis=1)
    bins_df = cooler.binnify(chrom_sizes, 500000)
    chrom_offsets = get_chrom_offsets(bins_df)

    mask = matrix_df["1st_chr"] == "human_chr22"
        
    # Update binId_1 and binId_2 using .loc
    matrix_df.loc[mask, "binId_1"] -= chrom_offsets["chr22"]
    matrix_df.loc[mask, "binId_2"] -= chrom_offsets["chr22"]

    matrix_df = matrix_df[['binId_1', 'binId_2', 'counts']]
    matrix_reflection = matrix_df.copy()
    matrix_reflection.columns = ['binId_2', 'binId_1', 'counts']
    additional_row = pd.DataFrame([[0, 0, 0]])
    additional_row.columns = ['binId_1', 'binId_2', 'counts']
    matrix_df = pd.concat([matrix_df, matrix_reflection, additional_row], ignore_index=True)
    matrix_df = matrix_df.drop_duplicates()
    pivot_df = matrix_df.pivot(index='binId_1', columns='binId_2', values='counts')
    pivot_df = pivot_df.fillna(0)

    # Clear the plot
    plt.clf()

     # Create a subplot with 1 row and 2 columns
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    im1 = ax1.imshow(cool_matrix, cmap='coolwarm', origin="lower")
    ax1.set_title(".cool")
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Genomic Position")
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.set_label('Intensity')

    ax2 = axes[1]
    im2 = ax2.imshow(pivot_df, cmap='coolwarm', origin="lower")
    ax2.set_title(".matrix")
    ax2.set_xlabel("Genomic Position")
    ax2.set_ylabel("Genomic Position")
    cbar2 = fig.colorbar(im2, ax=ax2)
    cbar2.set_label('Intensity')

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure
    plt.savefig(outputfile)

# human_9963_HFF_H1Esc-HFF.R1.cool
visualize_hic("/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/human_9963_HFF_H1Esc-HFF.R1.cool",
              "/project/compbio-lab/scHi-C/Kim2020/H1Esc-HFF.R1/human_9963_TGGAGAGG-GCTAACGA_500000.matrix",
              "/home/maa160/SnapHiC-D/experiments/2023-08-21/test_visualizations/human_9963_HFF_H1Esc-HFF.R1.png")

# human_9996_H1Esc_H1Esc.R2.cool
visualize_hic("/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/human_9996_H1Esc_H1Esc.R2.cool",
              "/project/compbio-lab/scHi-C/Kim2020/H1Esc.R2/human_9996_AAGCCGGT-CTTGGTTA_500000.matrix",
              "/home/maa160/SnapHiC-D/experiments/2023-08-21/test_visualizations/human_9996_H1Esc_H1Esc.R2.png")
