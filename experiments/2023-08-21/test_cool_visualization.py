import cooler 
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import cooler


def visualize_hic(coolfile, matrixfile, outputfile): 
    c = cooler.Cooler(coolfile)

    # Get the contact matrix
    chrom = "chr22"
    cool_matrix = c.matrix(balance=False).fetch(f"{chrom}")

    columns = ["binId_1", "binId_2", "counts","normalized_counts", "1st_chr", "2nd_chr"]
    matrix_df = pd.read_csv(matrixfile, delimiter='\t', header=None)
    matrix_df.columns = columns
    matrix_df = matrix_df[matrix_df["1st_chr"] == matrix_df["2nd_chr"]] # only include intra-chromosome pairs
    matrix_df = matrix_df[['binId_1', 'binId_2', 'counts']]
    pivot_df = matrix_df.pivot(index='binId_1', columns='binId_2', values='counts')
    pivot_df = pivot_df.fillna(0)

    # Clear the plot
    plt.clf()

     # Create a subplot with 1 row and 2 columns
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    im1 = ax1.imshow(cool_matrix, cmap='coolwarm', origin="lower")
    ax1.set_title("Cool Matrix")
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Genomic Position")
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.set_label('Intensity')

    ax2 = axes[1]
    im2 = ax2.imshow(pivot_df, cmap='coolwarm', origin="lower")
    ax2.set_title("Matrix DataFrame")
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
