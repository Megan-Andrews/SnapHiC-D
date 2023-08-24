import cooler 
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import cooler

chrom_size_path = "../../ext/hg19.chrom.sizes.ordered"

def visualize_hic(coolfile, coolRWRfile, outputfile): 
    c = cooler.Cooler(coolfile)
    c_rwr = cooler.Cooler(coolRWRfile)
    # Get the contact matrix
    chrom = "chr22"
    cool_matrix = c.matrix(balance=False).fetch(f"{chrom}")
    cool_rwr_matrix = c_rwr.matrix(balance=False).fetch(f"{chrom}")

    # Clear the plot
    plt.clf()

     # Create a subplot with 1 row and 2 columns
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    im1 = ax1.imshow(cool_matrix, cmap='coolwarm', origin="lower")
    ax1.set_title("Original scHi-C")
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Genomic Position")
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.set_label('Intensity')

    ax2 = axes[1]
    im2 = ax2.imshow(cool_rwr_matrix, cmap='coolwarm', origin="lower")
    ax2.set_title("RWR Imputed scHi-C")
    ax2.set_xlabel("Genomic Position")
    ax2.set_ylabel("Genomic Position")
    cbar2 = fig.colorbar(im2, ax=ax2)
    cbar2.set_label('Intensity')

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure
    plt.savefig(outputfile)

visualize_hic("/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/human_10089_GM12878_HFF-GM12878.R2.cool",
              "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/human_10089_GM12878_HFF-GM12878.R2_rwr.cool",
              "/home/maa160/SnapHiC-D/experiments/2023-08-24/test_visualizations/human_10089_GM12878_HFF-GM12878.R2_rwr.png")

visualize_hic("/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/human_9923_H1Esc_H1Esc.R1.cool",
               "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool_rwr/human_9923_H1Esc_H1Esc.R1_rwr.cool",
              "/home/maa160/SnapHiC-D/experiments/2023-08-24/test_visualizations//human_9923_H1Esc_H1Esc.R1_rwr.png")
