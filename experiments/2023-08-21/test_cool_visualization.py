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

    # Clear the plot
    plt.clf()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    ax1.imshow(cool_matrix, cmap='coolwarm', origin="lower")
    ax1.set_title("Cool Matrix")
    ax1.set_xlabel("Genomic Position")
    ax1.set_ylabel("Genomic Position")
    cbar1 = ax1.colorbar()
    cbar1.set_label('Intensity')

    # Plot the matrix_df on the second subplot
    ax2 = axes[1]
    ax2.imshow(matrix_df['counts'].values.reshape(cool_matrix.shape), cmap='coolwarm', origin="lower")
    ax2.set_title("Matrix DataFrame")
    ax2.set_xlabel("Genomic Position")
    ax2.set_ylabel("Genomic Position")
    cbar2 = ax2.colorbar()
    cbar2.set_label('Intensity')

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure
    plt.savefig(outputfile)


visualize_hic("/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/181218_21yr_2_A1_AD004_L23_100kb_contacts_imputed.cool","/home/maa160/CompBioRA2023/Hi-C_Visualizations/181218_21yr_2_A1_AD004_L23_100kb_contacts_imputed.png")
