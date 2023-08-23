import cooler 
import pandas as pd
import numpy as np
import os

# read .matrix file
libraries_list = ["H1Esc-HFF.R1",  "H1Esc.R1",  "HFF-GM12878.R1",  "IMR90-HAP1.R1", "GM12878_IMR90.R1",  "H1Esc-HFF.R2",  "H1Esc.R2",  "HFF-GM12878.R2",  "IMR90-HAP1.R2"]
kim2020_file_path = "/project/compbio-lab/scHi-C/Kim2020/"
chrom_size_path = "../../ext/hg19.chrom.sizes"
cool_directory = os.path.join(kim2020_file_path, "Kim2020_cool")
#os.mkdir(cool_directory)

def create_cooler(df, output_file):
    def get_chrom_offsets(bins_df):
        chrom_offset = {chrom: bins_df[bins_df['chrom'] == chrom].index[0]
                        for chrom in bins_df['chrom'].cat.categories}
        return chrom_offset

    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t',index_col=0, header=None).squeeze(axis=1)
    chrom_sizes = chrom_sizes[:24]
    bins_df = cooler.binnify(chrom_sizes, 500000)
    chrom_offsets = get_chrom_offsets(bins_df)

    for chr in chrom_sizes.index:
        # Define a boolean mask for the condition
        mask = df["1st_chr"] == chr
        
        # Update binId_1 and binId_2 using .loc
        df.loc[mask, "binId_1"] += chrom_offsets[chr]
        df.loc[mask, "binId_2"] += chrom_offsets[chr]

    data = pd.DataFrame({
            "bin1_id": df["binId_1"],
            "bin2_id": df["binId_2"],
            "count": df["counts"]
        }
    )

    cooler.create_cooler(output_file, 
                            bins=bins_df, 
                            pixels=data,
                            ordered=False,
                            dtypes={'count': np.float64}) 


def filter_matrices(file_path, output_file):
    columns = ["binId_1", "binId_2", "counts","normalized_counts", "1st_chr", "2nd_chr"]
    df = pd.read_csv(file_path, delimiter='\t', header=None)
    df.columns = columns
    df = df[df["1st_chr"] == df["2nd_chr"]] # only include intra-chromosome pairs
    df["1st_chr"] = df["1st_chr"].apply(lambda x: x.split("_")[1])
    temp_df = df[df["binId_1"] != df["binId_2"]] # only count interactions more than 500kb apart

    if temp_df["counts"].sum() > 2000:
        print("more than 2000 read pairs")
        create_cooler(df, output_file)
    else:
        print("less than 2000 read paris")
        

# Iterate through the list of directories
for library in libraries_list:
    directory = os.path.join(kim2020_file_path, library) 
    
    # Check if the item in the list is a directory
    if os.path.isdir(directory):
        # Perform an action, such as printing the directory name
        print(f"Directory: {directory}")
        print(library)
        labels_directory = kim2020_file_path + "labels/" + library + ".labeled"
        cell_type_df = pd.read_csv(labels_directory, header=None, delimiter="\t")
        cell_type_df.columns = ["file_name", "cell_type"]
        
        # Iterate through all files in the directory
        for filename in os.listdir(directory):
            file_no = filename.split("_")[1]
            cell_line = cell_type_df[cell_type_df["file_name"] == filename]["cell_type"]
            print(cell_line)
            cool_output_file = os.path.join(cool_directory, f"human_{file_no}_{cell_line}_{library}.cool")
            if os.path.isfile(os.path.join(directory, filename)):
                # Perform an action on each file
                print(f"File: {filename}")
                print(f"Output file: {cool_output_file}")
                filter_matrices(os.path.join(directory, filename), cool_output_file)
    else:
        print(f"Not a directory: {directory}")


