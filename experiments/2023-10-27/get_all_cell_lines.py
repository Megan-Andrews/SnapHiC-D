import os

# Directory where the files are located

directory = "/project/compbio-lab/scHi-C/Kim2020/kim2020_cool_rwr3"  # You can change this to the directory containing your files

# Get a list of all files in the directory
files = os.listdir(directory)

# Initialize lists for GM12878 and HFF files
gm12878_files = []
hff_files = []
imr90_files = []
k562_files = []
kbm7_files = []
h1hesc_files = []

#[ "GM12878" , "H1-hESC",  "HFF",  "IMR90" , "K562" , "KBM7" ]

# Iterate through the files and categorize them
for filename in files:
    if filename.endswith(".cool"):
        parts = filename.split("_")
        
        # Check if the third part of the filename contains 'GM12878' or 'HFF'
        if len(parts) > 2 and 'GM12878' in parts[2]:# and 'HFF-GM12878.R1' in parts[3]:
            gm12878_files.append(filename.replace("_rwr.cool", ".cool"))
        elif len(parts) > 2 and 'HFF' in parts[2]:# and 'HFF-GM12878.R1' in parts[3]:
            hff_files.append(filename.replace("_rwr.cool", ".cool"))
        elif len(parts) > 2 and 'H1Esc' in parts[2]:# and 'H1Esc.R2' in parts[3]:
            h1hesc_files.append(filename.replace("_rwr.cool", ".cool"))
        elif len(parts) > 2 and 'IMR90' in parts[2]:# and 'IMR90-HAP1.R2' in parts[3]:
            imr90_files.append(filename.replace("_rwr.cool", ".cool"))

# Save the lists to separate text files
with open("GM12878_files.txt", "w") as gm12878_file:
    gm12878_file.write("\n".join(gm12878_files))

with open("HFF_files.txt", "w") as HFF_file:
    HFF_file.write("\n".join(hff_files))

with open("IMR90_files.txt", "w") as IMR90_file:
    IMR90_file.write("\n".join(imr90_files))

with open("H1hESC_files.txt", "w") as H1hESC_file:
    H1hESC_file.write("\n".join(h1hesc_files))
print("Number of each type:")
print("GM12878: ", len(gm12878_files))
print("HFF: ", len(hff_files))
print("IMR90: ", len(imr90_files))
print("H1hESC: ", len(h1hesc_files))
print("Files categorization and lists creation completed.")
