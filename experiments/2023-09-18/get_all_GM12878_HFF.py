import os

# Directory where the files are located
directory = "/project/compbio-lab/scHi-C/Kim2020/kim2020_cool_rwr3"  # You can change this to the directory containing your files


# Get a list of all files in the directory
files = os.listdir(directory)

# Initialize lists for GM12878 and HFF files
gm12878_files = []
hff_files = []

# Iterate through the files and categorize them
for filename in files:
    if filename.endswith(".cool"):
        parts = filename.split("_")
        
        # Check if the third part of the filename contains 'GM12878' or 'HFF'
        if len(parts) > 2 and 'GM12878' in parts[2] and 'HFF-GM12878.R1' in parts[3]:
            gm12878_files.append(filename.replace("_rwr.cool", ".cool"))
        elif len(parts) > 2 and 'HFF' in parts[2] and 'HFF-GM12878.R1' in parts[3]:
            hff_files.append(filename.replace("_rwr.cool", ".cool"))

# Save the lists to separate text files
with open("GM12878_files_2.txt", "w") as gm12878_file:
    gm12878_file.write("\n".join(gm12878_files))

with open("HFF_files_2.txt", "w") as HFF_file:
    HFF_file.write("\n".join(hff_files))

print("Files categorization and lists creation completed.")
