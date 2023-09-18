import os

# Directory where the files are located
directory = "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool"  # You can change this to the directory containing your files

# Get a list of all files in the directory
files = os.listdir(directory)

# Initialize lists for GM12878 and H1Esc files
gm12878_files = []
h1esc_files = []

# Iterate through the files and categorize them
for filename in files:
    if filename.endswith(".cool"):
        parts = filename.split("_")
        
        # Check if the third part of the filename contains 'GM12878' or 'H1Esc'
        if len(parts) > 2 and 'GM12878' in parts[2]:
            gm12878_files.append(filename)
        elif len(parts) > 2 and 'H1Esc' in parts[2]:
            h1esc_files.append(filename)

# Save the lists to separate text files
with open("~/SnapHiC-D/experiments/2023-09-18/GM12878_files.txt", "w") as gm12878_file:
    gm12878_file.write("\n".join(gm12878_files))

with open("~/SnapHiC-D/experiments/2023-09-18/H1Esc_files.txt", "w") as h1esc_file:
    h1esc_file.write("\n".join(h1esc_files))

print("Files categorization and lists creation completed.")