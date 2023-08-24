import os
import pandas as pd
import numpy as np

file_array = []
for file in os.listdir("/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool"):
    file_array.append(file)

file_df = pd.DataFrame(file_array, columns=["file_name"])
file_df["cell_line"] = file_df["file_name"].apply(lambda x: x.split("_")[0])
file_df["library"] = file_df["file_name"].apply(lambda x: x.split("_")[1])

print(file_df["cell_line"].value_counts())
print(file_df["library"].value_counts())

