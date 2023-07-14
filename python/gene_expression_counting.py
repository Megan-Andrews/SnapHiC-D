import pandas as pd

# read csv file
df = pd.read_csv("/home/maa160/SnapHiC-D/experiments/2023-07-10/diffHiC_2x2/diffHiC_MG_Astro_results.csv", sep=",")

# count all the interactions
print("All interactions: ", df.count())

# count the values of logFC that are positive (MG specific interactions)
print("MG specific interactions: ", df[df["logFC"] > 0].count())

#count the values of logFC that are negative (Astro specific interactions)
print("Astro specific interactions: ", df[df["logFC"] < 0].count())