import pandas as pd
import math
import numpy as np
import os
import matplotlib.pyplot as plt



astro_MG = pd.read_excel('NIHMS739375-supplement-3.xlsx', skiprows=1, index_col=None)
astro_MG.columns = ['gene','B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z','AA', 'AB', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AJ', 'AK', 'AL', 'AM', 'AN', 'AO', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AU', 'AV', 'AW', 'AX', 'AY', 'AZ', 'BA', 'BB', 'BC', 'BD', 'BE', 'BF', 'BG', 'BH', 'BI', 'BJ', 'BK']
astro_MG['gene'] = astro_MG['gene'].str.upper().str.strip()
astro_columns = ['P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA']
MG_columns = ['AH', 'AI', 'AJ']
astro_MG['avg_astro'] = astro_MG[astro_columns].mean(axis=1)
astro_MG['avg_MG'] = astro_MG[MG_columns].mean(axis=1)
astro_MG['astro/MG'] = astro_MG['avg_astro']/astro_MG['avg_MG'] 
astro_MG['log(astro/MG)'] = astro_MG['astro/MG'].apply(lambda x: math.log(x))
astro_MG['|log(astro/MG)|'] = astro_MG['log(astro/MG)'].abs()
astro_MG = astro_MG[['gene', 'avg_astro','avg_MG', 'astro/MG', 'log(astro/MG)', '|log(astro/MG)|']]

gene_info = pd.read_csv('gene_info.txt', sep='\t')
gene_info.columns = ['0', 'chrom', 'start', 'end', '4', '5', 'gene', '7' ]
gene_info = gene_info.iloc[:,[1,2,3,6]]
gene_info['gene'] = gene_info['gene'].str.upper().str.strip()

result = pd.merge(astro_MG, gene_info,on='gene', how='left')
# result['floored_start'] = result['start'].fillna(-100000).apply(lambda x: math.floor(x / 100000) * 100000)
# result['floored_end'] = result['end'].fillna(-200000).apply(lambda x: math.ciel(x / 100000) * 100000)
# result['gap'] = (result['floored_end']-result['floored_start'])/100000
#result['gap/resolution'] = result['gap'].fillna(-100000).floordiv(100000).astype(int)

# result.to_csv("supplement_with_gene_info.csv")


# Assume the directory path where the files are located
directory = '/project/compbio-lab/scHi-C/Lee2019/results/2023-06-06/Astro_MG_diff/'


def within_range(row, start, end):
    # print(row)
    bin1_id = int(row['bin1_id'])
    bin2_id = int(row['bin2_id'])
    if start >= bin1_id and start <= bin1_id+100000 and end >= bin1_id and end <= bin1_id + 100000:
        return True
    elif start >= bin2_id and start <= bin2_id+100000 and end >= bin2_id and end <= bin2_id + 100000:
        return True
    else: return False

list_fold_change = list()
list_Tstat = list()

for index, row in result.iterrows():
    start = row['start']
    end = row['end']
    fold_change = row['astro/MG']


    # Iterate over the files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        diff = pd.read_csv(file_path, sep='\t', names=['bin1_id','bin2_id','mean.A','mean.B','Tstat','Ttest.Pvalue','fdr','significant'], skiprows=1)
        diff = diff[diff.apply(within_range, axis=1, args=(start, end))]
        # for diffIndex, diffRow in diff.iterrows():
        list_fold_change.extend([fold_change] * len(diff))
        list_Tstat.extend(diff['Tstat'])
        if len(diff) != 0:
            print(diff, start, end)
    print(index)


# Plot the scatterplot
plt.scatter(list_fold_change, list_Tstat, label='Fold Change Vs. Tstat')

# Calculate the trendline using numpy
slope, intercept = np.polyfit(list_fold_change, list_Tstat, deg=1)
trendline = slope * np.array(list_fold_change) + intercept

# Plot the trendline
plt.plot(list_fold_change, trendline, color='red', label='Trendline')

# Add labels and legend
plt.xlabel('Fold Change')
plt.ylabel('Tstat')
plt.legend()

# Show the plot
plt.show()
# result.to_csv("supplement_with_gene_info.csv")

