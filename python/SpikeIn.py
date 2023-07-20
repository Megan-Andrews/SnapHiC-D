import argparse
import cooler
import os
import pandas as pd
from tqdm import tqdm
from random import sample
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from numpy.random import normal

CHR = 'chr22'
th = 1.5

def main():
    parser = create_parser()
    args = parser.parse_args()

    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    cond1_coollist = list(pd.read_csv(args.cond1_coollist, header = None)[0])
    cond2_coollist = list(pd.read_csv(args.cond2_coollist, header = None)[0])
    keyCoord = ["bin1_id", "bin2_id"]
    keyC1 = ["cond1_" + str(i) for i in range(1, len(cond1_coollist) + 1)]
    keyC2 = ["cond2_" + str(i) for i in range(1, len(cond2_coollist) + 1)]

    total_df1 = None
    ind = 1
    print('reading condition 1 files.')
    for cp in tqdm(cond1_coollist):
        if total_df1 is None:
            total_df1 = read_cool(cp, 'cond1_{}'.format(ind))
        else:
            df = read_cool(cp, 'cond1_{}'.format(ind))
            total_df1 = pd.merge(total_df1, df, how = "outer", on = ["bin1_id", "bin2_id"]).fillna(0)
        ind += 1

    total_df2 = None
    ind = 1
    print('reading condition 2 files.')
    for cp in tqdm(cond2_coollist):
        if total_df2 is None:
            total_df2 = read_cool(cp, 'cond2_{}'.format(ind))
        else:
            df = read_cool(cp, 'cond2_{}'.format(ind))
            total_df2 = pd.merge(total_df2, df, how = "outer", on = ["bin1_id", "bin2_id"]).fillna(0)
        ind += 1
    total_df1 = pd.merge(total_df1, total_df2, how = "outer", on = ["bin1_id", "bin2_id"]).fillna(0)
    del(total_df2)

    print('learning the dependece between mean and variance in condition1 (intrinsic variability)')

    #non_diag_cond = total_df1['bin1_id'] != total_df1['bin2_id']
    means, vars = total_df1[keyC2].mean(axis=1), total_df1[keyC2].var(axis=1)
    mean_stdev_df = pd.DataFrame({'mean': means, 'stdev': np.sqrt(vars)})
    mean_stdev_df.sort_values('mean', inplace=True)
    mean_stdev_df.drop_duplicates('mean', inplace=True)
    us = UnivariateSpline(mean_stdev_df['mean'], mean_stdev_df['stdev'], s = 1000)
    mean_stdev_df['us_pred'] = us(mean_stdev_df['mean'])

    plt.scatter(mean_stdev_df['mean'], mean_stdev_df['stdev'])
    plt.plot(mean_stdev_df['mean'], mean_stdev_df['us_pred'], color='Red')
    plt.savefig(os.path.join(output_dir, 'cond1_mean_stdev.png'))

    print('finding DCC candidates')

    avg_C1 = total_df1[keyC1].mean(axis=1)
    avg_C2 = total_df1[keyC2].mean(axis=1)
    logFCs = np.log(avg_C2/avg_C1)
    #total_df1['logFCs'] = logFCs
    candidates_idx = (np.isfinite(logFCs)) & (abs(logFCs)>=th) & ((avg_C1 > 0.1) | (avg_C2 > 0.1))
    #total_df1['is_DCC'] = candidates_idx
    total_df1['FC'] = [np.exp(logFC) if is_DCC else 1 for logFC, is_DCC in zip(logFCs, candidates_idx)]
    DCC_candidates = pd.concat([total_df1[keyCoord][candidates_idx], logFCs[candidates_idx]], axis = 1)
    DCC_candidates.to_csv(os.path.join(output_dir, 'DCC.txt'), sep = "\t", header = None, index = False)

    print('writing simulated data.')

    cond1_dir = os.path.join(output_dir, 'cond1')
    cond2_dir = os.path.join(output_dir, 'cond2')
    if not os.path.exists(cond1_dir):
        os.mkdir(cond1_dir)
    if not os.path.exists(cond2_dir):
        os.mkdir(cond2_dir)
    for k, key in tqdm(enumerate(keyC1)):
        df = total_df1[keyCoord+[key, 'FC']]
        df = df[df[key]!=0]
        #diag_cond = df['bin1_id'] == df['bin2_id']
        #intrinsic_noise = [0 if is_diag else normal(0, np.sqrt(us(m))) for m, is_diag in zip(df[key], diag_cond)]
        #print(us(df[key]))
        intrinsic_noise = [normal(0, us(m)) for m in df[key]]
        df['cond2'] = np.maximum(np.round((df[key] + intrinsic_noise)*df['FC']), 0)
        df[keyCoord + [key]].to_csv(os.path.join(cond1_dir, 'cell{}.txt'.format(k)), sep = "\t", header = None, index = False)
        df = df[df['cond2']!=0]
        df[keyCoord + ['cond2']].to_csv(os.path.join(cond2_dir, 'cell{}.txt'.format(k)), sep = "\t", header = None, index = False)

def read_cool(coolpath, count_colname):

    c = cooler.Cooler(coolpath)
    df = c.matrix(balance = False, as_pixels = True).fetch(CHR)
    df.loc[:,['bin1_id','bin2_id']] = df.loc[:,['bin1_id','bin2_id']] - c.offset(CHR)
    #df = df[df['bin2_id'] < df['bin1_id'] + 200]
    df.rename(columns = {'count': count_colname}, inplace = True)
    return df

def create_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cond1-coollist', action = 'store', required = True, help = 'path of a file including coolpaths from condition 1')
    parser.add_argument('--cond2-coollist', action = 'store', required = True, help = 'path of a file including coolpaths from condition 2')
    parser.add_argument('-o', '--output-dir', action = 'store', required = True, help = 'output directory to save simulated spike-ined data')

    return parser


if __name__ == "__main__":
    main()
