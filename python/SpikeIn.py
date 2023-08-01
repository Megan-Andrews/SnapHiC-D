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
th = 2


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
    keyNoise1 = ["sim_cond1_" + str(i) for i in range(1, len(cond1_coollist) + 1)]
    keyNoise2 = ["sim_cond2_" + str(i) for i in range(1, len(cond1_coollist) + 1)]

    total_df = read_conditions(cond1_coollist, cond2_coollist)

    print('learning the dependece between mean and variance in condition1 (intrinsic variability)')
    cond1_mean_stdevs = get_mean_stdev(total_df, keyC1)
    diag_us = UnivariateSpline(cond1_mean_stdevs['diag']['mean'], cond1_mean_stdevs['diag']['stdev'], s = 100000)
    nondiag_us = UnivariateSpline(cond1_mean_stdevs['nondiag']['mean'], cond1_mean_stdevs['nondiag']['stdev'], s = 100000)
    cond1_mean_stdevs['diag']['us_pred'] = diag_us(cond1_mean_stdevs['diag']['mean'])
    cond1_mean_stdevs['nondiag']['us_pred'] = nondiag_us(cond1_mean_stdevs['nondiag']['mean'])

    plt.scatter(cond1_mean_stdevs['diag']['mean'], cond1_mean_stdevs['diag']['stdev'])
    plt.plot(cond1_mean_stdevs['diag']['mean'], cond1_mean_stdevs['diag']['us_pred'], color='Red', label = 'diagonal IFs')
    plt.scatter(cond1_mean_stdevs['nondiag']['mean'], cond1_mean_stdevs['nondiag']['stdev'])
    plt.plot(cond1_mean_stdevs['nondiag']['mean'], cond1_mean_stdevs['nondiag']['us_pred'], color='Green', label = 'nondiagonal IFs')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'cond1_mean_stdev.png'))

    print('finding DCC candidates.')
    FCs, DCC_df = get_DCC_candidates(total_df, keyCoord, keyC1, keyC2, th)
    DCC_df.to_csv(os.path.join(output_dir, 'DCC.txt'), sep = "\t", header = None, index = False)

    print('injecting noise and DCCs')
    noisy_total_df = inject_noise_and_DCC(total_df, keyCoord, keyC1, keyNoise1, keyNoise2, diag_us, nondiag_us, FCs)

    print('writing simulated data.')
    bins = cooler.Cooler(cond1_coollist[0]).bins().fetch(CHR)
    bins.index = range(bins.shape[0])

    subdirs = [os.path.join(output_dir, 'cond1'), os.path.join(output_dir, 'cond2'), os.path.join(output_dir, 'sim_cond1'), os.path.join(output_dir, 'sim_cond2')]

    for subdir in subdirs:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

    column_keys = [keyC1, keyC2, keyNoise1, keyNoise2]

    for k in range(4):
        c = 1
        for col_key in column_keys[k]:
            pixels = noisy_total_df[keyCoord + [col_key]].copy()
            pixels.fillna(0, inplace = True)
            pixels = pixels[pixels[col_key] != 0]
            pixels.rename(columns = {col_key: 'count'}, inplace = True)
            cooler.create_cooler(os.path.join(subdirs[k], 'cell{}.cool'.format(c)),
                        bins = bins,
                        pixels = pixels)
            c += 1

    for cond, subdir in zip(['cond1', 'cond2', 'sim_cond1', 'sim_cond2'], subdirs):
        coollist = pd.DataFrame({'cp': [os.path.join(subdir, cp) for cp in os.listdir(subdir)]})
        coollist_path = os.path.join(output_dir, '{}_coollist.txt'.format(cond))
        coollist.to_csv(coollist_path, index = False, header = None)

def read_conditions(cond1_coollist, cond2_coollist):

    total_df = None
    ind = 1
    print('reading condition 1 files.')
    for cp in tqdm(cond1_coollist):
        if total_df is None:
            total_df = read_cool(cp, 'cond1_{}'.format(ind))
        else:
            df = read_cool(cp, 'cond1_{}'.format(ind))
            total_df = pd.merge(total_df, df, how = "outer", on = ["bin1_id", "bin2_id"]).fillna(0)
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
    total_df = pd.merge(total_df, total_df2, how = "outer", on = ["bin1_id", "bin2_id"]).fillna(0)
    return total_df

def read_cool(coolpath, count_colname):

    c = cooler.Cooler(coolpath)
    df = c.matrix(balance = False, as_pixels = True).fetch(CHR)
    df.loc[:,['bin1_id','bin2_id']] = df.loc[:,['bin1_id','bin2_id']] - c.offset(CHR)
    #df = df[df['bin2_id'] < df['bin1_id'] + 200]
    df.rename(columns = {'count': count_colname}, inplace = True)
    return df

def get_mean_stdev(total_df, key):
    valid_IFs = (total_df[key]!=0).sum(axis = 1) >= 2
    diag_cond = (total_df['bin1_id'] == total_df['bin2_id'])
    diag_means, diag_stdevs = total_df[valid_IFs & diag_cond][key].mean(axis=1), total_df[valid_IFs & diag_cond][key].std(axis=1)
    nondiag_means, nondiag_stdevs = total_df[valid_IFs & ~diag_cond][key].mean(axis=1), total_df[valid_IFs & ~diag_cond][key].std(axis=1)
    def make_mean_stdev_df(means, stdevs):
        mean_stdev = pd.DataFrame({'mean': means, 'stdev': stdevs})
        mean_stdev.sort_values('mean', inplace=True)
        mean_stdev.drop_duplicates('mean', inplace=True)
        return mean_stdev
    mean_stdevs = {'diag': make_mean_stdev_df(diag_means, diag_stdevs),
                  'nondiag': make_mean_stdev_df(nondiag_means, nondiag_stdevs)}
    return mean_stdevs

def get_DCC_candidates(total_df, keyCoord, keyC1, keyC2, th):
    seqDep = total_df[keyC1 + keyC2].sum(axis = 0)
    norm_total_df = pd.concat([total_df[keyCoord], total_df[keyC1 + keyC2]/seqDep * seqDep.max()], axis = 1)
    avg_C1 = norm_total_df[keyC1].mean(axis=1)
    avg_C2 = norm_total_df[keyC2].mean(axis=1)
    logFCs = np.log2(avg_C2/avg_C1)
    candidates_idx = (np.isfinite(logFCs)) & (abs(logFCs)>=th) & ((avg_C1 > 0.1) | (avg_C2 > 0.1))
    FCs = np.power(2, logFCs)
    FCs[~candidates_idx] = 1
    DCC_candidates = pd.concat([total_df[keyCoord][candidates_idx], logFCs[candidates_idx]], axis = 1)
    return FCs, DCC_candidates

def inject_noise_and_DCC(df, coordKeys, orig_keys, noisy_keys1, noisy_keys2, diag_us, nondiag_us, FCs):
    diag_idx = df[coordKeys[0]] == df[coordKeys[1]]
    orig_counts = df[orig_keys].copy()
    noisy_counts = pd.DataFrame(0, index = df.index, columns = noisy_keys1)
    noisy_counts[diag_idx] = np.maximum(0, np.round((orig_counts[diag_idx] +
                                                     normal(0,diag_us(orig_counts[diag_idx])))))

    noisy_counts[~diag_idx] = np.maximum(0, np.round((orig_counts[~diag_idx] +
                                                      normal(0,nondiag_us(orig_counts[~diag_idx])))))

    noisy_DCC_counts = pd.DataFrame(0, index = df.index, columns = noisy_keys2)
    noisy_DCC_counts[diag_idx] = np.maximum(0, np.round((orig_counts[diag_idx] +
                                                     normal(0,diag_us(orig_counts[diag_idx]))).mul(FCs, axis = 0)))

    noisy_DCC_counts[~diag_idx] = np.maximum(0, np.round((orig_counts[~diag_idx] +
                                                      normal(0,nondiag_us(orig_counts[~diag_idx]))).mul(FCs, axis = 0)))

    df = pd.concat([df, noisy_counts, noisy_DCC_counts], axis = 1)
    return df




def create_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cond1-coollist', action = 'store', required = True, help = 'path of a file including coolpaths from condition 1')
    parser.add_argument('--cond2-coollist', action = 'store', required = True, help = 'path of a file including coolpaths from condition 2')
    parser.add_argument('-o', '--output-dir', action = 'store', required = True, help = 'output directory to save simulated spike-ined data')

    return parser


if __name__ == "__main__":
    main()
