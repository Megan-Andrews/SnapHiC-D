import argparse
import cooler
import sys
import os
import numpy as np
from scipy.sparse import coo_matrix
import seaborn as sns
import pandas as pd
from tqdm import tqdm
import qnorm
from scipy import stats
import statsmodels.stats.multitest as multi
from mpi4py import MPI


def main():

    parser = create_parser()
    args = parser.parse_args()

    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()
    mpi_rank = mpi_comm.Get_rank()

    A_cools, B_cools = [], []
    for coolname in os.listdir(args.indir_A):
        A_cools.append(cooler.Cooler(os.path.join(args.indir_A, coolname)))
    for coolname in os.listdir(args.indir_B):
        B_cools.append(cooler.Cooler(os.path.join(args.indir_B, coolname)))

    chrom_tasks = [['chr{}'.format(i),'chr{}'.format(23-i)] for i in range(1,12)]
    if mpi_size > 10:
        my_chroms = chrom_tasks[mpi_rank]
        for chrom in my_chroms:
            print('rank {} is making candidates df for chromosome {}'.format(mpi_rank,chrom))
            candidates_outfile_path = os.path.join(args.output_dir, '{}_candidates.txt'.format(chrom))
            if not os.path.exists(candidates_outfile_path):
                candidates = make_candidates(A_cools, B_cools, chrom, args.resolution, args.snapHiC_dir, args.gen, min_distance = args.min_distance, max_distance = args.max_distance)
                candidates.to_csv(candidates_outfile_path, sep = "\t", index = False)
            else:
                print('candidates for chromosome {} exists.'.format(chrom))

    else:
        print('mpi size should be 11')

    mpi_comm.barrier()

    if mpi_rank == 0:
        groups = np.array(['A']*len(A_cools) + ['B']*len(B_cools))
        result_table = pd.DataFrame()
        for chr_num in np.arange(1,23):
            chrom = 'chr{}'.format(chr_num)
            candidates_path = os.path.join(args.output_dir, '{}_candidates.txt'.format(chrom))
            candidates = pd.read_csv(candidates_path, sep = "\t")
            gaps = np.unique((candidates['bin2_id']-candidates['bin1_id']).values)
            gaps = np.sort(gaps)
            for g in gaps:
                gap_candidates = candidates[candidates['bin2_id']-candidates['bin1_id']==g]
                gap_candidates.fillna(gap_candidates.min(), inplace=True)
                if result_table.empty:
                    result_table = get_result_table(gap_candidates, groups, 10000, 0.05, chrom)
                else:
                    result_table = pd.concat([result_table, get_result_table(gap_candidates, groups, 10000, 0.05, chrom)])
        result_path = os.path.join(args.output_dir, 'results.txt')
        result_table.to_csv(result_path, sep = "\t", index = False)


def upper_trim(x, trim):
    x = np.sort(x)
    trim_length = int(np.ceil(len(x) * (1 - trim)))
    trimmed_x = x[:trim_length]
    if len(trimmed_x) == 0:
        return 0,0
    trimmed_mean = np.mean(trimmed_x)
    trimmed_std = np.std(trimmed_x)
    return trimmed_mean, trimmed_std

def valid_bin_per_gap(num_bins, valid_bins, gap):
    cnt = 0
    for b in range(num_bins-gap):
        if (b in valid_bins) & ((b+gap) in valid_bins):
            cnt = cnt + 1
    return cnt

def get_diagonal_mean_sd(df, gap, num_bins, valid_bins):
    gap_valid_bins = valid_bin_per_gap(num_bins, valid_bins, gap)
    gap_values = np.zeros(gap_valid_bins)
    gap_nz_values = df.loc[df['gap']==gap,'count']
    gap_values[:gap_nz_values.shape[0]] = gap_nz_values
    mean, std = upper_trim(gap_values, 0.01)
    return mean, std

def diagonal_norm(cool, chrom, resolution, min_distance = 20000, max_distance = 2000000):
    df = cool.matrix(balance=False, as_pixels=True).fetch(chrom)
    df[['bin1_id','bin2_id']] = df[['bin1_id','bin2_id']] - cool.offset(chrom)
    df['gap'] = df['bin2_id'] - df['bin1_id']
    valid_bins = np.unique(list(df['bin1_id']) + list(df['bin2_id']))
    chr_size = cool.chromsizes[chrom]
    num_bins = np.ceil(chr_size/resolution).astype(int)
    min_gap, max_gap = int(min_distance/resolution), int(max_distance/resolution)
    gap_means, gap_stds = np.zeros(max_gap+1), np.zeros(max_gap+1)
    for g in range(min_gap,max_gap+1):
        gap_means[g], gap_stds[g] = get_diagonal_mean_sd(df, g, num_bins, valid_bins)
    std_th = 0.0001
    df = df[(df['gap']>=min_gap) & (df['gap']<=max_gap)]
    df['norm_count'] = [(count-gap_means[d])/gap_stds[d] if gap_stds[d]>std_th else count for count, d in zip(df['count'],df['gap'])]
    return df

def get_ID(CHR,BINSIZE,SnapHiC_dir,genome):
    if(genome == "hg19"):
        b = pd.read_csv(os.path.join(SnapHiC_dir, "ext/hg19_filter_regions.txt"), sep="\t", header=None)
        g0 = pd.read_csv(os.path.join(SnapHiC_dir, "ext/hg19.refGene.transcript.TSS.061421.txt"), sep="\t", header='infer')
    if(genome == "mm10"):
        b = pd.read_csv(os.path.join(SnapHiC_dir, "ext/mm10_filter_regions.txt"), sep="\t", header=None)
        g0 = pd.read_csv(os.path.join(SnapHiC_dir, "ext/mm10.refGene.transcript.TSS.061421.txt"), sep="\t", header='infer')

    b.columns = ['chr_name', 'x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6']
    b1 = b.loc[b['chr_name'] == CHR]
    # is a bin containing unmappable regions but not all unmappable valid?
    b_ID = ((b1.iloc[:, 2].values)/BINSIZE).astype(int)


    g = g0[g0['chr'] == CHR]
    g['TSS.Bin'] = np.ceil(g['TSS']/BINSIZE)
    gID = np.unique(g['TSS.Bin'])
    return b_ID, gID

def make_candidates(A_cools, B_cools, chrom, resolution, snapHiC_dir, genome, min_distance = 20000, max_distance = 2000000):
    b_id, g_id = get_ID(chrom, resolution, snapHiC_dir, genome)
    all_df = pd.DataFrame()
    for c in tqdm(A_cools+B_cools):
        norm_df = diagonal_norm(c, chrom, resolution)
        norm_df = norm_df[(~norm_df['bin1_id'].isin(b_id)) & (~norm_df['bin2_id'].isin(b_id))]
        norm_df = norm_df[(norm_df['bin1_id'].isin(g_id)) | (norm_df['bin2_id'].isin(g_id))]
        if all_df.empty:
            all_df = norm_df[['bin1_id','bin2_id','norm_count']]
        else:
            all_df = pd.merge(all_df, norm_df[['bin1_id','bin2_id','norm_count']], on = ['bin1_id','bin2_id'], how='outer')
    all_df.columns = ['bin1_id','bin2_id'] + ['cell{}'.format(i) for i in np.arange(1,all_df.shape[1]-1)]
    return all_df
# todo: add value filtering to make_candidates
def get_result_table(candidates, groups, resolution, fdr_t, chrom):
    A_scores = candidates.iloc[:,2:].iloc[:, np.where(groups=='A')[0]]
    B_scores = candidates.iloc[:,2:].iloc[:, np.where(groups=='B')[0]]
    num_A, num_B = np.where(groups=='A')[0].shape[0], np.where(groups=='B')[0].shape[0]
    sigA = np.count_nonzero(A_scores > 1.96, axis=1) > num_A*.1
    sigB = np.count_nonzero(B_scores > 1.96, axis=1) > num_B*.1
    sig = np.array([a or b for a, b in zip(sigA,sigB)]).astype(int)
    candidates = pd.concat([candidates.iloc[:,:2], A_scores, B_scores], axis = 1)[sig==1]
    candidates_scores = candidates.iloc[:,2:]
    candidates_scores = qnorm.quantile_normalize(candidates_scores, axis = 1)
    A_scores = candidates_scores.iloc[:, np.where(groups=='A')[0]]
    B_scores = candidates_scores.iloc[:, np.where(groups=='B')[0]]
    if all(np.var(A_scores, ddof=1, axis=0) > 0) is True and all(np.var(B_scores, ddof=1, axis=0) > 0) is True:
        t_stat, p_val = stats.ttest_ind(A_scores, B_scores, axis=1, equal_var=False)
        out = pd.DataFrame(columns=['chrom', 'bin1_id','bin2_id', 'mean.A','mean.B', 'Tstat', 'Ttest.Pvalue'])
        out['chrom'] = [chrom]* candidates_scores.shape[0]
        out[['bin1_id','bin2_id']] = np.array(candidates.iloc[:,:2]*resolution)
        out['mean.A'] = np.array(np.mean(A_scores, axis=1)) # should we save mean of quantile-normalized scores?
        out['mean.B'] = np.array(np.mean(B_scores, axis=1))
        out['Tstat'] = t_stat
        out['Ttest.Pvalue'] = p_val
        rej, pval_corr = multi.multipletests(out['Ttest.Pvalue'], method='fdr_bh')[:2]
        out['fdr'] = pval_corr
        out['significant'] = (out['fdr'] < fdr_t) & (abs(out['Tstat']) > 2)
    else:
        out = pd.DataFrame(columns=['chrom','bin1_id','bin2_id', 'mean.A','mean.B', 'Tstat', 'Ttest.Pvalue','fdr','significant'])
    return out


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir-A', action = 'store', required = True, \
                        help = 'path of the folder including imputed cools for the first cell type')
    parser.add_argument('--indir-B', action = 'store', required = True, \
                        help = 'path of the folder including imputed cools for the second cell type')
    parser.add_argument('-o', '--output-dir', action = 'store', required = True, help = 'output directory to save candidate data frames for all chromosomes and results')
    parser.add_argument('-r', '--resolution', action = 'store', required = True, type = int, help = "resolution of contact maps")
    parser.add_argument('--min-distance', action = 'store', required = True, type = int, help = "minimum distance between loop anchors")
    parser.add_argument('--max-distance', action = 'store', required = True, type = int, help = "maximum distance between loop anchors")
    parser.add_argument('--snapHiC-dir', action = 'store', required = True, help = 'path of snapHiC directory')
    parser.add_argument('-g', '--gen', action = 'store', help = 'genome - mm10 or hg19', \
                        required = True)
    return parser


if __name__ == "__main__":
    main()
