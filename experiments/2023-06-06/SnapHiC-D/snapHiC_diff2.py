import numpy as np
import pandas as pd
import qnorm
from scipy import stats
import statsmodels.stats.multitest as multi
import concurrent.futures
from functools import partial
import argparse
import cooler
import sys
import os



def main():
    parser = create_parser()
    args = parser.parse_args()

    CHR =  args.chrom
    genome = args.gen
    SnapHiC_dir = args.indirSHD
    dir_in_A = args.indirA
    dir_in_B = args.indirB
    out_dir = args.outdir

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    min_gap = args.mini_gap
    max_gap = args.maxi_gap

    BINSIZE = args.binsize

    fdr_t = args.fdr_threshold

    max_worker = args.num_cpus

    A_coolnames = os.listdir(dir_in_A)
    B_coolnames = os.listdir(dir_in_B)
    A_coolpaths = [os.path.join(dir_in_A, A_coolname) for A_coolname in A_coolnames]
    B_coolpaths = [os.path.join(dir_in_B, B_coolname) for B_coolname in B_coolnames]
    all_pixels = pd.DataFrame()

    print('Start reading pixels...')

    for coolname,coolpath in zip(A_coolnames+B_coolnames,A_coolpaths+B_coolpaths):
        c = cooler.Cooler(coolpath)
        pixels = c.matrix(balance=False, as_pixels=True).fetch(CHR)
        pixels.columns = ['bin1_id', 'bin2_id', coolname.split('.')[0]]
        if all_pixels.empty:
            all_pixels = pixels
        else:
            all_pixels = pd.merge(all_pixels, pixels, on = ['bin1_id','bin2_id'], how = 'outer')

    print('Finished reading pixels...')
    global num_A
    global num_B
    global cell_types
    num_A, num_B = len(A_coolnames), len(B_coolnames)
    cell_types = np.array(['A']*num_A + ['B']*num_B)

    b_ID, gID = get_ID(CHR,BINSIZE,SnapHiC_dir,genome)
    GAP=list(range(min_gap, max_gap))
    with concurrent.futures.ProcessPoolExecutor(max_worker) as executor:
        for result in executor.map(partial(gap_diff_analysis, pixels=all_pixels, out_dir = out_dir, b_ID = b_ID, gID = gID, fdr_t = fdr_t, binsize=BINSIZE) , GAP):

            print("..")

    print("Exiting program")

def g_filter(x):
    counts = np.count_nonzero(x > 1.96, axis=1)
    return counts



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
    b_ID = b1.iloc[:, 2]/BINSIZE


    g = g0[g0['chr'] == CHR]
    g['TSS.Bin'] = np.ceil(g['TSS']/BINSIZE)
    gID = np.unique(g['TSS.Bin'])
    return b_ID, gID


def gap_diff_analysis(GAP, pixels, out_dir, b_ID, gID, fdr_t, binsize):

    print("Starting with GAP=", GAP)

    candidates = pixels[pixels['bin2_id']-pixels['bin1_id']==GAP]
    candidates.fillna(candidates.min(), inplace=True) # check if it makes sense
    candidates = candidates[(~candidates['bin1_id'].isin(b_ID)) & (~candidates['bin2_id'].isin(b_ID))]
    candidates = candidates[(candidates['bin1_id'].isin(gID)) | (candidates['bin2_id'].isin(gID))]
    A_scores = candidates.iloc[:,2:].iloc[:, np.where(cell_types=='A')[0]]
    B_scores = candidates.iloc[:,2:].iloc[:, np.where(cell_types=='B')[0]]
    sigA = g_filter(A_scores) > num_A*.1
    sigB = g_filter(B_scores) > num_B*.1
    sig = np.array([a or b for a, b in zip(sigA,sigB)]).astype(int)
    candidates = pd.concat([candidates.iloc[:,:2], A_scores, B_scores], axis = 1)[sig==1]
    candidates_scores = candidates.iloc[:,2:]
    candidates_scores = qnorm.quantile_normalize(candidates_scores, axis = 1)
    A_scores = candidates_scores.iloc[:, np.where(cell_types=='A')[0]]
    B_scores = candidates_scores.iloc[:, np.where(cell_types=='B')[0]]

    if all(np.var(A_scores, ddof=1, axis=0) > 0) is True and all(np.var(B_scores, ddof=1, axis=0) > 0) is True:
        t_stat, p_val = stats.ttest_ind(A_scores, B_scores, axis=1, equal_var=False)
        out = pd.DataFrame(columns=['bin1_id','bin2_id', 'mean.A','mean.B', 'Tstat', 'Ttest.Pvalue'])
        out[['bin1_id','bin2_id']] = candidates.iloc[:,:2]*binsize
        out['mean.A'] = np.mean(A_scores, axis=1) # should we save mean of quantile-normalized scores?
        out['mean.B'] = np.mean(B_scores, axis=1)
        out['Tstat'] = t_stat
        out['Ttest.Pvalue'] = p_val
        rej, pval_corr = multi.multipletests(out['Ttest.Pvalue'], method='fdr_bh')[:2]
        out['fdr'] = pval_corr
        out['significant'] = (out['fdr'] < fdr_t) & (abs(out['Tstat']) > 2)
        out_filepath = os.path.join(out_dir, 'gap{}_diffs.txt'.format(GAP))
        out.to_csv(out_filepath, header=True, index=None, sep='\t')

    else:
        print("Variance is not >0")
    return 0
   # print("Done with GAP=", GAP, flush=True)
    #return out, final_out

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--indirSHD', action = 'store', required = True, \
                        help = 'SnapHiC-D directory')
    parser.add_argument('-i', '--indirA', action = 'store', required = True, \
                        help = 'input group A directory')
    parser.add_argument('-j', '--indirB', action = 'store', required = True, \
                        help = 'input group B directory')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-c', '--chrom', action = 'store', help = 'chromosome to process', \
                        required = True)
    parser.add_argument('-g', '--gen', action = 'store', help = 'genome - mm10 or hg19', \
                        required = True)
    parser.add_argument('--binsize', type = int, help = 'bin size used for binning the reads', \
                        required = False, default = 1e4)
    parser.add_argument('-n', '--num_cpus', type = int , action = 'store', required = True, \
                        help = 'number of CPUS used for parallel computing')
    parser.add_argument('--fdr_threshold', default = 0.1, type = float, required = False, \
                        help = 'FDR threshold used for candidate peak detection')
    parser.add_argument('--maxi_gap', type = int, help = 'maximum gap between the read pairs', \
                        default = 100, required = False)
    parser.add_argument('--mini_gap', type = int, help = 'minimum gap between the read pairs', \
                        default = 2, required = False)

    return parser




if __name__ == "__main__":
    main()
