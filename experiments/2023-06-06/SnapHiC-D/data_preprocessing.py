
import time
import numpy as np
import pandas as pd
from scipy.sparse import  eye, csr_matrix, diags
from scipy.sparse.linalg import norm
from scipy.ndimage import gaussian_filter
import cooler
import logging
import os


# -------------------RWR-----------------------------------------------------
def calc_sparsity(matrix):
    row, col = matrix.shape
    sparsity = matrix.nnz / row / col
    return sparsity


def random_walk_cpu(P, rp, tol):
    if rp == 1:
        return P

    _start_time = time.time()
    n_genes = P.shape[0]
    I = eye(n_genes, dtype=np.float32)
    Q = P.copy()
    for i in range(30):
        Q_new = P.dot(Q * (1 - rp) + rp * I)
        delta = norm(Q - Q_new)
        Q = Q_new.copy()
        sparsity = calc_sparsity(Q)
        _end_time = time.time()
        logging.debug(
            f"Iter {i + 1} takes {(_end_time - _start_time):.3f} seconds. "
            f'Loss: {delta:.3f}; Sparsity: {sparsity:.3f}', P.dtype, Q.dtype)
        if delta < tol:
            break
    return Q


def impute_chromosome(scool_url,
                      chrom,
                      resolution,
                      logscale=False,
                      pad=1,
                      std=1,
                      rp=0.5,
                      tol=0.01,
                      window_size=500000000,
                      step_size=10000000,
                      output_dist=500000000,
                      min_cutoff=0):
    cell_cool = cooler.Cooler(scool_url)
    A = cell_cool.matrix(balance=False, sparse=True).fetch(chrom)
    #print(A)
    n_bins = A.shape[0]
    ws = int(window_size // resolution)
    ss = int(step_size // resolution)

    # log transform
    if logscale:
        A.data = np.log2(A.data + 1)

    # Remove diagonal before convolution
    A = A - diags(A.diagonal())

    # Gaussian convolution and
    start_time = time.time()
    if pad > 0:
        # full matrix step
        A = gaussian_filter((A + A.T).astype(np.float32).toarray(),
                            std, order=0, mode='mirror', truncate=pad)
        A = csr_matrix(A)
    else:
        A = A + A.T
    end_time = time.time()
    logging.debug(f'Convolution takes {end_time - start_time:.3f} seconds')

    # Remove diagonal before RWR
    A = A - diags(A.diagonal())

    # Random Walk with Restart
    start_time = time.time()
    if ws >= n_bins or rp == 1:
        B = A + diags((A.sum(axis=0).A.ravel() == 0).astype(int))
        d = diags(1 / B.sum(axis=0).A.ravel())
        P = d.dot(B).astype(np.float32)
        E = random_walk_cpu(P, rp, tol)
    else:
        # if the chromosome is too large, compute by chunks
        idx = (np.repeat(np.arange(ws), ws), np.tile(np.arange(ws), ws))
        idxfilter = (np.abs(idx[1] - idx[0]) < (output_dist // resolution + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
        # first filter
        idxfilter = ((idx[0] + idx[1]) < (ws + ss))
        idx1 = (idx[0][idxfilter], idx[1][idxfilter])
        mask1 = csr_matrix((np.ones(len(idx1[0])), (idx1[0], idx1[1])),
                           (ws, ws))
        # last filter
        idxfilter = ((idx[0] + idx[1]) >= (
                (n_bins - ws) // ss * 2 + 1) * ss + 3 * ws - 2 * n_bins)
        idx2 = (idx[0][idxfilter], idx[1][idxfilter])
        mask2 = csr_matrix((np.ones(len(idx2[0])), (idx2[0], idx2[1])),
                           (ws, ws))
        # center filter
        idxfilter = np.logical_and((idx[0] + idx[1]) < (ws + ss),
                                   (idx[0] + idx[1]) >= (ws - ss))
        idx0 = (idx[0][idxfilter], idx[1][idxfilter])
        mask0 = csr_matrix((np.ones(len(idx0[0])), (idx0[0], idx0[1])),
                           (ws, ws))

        start_time = time.time()
        E = csr_matrix(A.shape, dtype=np.float32)
        for ll in [x for x in range(0, n_bins - ws, ss)] + [n_bins - ws]:
            B = A[ll:(ll + ws), ll:(ll + ws)]
            B = B + diags((B.sum(axis=0).A.ravel() == 0).astype(int))
            d = diags(1 / B.sum(axis=0).A.ravel())
            P = d.dot(B).astype(np.float32)
            Etmp = random_walk_cpu(P, rp, tol)
            if ll == 0:
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask1)
            elif ll == (n_bins - ws):
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask2)
            else:
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask0)
    logging.debug(f'RWR takes {time.time() - start_time:.3f} seconds')

    # Normalize
    start_time = time.time()
    E += E.T
    d = E.sum(axis=0).A.ravel()
    d[d == 0] = 1
    b = diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    logging.debug(f'SQRTVC takes {time.time() - start_time:.3f} seconds')

    start_time = time.time()
    # mask the lower triangle of E
    # TODO This part is MEM intensive, the mask below can be combined with the chunk mask above
    idx = np.triu_indices(E.shape[0], 0)
    if (output_dist // resolution + 1) < n_bins:
        # longest distance filter mask
        idxfilter = ((idx[1] - idx[0]) < (output_dist // resolution + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
    mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])),
                      E.shape,
                      dtype=np.float32)
    E = E.tocsr().multiply(mask)
    logging.debug(f'Filter takes {time.time() - start_time:.3f} seconds')

    # TODO put this part inside RWR, before normalize
    # min_cutoff = tol/
    # Make values < min_cutoff to 0
    if min_cutoff > 0:
        s_before = calc_sparsity(E)
        E = E.multiply(E > min_cutoff)
        s_after = calc_sparsity(E)
        logging.debug(f'Mask values smaller than {min_cutoff}. Sparsity before {s_before:.3f}, after {s_after:.3f}')

    # save to file
    #write_coo(output_path, E)

    # Replaced write_coo with cooler function
     # Determine the number of bins
    num_bins = E.shape[0]

    # Create a DataFrame with bin coordinates
    # bins = pd.DataFrame({
    #     'chrom': [chrom] * num_bins,
    #     'start': np.arange(0, num_bins * resolution, resolution),
    #     'end': np.arange(resolution, num_bins * resolution + resolution, resolution)
    # })
    coo_matrix = E.tocoo()
    
    # Create DataFrame from the coo_matrix data
    pixels = pd.DataFrame({
        'bin1_id': coo_matrix.row,
        'bin2_id': coo_matrix.col,
        'count': coo_matrix.data
    })
    #print(E)
    #print(bins)
    #print(pixels)
    
    ## Check if cool file looks correct
    #c= cooler.Cooler(output_path)
    #T = c.matrix(balance=False, sparse=True).fetch(chrom)
    #print(T)

   
    return pixels, num_bins


def impute(coolfile, chrom, resolution):
    # Provide the necessary parameters
    scool_url = coolfile
    # Call the impute_chromosome function
    return impute_chromosome(
        scool_url=scool_url,
        chrom=chrom,
        resolution=resolution,
        logscale=False,
        pad=1,
        std=1,
        rp=0.5,
        tol=0.01,
        window_size=500000000,
        step_size=10000000,
        output_dist=500000000,
        min_cutoff=0
    )

def coarsen(inputfile, outputfile, chrom, factor):
    # c = cooler.Cooler(inputfile)
    # m = c.matrix(balance=False, as_pixels=True).fetch(chrom)
    base_uri = inputfile
    output_uri = outputfile
    chunksize = 10000
    cooler.coarsen_cooler(base_uri, output_uri, factor, chunksize)


"""def temp_normalize():
    def valid_bin_per_gap(num_bins, valid_bins, gap):
        cnt = 0
        for b in range(num_bins-gap):
            if (b in valid_bins) & ((b+gap) in valid_bins):
                cnt = cnt + 1
        return cnt

    def upper_trim(x, trim, n_counts):
        #print(type(x))
        
        x = np.sort(x)
        trim_length = int(np.ceil(len(x) * (1 - trim)))
        trimmed_x = x[:trim_length]
        if len(trimmed_x) == 0:
            return 0,0
        trimmed_mean = np.sum(trimmed_x)/n_counts
        trimmed_std = np.sum(np.power(trimmed_x-trimmed_mean, 2))/n_counts
        return np.mean(trimmed_x), np.std(trimmed_x)

    #filtered_regions = pd.read_csv('../ext/hg19_filter_regions.txt', sep = "\t", header = None)
    #filtered_regions = filtered_regions[filtered_regions.iloc[:,0]=='chr22']
    #filtered_regions['bin_num'] = (filtered_regions.iloc[:,1]/100000).astype(int)
    #filtered_bins = np.unique(filtered_regions['bin_num'])
    raw_cool = cooler.Cooler('/project/compbio-lab/scHi-C/100kb_imputed_cool/181218_21yr_2_A12_AD004_OPC_100kb_contacts/181218_21yr_2_A12_AD004_OPC_100kb_contacts_coarse.cool')
    raw_pixel = raw_cool.matrix(balance = False, as_pixels = True).fetch('chr22')
    raw_pixel = raw_pixel[raw_pixel['bin1_id']!=raw_pixel['bin2_id']]
    chr22_offset = raw_cool.offset('chr22')
    valid_bins = np.union1d(raw_pixel['bin1_id'],raw_pixel['bin2_id']) - chr22_offset
    c = cooler.Cooler('/project/compbio-lab/scHi-C/100kb_imputed_cool/181218_21yr_2_A12_AD004_OPC_100kb_contacts/181218_21yr_2_A12_AD004_OPC_100kb_contacts_imputed.cool')
    num_bins = int(c.chromsizes['chr22']/c.binsize)
    #valid_bins = [b for b in range(num_bins) if not b in filtered_bins]
    m = c.matrix(balance=False, as_pixels=True).fetch('chr22')
    m = m[(m['bin1_id'].isin(valid_bins)) & (m['bin2_id'].isin(valid_bins))]
    m['gap'] = m['bin2_id'] - m['bin1_id']  
    gap_means, gap_stds = [], []
    for gap in range(num_bins):
        gap_counts = m[m['gap']==gap]['count']
        n_counts = valid_bin_per_gap(num_bins, valid_bins, gap)
        mean,std = upper_trim(gap_counts, 0.01, n_counts)
        gap_means.append(mean)
        gap_stds.append(std)
    gap_means, gap_stds = np.array(gap_means), np.array(gap_stds)
    m['gap_mean'] = gap_means[m['gap']]
    m['gap_std'] = gap_stds[m['gap']]
    m['z_score'] = ((m['count'] - m['gap_mean'])/m['gap_std'])*(m['gap_std']>0.000001)
"""
def normalize(coarsefile, imputed_pixel, num_bins, chrom):
    def valid_bin_per_gap(num_bins, valid_bins, gap):
        cnt = 0
        for b in range(num_bins-gap):
            if (b in valid_bins) & ((b+gap) in valid_bins):
                cnt = cnt + 1
        return cnt

    # x: the array being trimmed
    # trim: the trim percentage
    def upper_trim(x, trim, n_counts):
        x = np.sort(x)
        trim_length = int(np.ceil(len(x) * (1 - trim)))
        trimmed_x = x[:trim_length]
        if len(trimmed_x) == 0:
            return 0,0
        trimmed_mean = np.sum(trimmed_x)/n_counts
        trimmed_std = np.sqrt(np.sum(np.power(trimmed_x-trimmed_mean, 2))/n_counts)
        return trimmed_mean, trimmed_std

    # filtered_regions = pd.read_csv('./ext/hg19_filter_regions.txt', sep = "\t", header = None)
    # filtered_regions = filtered_regions[filtered_regions.iloc[:,0]==chrom]
    # filtered_regions['bin_num'] = (filtered_regions.iloc[:,1]/100000).astype(int)
    # filtered_bins = np.unique(filtered_regions['bin_num'])

    raw_cool = cooler.Cooler(coarsefile)
    raw_pixel = raw_cool.matrix(balance = False, as_pixels = True).fetch('chr22')
    raw_pixel = raw_pixel[raw_pixel['bin1_id'] != raw_pixel['bin2_id']]
    offset = raw_cool.offset(chrom)
    valid_bins = np.union1d(raw_pixel['bin1_id'],raw_pixel['bin2_id']) - offset

    m = imputed_pixel
    m = m[(m['bin1_id'].isin(valid_bins)) & (m['bin2_id'].isin(valid_bins))]
    m['gap'] = m['bin2_id'] - m['bin1_id']
    
    gap_means, gap_stds = [], []
    for gap in range(num_bins):
        gap_counts = m[m['gap']==gap]['count']
        n_counts = valid_bin_per_gap(num_bins, valid_bins, gap)
        mean,std = upper_trim(gap_counts, 0.01, n_counts)
        gap_means.append(mean)
        gap_stds.append(std)
    
    gap_means, gap_stds = np.array(gap_means), np.array(gap_stds)

    m['gap_mean'] = gap_means[m['gap']]
    m['gap_std'] = gap_stds[m['gap']]
    m['count'] = ((m['count'] - m['gap_mean'])/m['gap_std'])*(m['gap_std']>0.000001)
    
    #print(m)
    return m
    

def main():
    chrom = 'chr22'
    resolution = 100000
    curr_resolution = 10000

    factor = resolution/curr_resolution

    inputdir = "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_cool/"
    coarsedir = "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"
    outdir = "/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/ODC_100kb_imputed_cool/"

    file_list = "/home/maa160/SnapHiC-D/experiments/2023-06-06/SnapHiC-D/file_lists/ODC_file_listA.txt"
    
    file_no = 1
    with open(file_list, "r") as file:
    # Iterate through each line in the file
        for line in file:
            print(file_no, line.strip())
            run(filename = line.strip(), chrom=chrom, resolution=resolution, factor= factor, inputdir=inputdir, coarsedir=coarsedir, outdir=outdir)
            file_no += 1

        

def run(filename, chrom, resolution, factor, inputdir, coarsedir, outdir):
    
    inputfile = inputdir + filename
    coarsefile = coarsedir + filename.replace("_10kb_", "_100kb_")
    outputfile = outdir + filename.replace("_10kb_", "_100kb_").replace(".cool","_imputed.cool")
   
    if coarsefile not in os.listdir(coarsedir):
        coarsen(inputfile, coarsefile, chrom, factor)

    if outputfile not in os.listdir(outdir):
        imputed_pixels, num_bins= impute(coarsefile, chrom, resolution)
        normalized_pixels = normalize(coarsefile, imputed_pixels, num_bins, chrom)

        bins = pd.DataFrame({
            'chrom': [chrom] * num_bins,
            'start': np.arange(0, num_bins * resolution, resolution),
            'end': np.arange(resolution, num_bins * resolution + resolution, resolution)
        })

        cooler.create_cooler(outputfile, 
                            bins=bins, 
                            pixels=normalized_pixels,
                            ordered=True,
                            # columns = ['z_score'],
                            dtypes={'count': np.float64}) # TODO: should this be float32? and is the column 'z_score'?


if __name__ == "__main__":
    main() 