import argparse
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
import gc
import cooler
from scipy.sparse import eye, coo_matrix, csc_matrix
import re
import time
import sys

def main():
    parser = create_parser()
    args = parser.parse_args()
    #file_paths = []
    #with open(args.filelist, 'r') as f:
    #    for line in f:
    #        file_paths.append(line.rstrip('\n'))
    window_size = args.window_size
    step_size = int(window_size/2)
    chromosome_lengths = read_chrom_lens(args.chrom_lens)
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    #for filepath in file_paths:
    log = open(args.log_path, "w")
    process_cell(args.input_filepath, args.format, args.extension, args.output_dir, args.resolution, chromosome_lengths, args.chrom_columns, args.pos_columns, window_size, step_size, args.upper_distance, args.rp, log)

def process_cell(filepath, format, extension, output_dir, resolution, chromosome_lengths, chrom_columns, pos_columns, window_size, step_size, upper_distance, rp, log):


    start_time = time.time()
    bins, offset = make_bins(chromosome_lengths, resolution)
    all_imp_edgelist = pd.DataFrame()
    if format == 'cooler':
        c = cooler.Cooler(filepath)
        for chrom in chromosome_lengths.keys():
            edgelist = c.matrix(balance = False, as_pixels=True).fetch(chrom)
            edgelist.iloc[:,[0,1]] -= c.offset(chrom)
            log.write('imputing chromosome {} from {}\n'.format(chrom, filepath))
            imp_edgelist = rwr(edgelist, chromosome_lengths[chrom], resolution, window_size, step_size, upper_distance, rp)
            imp_edgelist.iloc[:,[0,1]] += offset[chrom]
            all_imp_edgelist = pd.concat([all_imp_edgelist, imp_edgelist], axis = 0)

    elif format == 'tab-sep':
        all_edgelist = pd.read_csv(filepath, sep = "\t", header = None)
        for chrom in chromosome_lengths.keys():
            edgelist = all_edgelist[(all_edgelist.iloc[:,chrom_columns[0]] == chrom) & (all_edgelist.iloc[:,chrom_columns[1]] == chrom)]
            edgelist = edgelist.iloc[:,pos_columns]
            edgelist = (edgelist/resolution).astype(int)
            edgelist.columns = ['bin1_id', 'bin2_id']
            edgelist = edgelist.groupby(['bin1_id','bin2_id']).size().to_frame('count').reset_index()
            log.write('imputing chromosome {} from {}\n'.format(chrom, filepath))
            imp_edgelist = rwr(edgelist, chromosome_lengths[chrom], resolution, window_size, step_size, upper_distance, rp)
            imp_edgelist.iloc[:,[0,1]] += offset[chrom]
            all_imp_edgelist = pd.concat([all_imp_edgelist, imp_edgelist], axis = 0)
    else:
        log.write('cool and tab-sep input formats are available.\n')
        return

    reg = r"(?P<fname>.*)(%s)" % extension
    filename = re.match(reg, os.path.basename(filepath)).groupdict()['fname']
    out_filepath = os.path.join(output_dir, '{}_rwr.cool'.format(filename))
    log.write('writing results to {}\n'.format(out_filepath))
    cooler.create_cooler(out_filepath,
                    bins=bins,
                    pixels=all_imp_edgelist,
                    ordered=True,
                    dtypes={'count': np.float64})
    log.write("--- it took %s seconds to process one single cell ---\n" % (time.time() - start_time))
    return


def rwr(edges, chrom_len, resolution, window_size, step_size, distance_threshold, rp):

    edges.columns = ['bin1_id', 'bin2_id', 'count']
    valid_bins = np.unique(list(edges['bin1_id'].values) + list(edges['bin2_id'].values))
    chrom_size = int(np.ceil(chrom_len/resolution))
    ws, ss = int(window_size/resolution), int(step_size/resolution)
    distance_threshold = int(distance_threshold/resolution)
    #r = coo_matrix(((0,), ((0,), (0,))), shape = (chrom_size,chrom_size)).todense()
    rwr_edges = pd.DataFrame()
    for window_start in tqdm(range(0, chrom_size, ss)):
        window_end = window_start + ws
        window_edges = edges[(edges['bin1_id'] >= window_start) & (edges['bin2_id'] < window_end)]
        neighbor_edges = pd.DataFrame({'bin1_id':list(range(window_start, min(chrom_size-1, window_end - 1))), 'bin2_id': list(range(window_start+1, min(chrom_size, window_end)))})
        window_edges = pd.concat([neighbor_edges[['bin1_id', 'bin2_id']], window_edges[['bin1_id', 'bin2_id']]], axis = 0)
        if window_edges.shape[0] == 0:
            continue
        window_edges.loc[:,'count'] = 1
        g = get_stochastic_matrix_from_edgelist(window_edges)
        window_edges_imp = solve_rwr_iterative(g, rp)
        window_edges_imp = window_edges_imp[(window_edges_imp['i'] + window_edges_imp['j'] > ss) & (window_edges_imp['i'] + window_edges_imp['j'] < ws + ss)]
        window_edges_imp = window_edges_imp[(window_edges_imp['j'] - window_edges_imp['i'] <= distance_threshold)]
        window_edges_imp['i'] += window_start
        window_edges_imp['j'] += window_start
        window_edges_imp.columns = ['bin1_id', 'bin2_id', 'count']
        window_edges_imp = window_edges_imp[window_edges_imp['bin2_id'] >= window_edges_imp['bin1_id']]
        rwr_edges = pd.concat([rwr_edges, window_edges_imp], axis = 0)
        #partial_r = coo_matrix((window_edges_imp['v'], (window_edges_imp['i'], window_edges_imp['j'])), shape = (chrom_size, chrom_size))
        #r += partial_r
    #r = coo_matrix(r)
    #rwr_edges = pd.DataFrame({'bin1_id':r.row, 'bin2_id': r.col, 'count': r.data})
    #rwr_edges = rwr_edges[rwr_edges['bin2_id'] >= rwr_edges['bin1_id']]
    rwr_edges = rwr_edges[(rwr_edges['bin1_id'].isin(valid_bins)) & (rwr_edges['bin2_id'].isin(valid_bins))]
    return rwr_edges

def get_stochastic_matrix_from_edgelist(edgelist):
    g = nx.from_pandas_edgelist(edgelist, source = 'bin1_id', target = 'bin2_id', edge_attr = ['count'], create_using = nx.Graph())
    degrees = np.array([g.degree(i) for i in g.nodes()])
    m = csc_matrix(nx.adjacency_matrix(g).astype(float))
    m.data /= degrees[m.indices] #stochastic matrix
    del g, degrees
    return m

def solve_rwr_iterative(stoch_matrix, alpha = 0.05, max_iter = None):
    max_iter = max_iter if max_iter else float("inf")
    # what is this for?
    gc.collect()
    I = eye(stoch_matrix.shape[0], dtype=np.float32)
    delta = float("inf")
    A = I.copy()
    counter = 0
    while delta > 1e-6 and counter < max_iter:
        counter += 1
        Aold = A.copy()
        A = (1-alpha) * stoch_matrix * Aold + (alpha) * I
        delta = (abs(A - Aold)).max()

    A += A.transpose()
    A = coo_matrix(A)
    df = pd.DataFrame({'i':A.row, 'j': A.col, 'count': A.data})
    del A
    return df

'''
def process_chromosome(all_edgelist, chrom, chroms_length, offset, resolution, chrom_columns, pos_columns,
                      window_size, step_size, upper_distance, rp):
    edgelist = all_edgelist[(all_edgelist.iloc[:,chrom_columns[0]] == chrom) & (all_edgelist.iloc[:,chrom_columns[1]] == chrom)]
    edgelist = edgelist.iloc[:,pos_columns]
    edgelist = (edgelist/resolution).astype(int)
    edgelist.columns = ['bin1_id', 'bin2_id']
    edgelist = edgelist.groupby(['bin1_id','bin2_id']).size().to_frame('count').reset_index()
    print('imputing chromosome {}'.format(chrom))
    imp_edgelist = snapHiC_preprocessing.rwr(edgelist, chroms_length[chrom], resolution, window_size, step_size, upper_distance, rp)
    imp_edgelist.iloc[:,[0,1]] += offset[chrom]
    return imp_edgelist
'''


def read_chrom_lens(f):
    chrom_lens = {}
    with open(f, 'r') as f:
        for line in f:
            chrom, size = line.rstrip('\n').split('\t')
            chrom_lens[chrom] = int(size)
    return chrom_lens

def make_bins(chromosome_lengths, resolution):
    bins = pd.DataFrame()
    offset = {}
    for chrom in chromosome_lengths.keys():
        offset[chrom] = bins.shape[0]
        starts = np.arange(0,chromosome_lengths[chrom],resolution)
        ends = starts[1:]
        ends = np.insert(ends, ends.shape, chromosome_lengths[chrom])
        bins = pd.concat([bins, pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})], axis = 0)
    bins.index = range(bins.shape[0])
    return bins, offset

def create_parser():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-f', '--filelist', action = 'store', required = True, \
    #                    help = 'path of the file including filepaths to impute')
    parser.add_argument('-i', '--input-filepath', action = 'store', required = True, help = 'path of input file to impute')
    parser.add_argument('-o', '--output-dir', action = 'store', required = True, help = 'output directory to save imputed contact maps')
    parser.add_argument('-r', '--resolution', action = 'store', required = True, type = int, help = "resolution of contact maps")
    parser.add_argument('-l', '--chrom-lens', action = 'store', required = True, help = 'path to the chromosome lengths file')
    parser.add_argument('--format', action = 'store', required = True, help = 'format of input files: cooler or tab-sep')
    parser.add_argument('--window-size', action = 'store', default = 10000000, type = int, help = 'window size for window-slicing rwr')
    parser.add_argument('--rp', action = 'store', default = 0.05, type = float, help = 'restart probability')
    parser.add_argument('--chrom-columns', action = 'store', required = False, nargs = 2, type = int, help = 'two integer column numbers for chromosomes')
    parser.add_argument('--pos-columns', action = 'store', required = False, nargs = 2, type = int, help = 'two integer column numbers for positions')
    parser.add_argument('--extension', action = 'store', required = True, help = 'extension of input files like .cool and .txt.gz')
    parser.add_argument('--upper-distance', action = 'store', default = 2000000, type = int, help = 'maximum distance between bin pairs to store')
    parser.add_argument('--log-path', action = 'store', required = True, help = 'log path')

    return parser


if __name__ == "__main__":
    main()
