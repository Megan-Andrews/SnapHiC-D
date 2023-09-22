import argparse
import h5py
import os
import numpy as np
import pandas as pd
import cooler
import json
from tqdm import tqdm


def main():
    parser = create_parser()
    args = parser.parse_args()
    
    config_path = args.config
    outdir = args.output_dir
    use_k_nbr = args.neighbor_imputation
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    config = open(config_path)
    config = json.load(config)
    temp_dir = config['temp_dir']
    imputed_chromosomes = config['impute_list']
    neighbor_num = config['neighbor_num']
    neighbor_num = neighbor_num if use_k_nbr else 0
    embedding_name = config['embedding_name']
    chr_size_file = config['genome_reference_path']
    resolution = config['resolution']
    

    filelist_path = os.path.join(config['data_dir'], 'filelist.txt')
    filelist = list( pd.read_csv(filelist_path, header = None)[0] )
    outfile_paths = [in2out_filepath(f, outdir) for f in filelist]
    
    bins, offsets = make_bins(chr_size_file, imputed_chromosomes, resolution)
    
    imp_filepaths = {chrom : 
        os.path.join(temp_dir, 
                     '{}_{}_nbr_{}_impute.hdf5'.format(chrom, embedding_name, neighbor_num)) 
        for chrom in offsets.keys()}
    
    imp_files = {chrom : h5py.File(imp_filepaths[chrom], "r") 
                 for chrom in imp_filepaths.keys()}
    
    coords = {}
    for chrom in offsets.keys():
        coords[chrom] = np.array(imp_files[chrom]['coordinates']) + offsets[chrom]
    
    for cell_num in tqdm(range( len(outfile_paths) )):
        pixels = pd.DataFrame()
        for chrom in offsets.keys():
            imp_counts = imp_files[chrom]['cell_{}'.format(cell_num)]
            pixel = pd.DataFrame({'bin1_id': coords[chrom][:,0],
                                'bin2_id': coords[chrom][:,1],
                                'count': imp_counts})
            if pixels.empty:
                pixels = pixel
            else:
                pixels = pd.concat([pixels, pixel], axis = 0)
        cooler.create_cooler(outfile_paths[cell_num],
                        bins = bins,
                        pixels = pixels,
                        dtypes={'count': np.float64})
        
    

def in2out_filepath(infile_path, outdir):
    fn = os.path.basename(infile_path)
    fn = os.path.splitext(fn)[0]
    new_fn = '{}.cool'.format(fn)
    outfile_path = os.path.join(outdir, new_fn)
    return outfile_path
   
def make_bins(chr_size_file, valid_chroms, resolution):
    
    bins = pd.DataFrame()
    offsets = {}
    curr_offset = 0
    with open(chr_size_file, "r") as f:
        for line in f:
            chrom, chrom_size = line.rstrip("\n").split("\t")
            chrom_size = int(chrom_size)
            if chrom in valid_chroms:
                starts = np.arange(0, chrom_size, resolution)
                ends = starts + resolution 
                ends[-1] = chrom_size
                if bins.empty:
                    bins = pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})
                else:
                    bins = pd.concat([bins, pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})], axis = 0)
                offsets[chrom] = curr_offset 
                curr_offset += starts.shape[0]
    bins.index = range(bins.shape[0])
    return bins, offsets
    
def create_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--config', action = 'store', required = True, help = 'path of higashi config')
    parser.add_argument('-o', '--output-dir', action = 'store', required = True, help = 'output directory to save higashi imputed cools')
    parser.add_argument('--neighbor-imputation', action = 'store_true', default = False, \
                        help = 'if set, convert imputations with consideration of neighbors', required = False)
    return parser


if __name__ == "__main__":
    main()
