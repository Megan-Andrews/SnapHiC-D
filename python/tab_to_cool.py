# input: config with parameters: filelist path including paths of non-binned tab-sep files, 
# columns names (including chrom1, pos1, chrom2, pos2), chromosome size file, 
# list of chromosomes, resolution

import argparse
import os
import cooler
import json 
import pandas as pd
import numpy as np
from tqdm import tqdm

def main():
    
    parser = create_parser()
    args = parser.parse_args()
    config = json.load( open(args.config_path, 'r') )
    if not os.path.exists(config['outdir']):
        os.mkdir(config['outdir'])
    bins, offsets = make_bins(config['chr_size_file'], config['chromosomes'], config['resolution'])
    with open( config['filelist'], 'r' ) as f:
        for line in tqdm(f):
            tab2cool(line.rstrip('\n'), config['outdir'], config['colnames'],
                     config['resolution'], bins, offsets)

def tab2cool(tab_filepath, outdir, colnames, resolution, bins, offsets):
    
    tab_df = pd.read_csv(tab_filepath, sep = "\t", header = None, names = colnames)
    tab_df = tab_df[['chrom1', 'pos1', 'chrom2', 'pos2']]
    tab_df = tab_df[tab_df['chrom1'] == tab_df['chrom2']]
    pixels = pd.DataFrame()
    for chrom in offsets.keys():
        tab_df_ = tab_df[tab_df['chrom1'] == chrom].copy()
        tab_df_['bin1_id'] = (tab_df_['pos1'] / resolution).astype(int)
        tab_df_['bin2_id'] = (tab_df_['pos2'] / resolution).astype(int)
        bin1_id, bin2_id = tab_df_['bin1_id'], tab_df_['bin2_id']
        tab_df_['bin1_id'] = np.minimum(bin1_id, bin2_id)
        tab_df_['bin2_id'] = np.maximum(bin1_id, bin2_id)
        tab_df_ = tab_df_[['bin1_id', 'bin2_id']]
        pixel = tab_df_.groupby(['bin1_id', 'bin2_id']).size().reset_index()
        pixel.rename(columns = {0: 'count'}, inplace = True)
        pixel[['bin1_id', 'bin2_id']] += offsets[chrom]
        if pixels.empty:
            pixels = pixel
        else:
            pixels = pd.concat([pixels,pixel], axis = 0)
    filename = os.path.splitext( os.path.basename(tab_filepath) )[0]
    out_filepath = os.path.join(outdir, '{}.cool'.format(filename))
    cooler.create_cooler(out_filepath,
                    bins=bins,
                    pixels=pixels,
                    ordered=True,
                    dtypes={'count': np.int32})
    
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

    parser.add_argument('-c', '--config-path', action = 'store', required = True, help = 'config path')
   
    return parser


if __name__ == "__main__":
    main()