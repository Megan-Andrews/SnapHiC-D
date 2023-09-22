import argparse
import sys
import data
import scVI_3D
import numpy as np
from joblib import Parallel, delayed
import pandas as pd
import os
import cooler
import torch
import scanpy as sc
import json


def main():
    
    parser = create_parser()
    args = parser.parse_args()
    config = read_config(args.config_path)
    dataset = data.scHiCDataset(config['filepaths_path'], config['format'], 
                                config['cellInfo_path'], config['resolution'], 
                                config['chromosomes'], config['chrom_size_filepath'], 
                                config['max_distance'], config['colnames'], 
                                config['data_save_dir'], config['save_bands'], config['load_bands'])
    dataset.create_band_mat(config['n_cpu'])
    model = scVI_3D.scVI_3D(dataset, config['model_save_dir'], 
                            config['n_latent'], 1, config['gpu_flag'], 
                            config['pool_normalization'])

def read_config(config_path):
    
    config = json.load( open(config_path) )
    assert all([key in config for key in ['format', 'filepaths_path', 'cellInfo_path', 'chrom_size_filepath',
                              'chromosomes', 'resolution', 'model_save_dir', 'n_cpu',
                              'n_latent', 'gpu_flag']]), '"format", "filepaths_path", "cellInfo_path", "chrom_size_filepath", "chromosomes", "resolution", "save_dir" are required parameters.'
    
    if 'save_bands' in config.keys() or 'load_bands' in config.keys():
        assert 'data_save_dir' in config.keys(), '"data_save_dir" should be specified if load from or save to it.'
    if not 'save_bands' in config.keys():
        config['save_bands'] = False
    if not 'load_bands' in config.keys():
        config['load_bands'] = False
    if not 'data_save_dir' in config.keys():
        config['data_save_dir'] = None
    if not 'max_distance' in config.keys():
        config['max_distance'] = None
    if not 'pool_normalization' in config.keys():
        config['pool_normalization'] = False
    if config['format'].endswith('tab') and not 'colnames' in config.keys():
        raise ValueError('"colnames" should be specified for tab formats.')
    if not 'colnames' in config.keys():
        config['colnames'] = None
    
    return config
    
def create_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--config-path', action = 'store', required = True, \
                        help = 'path of a config')
    return parser


if __name__ == "__main__":
    main()