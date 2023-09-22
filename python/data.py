import os
from typing import Literal
from joblib import Parallel, delayed
from tqdm import tqdm 
import time

import pandas as pd
import numpy as np

import cooler

class scHiCDataset:

    def __init__(self, 
                 filelist_path: str,
                 format: Literal["tab", "binned_tab", "cool"], 
                 cellinfo_path: str, 
                 resolution: int, 
                 chromosomes: list, 
                 chrom_sizes_filepath: str,
                 max_distance: int = None,
                 colnames: list = None,
                 save_dir: str = None,
                 save_band_mats: bool = False,
                 load_existing_band_mats: bool = False):

        self.filelist_path = filelist_path
        self.filepaths = []
        with open(self.filelist_path, 'r') as f:
            for line in f:
                self.filepaths.append(line.rstrip("\n"))
                
        self.cell_num = len(self.filepaths)
        
        self.format = format 
        if not self.format in ["tab", "binned_tab", "cool"]:
            raise ValueError("'tab', 'binned_tab', and 'cool' formats are available.")
        if self.format in ["tab", "binned_tab"]:
            if colnames is None:
                raise ValueError("column names including 'chr1', 'pos1', 'chr2', 'pos2', and optional 'count' should be specified for 'tab' and 'binned_tab' formats.")
            else:
                coord_colnames = ['chr1', 'pos1', 'chr2', 'pos2']
                if all([n in colnames for n in coord_colnames]):
                    self.colnames = colnames
                else:
                    raise ValueError("column names including 'chr1', 'pos1', 'chr2', 'pos2', and optional 'count' should be specified for 'tab' and 'binned_tab' formats.")
        
        self.cellinfo_path = cellinfo_path
        self.resolution = resolution
        self.chromosomes = chromosomes
        self.chrom_sizes_filepath = chrom_sizes_filepath
        self.chrom_sizes = read_chr_sizes(self.chrom_sizes_filepath)
        self.num_bins = {}
        for chrom in self.chromosomes:
            self.num_bins[chrom] = int(np.ceil(self.chrom_sizes[chrom] / self.resolution))
            
        self.max_distance = max_distance
        self.max_diags = {}
        for chrom in self.chromosomes:
            if self.max_distance == None:
                self.max_diags[chrom] = self.num_bins[chrom]
            else:
                self.max_diags[chrom] = min(int(self.max_distance / self.resolution) + 1,
                                            self.num_bins[chrom])
       
        self.save_dir = save_dir
        self.load_existing_band_mats = load_existing_band_mats
        self.save_band_mats = save_band_mats
        if (self.load_existing_band_mats):
            for chrom in self.chromosomes:
                for diag in np.arange(1, self.max_diags[chrom]):
                    band_mat_filepath = os.path.join(self.save_dir, '{}_band{}.npy'.format(chrom, diag))
                    assert os.path.exists(band_mat_filepath), 'Stored band matrices do not exist. Do not set load_existing_band_mats True.'
            print('Found all existing band matrices.')
            self.save_band_mats = False
        if not os.path.exists(self.save_dir):
            os.mkdir(self.save_dir)
        
        self.all_band_mats = None
        
        self.batch_names, self.batch_nums, self.cell_types = self.read_cellinfo()
        self.batchFlag = False
        if isinstance(self.batch_names, np.ndarray):
            self.batchFlag = True
            self.batch_types = list(set(self.batch_names))
        
    def load_cell_bands(self, cell_id, filepath):
      
        if self.format == 'tab':
            pixel = pd.read_csv(filepath, sep = "\t", header = None, names = self.colnames)
            pixel = pixel[pixel['chr1'] == pixel['chr2']]
            pixel = pixel[pixel['chr1'].isin(self.chromosomes)]
            pixel[['bin1_id', 'bin2_id']] = (pixel[['pos1', 'pos2']] / self.resolution).astype(int)
            if not 'count' in pixel.columns:
                pixel['count'] = 1
            pixel = pixel.groupby(['chr1', 'bin1_id', 'bin2_id']).aggregate({'count': 'sum'}).reset_index()
            
        elif self.format == 'binned_tab':
            pixel = pd.read_csv(filepath, sep = "\t", header = None, names = self.colnames)
            pixel = pixel[pixel['chr1'] == pixel['chr2']]
            pixel = pixel[pixel['chr1'].isin(self.chromosomes)]
            pixel[['bin1_id', 'bin2_id']] = pixel[['pos1', 'pos2']]
            if not 'count' in pixel.columns:
                pixel['count'] = 1
            
        else:
            c = cooler.Cooler(filepath)
            pixel = []
            for chrom in self.chromosomes:
                pixel_tmp = c.matrix(balance = False, as_pixels = True).fetch(chrom)
                pixel_tmp[['bin1_id', 'bin2_id']] -= c.offset(chrom) 
                pixel_tmp['chr'] = chrom
                pixel.append(pixel_tmp)
            pixel = pd.concat(pixel, axis = 0)
            
        
        pixel['cell_id'] = cell_id
        if not 'chr' in pixel.columns:
            pixel.rename(columns = {'chr1': 'chr'}, inplace = True)
        pixel = pixel[['cell_id', 'chr', 'bin1_id', 'bin2_id', 'count']]
        pixel['diag'] = abs(pixel['bin2_id'] - pixel['bin1_id'])
        
        chrom_band_vecs = {}
        for chrom in self.chromosomes:
            chrom_band_vecs[chrom] = {}
            chrom_pixel = pixel[pixel['chr'] == chrom]
            for diag in np.arange(1, self.max_diags[chrom]):
                chrom_diag_pixel = chrom_pixel[chrom_pixel['diag'] == diag]
                counts = np.zeros(self.num_bins[chrom] - diag)
                counts[chrom_diag_pixel['bin1_id']] = chrom_diag_pixel['count']
                chrom_band_vecs[chrom][diag] = counts
                
        return chrom_band_vecs

    
    def create_band_mat(self, n_cpus):
        
        self.all_band_mats = {}
        if self.load_existing_band_mats:
            
            st = time.time()
            print('Start loading existing band matrices.')
            for chrom in self.chromosomes:
                self.all_band_mats[chrom] = {}
                for diag in np.arange(1, self.max_diags[chrom]):
                    band_mat_filepath = os.path.join(self.save_dir, '{}_band{}.npy'.format(chrom, diag))
                    self.all_band_mats[chrom][diag] = np.load(band_mat_filepath)
            print('It took {} s to load band matrices'.format(time.time() - st))
            
        else:
            
            st = time.time()
            print('Start loading cell bands.')
            cells_bands = Parallel(n_jobs = n_cpus)(delayed(self.load_cell_bands)
                                                    (cell_id, filepath)
                                                    for (cell_id, filepath) in 
                                                    tqdm(enumerate(self.filepaths)))
            print('It took {} s to load cell bands'.format(time.time() - st))

            st = time.time()
            print('Start creating band matrices.')
            
            for chrom in self.chromosomes:
                
                self.all_band_mats[chrom] = {}
                
                for diag in np.arange(1, self.max_diags[chrom]):
                    band_vecs = []
                    for i in range(self.cell_num):
                        band_vecs.append(cells_bands[i][chrom][diag])
                    self.all_band_mats[chrom][diag] = np.vstack(band_vecs)
            print('It took {} s to create band matrices.'.format(time.time() - st))
        
            if self.save_band_mats:
                print('Saving band matrices.')
                for chrom in self.all_band_mats:
                    for diag in self.all_band_mats[chrom]:
                        outfile_path = os.path.join(self.save_dir, '{}_band{}'.format(chrom, diag))
                        np.save(outfile_path, self.all_band_mats[chrom][diag])
        
   
        
    def load_HiC_maps(self):

        filepaths = []
        with open(self.filelist_path, "r") as f:
            for line in f:
                filepaths.append(line.rstrip("\n"))
                            
        contact_map_list = []
        library_sizes = []
        print("reading contact maps")
        for filepath in tqdm(filepaths):
            library_size, adj_mat = read_binned_tab_hic(filepath, self.chrom, self.chrom_size, 
                                                        self.resolution, self.filter_bin_id)
            adj_mat = torch.from_numpy(adj_mat)
            #adj_mat = torch.from_numpy(adj_mat[None,:])
            contact_map_list.append(adj_mat)
        print('finished reading contact maps')
        return contact_map_list, library_sizes

    def read_cellinfo(self):

        cellinfo = pd.read_csv(self.cellinfo_path, sep = "\t")
        
        if 'batch_name' in cellinfo.columns:
            batches = cellinfo['batch_name'].values
            unique_batches = np.unique(batches)
            batch_dict = {}
            for b, batch in enumerate(unique_batches):
                batch_dict[batch] = b
            batches_num = [batch_dict[b] for b in batches]
        else:
            batches = None
            batches_num = None
        
        if 'cell_type' in cellinfo.columns:
            cell_types = cellinfo['cell_type'].values
        else:
            cell_types = None
            
        return batches, batches_num, cell_types

    def get_band_lps(self, chrom, band):
        
        bin1_id = np.arange(0, self.num_bins[chrom] - band)
        bin2_id = np.arange(band, self.num_bins[chrom])
        return np.array([bin1_id, bin2_id])
        
        
def read_chr_sizes(chr_size_file):
    chr_sizes = {}
    with open(chr_size_file) as f:
        for line in f:
            chr_name, chr_size = line.rstrip("\n").split("\t")
            chr_sizes[chr_name] = int(chr_size)
    return chr_sizes

