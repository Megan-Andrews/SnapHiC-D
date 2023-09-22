# this script is based on https://github.com/yezhengSTAT/scVI-3D/tree/master


import os
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import cooler
import shutil


class scVI_3D:
    
    def __init__(self,
                 scHiCDataset,
                 out_dir: str,
                 nLatent: int,
                 n_cpus: int,
                 gpuFlag: bool,
                 pool_normalization: bool = False,
                 ):
        
        self.dataset = scHiCDataset
        if self.dataset.all_band_mats == None:
            raise ValueError('band matrices should be loaded to scHiCDataset.')
        self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        else:
            shutil.rmtree(self.out_dir)
            os.mkdir(self.out_dir)
        
        self.nLatent = nLatent
        self.n_cpus = n_cpus
        self.gpuFlag = gpuFlag
        
        
        self.pool_normalization = pool_normalization
        self.pools_idx = {}
        if self.pool_normalization:
            for chrom in self.dataset.chromosomes:
                self.pools_idx[chrom] = get_pool_idx(self.dataset.max_diags[chrom])
        else:
            for chrom in self.dataset.chromosomes:
                self.pools_idx[chrom] = [[b] for b in np.arange(1,self.dataset.max_diags[chrom])]
                
        
        self.latent_dir = os.path.join(self.out_dir, 'latents')
        os.mkdir(self.latent_dir)
        
        self.unique_cell_types = list(set(self.dataset.cell_types))
        
        DCC_dir = os.path.join(self.out_dir, 'lvm_DCC')
        os.mkdir(DCC_dir)
        self.lvm_DCC_filepaths = {}
        for c, ct1 in enumerate(self.unique_cell_types):
            for ct2 in self.unique_cell_types[(c+1):]:
                self.lvm_DCC_filepaths[ct1,ct2] = os.path.join(DCC_dir, 
                                                               '{}_{}.txt'.format(ct1,ct2))
                self.lvm_DCC_filepaths[ct2,ct1] = os.path.join(DCC_dir, 
                                                               '{}_{}.txt'.format(ct1,ct2))
        
        self.jobs_params = [[chrom, pool_id, pool_bands] 
                            for chrom in self.pools_idx
                            for (pool_id, pool_bands) in enumerate(self.pools_idx[chrom])
                            ]
        print('number of jobs: {}.'.format(len(self.jobs_params)))
        # TODO: temporary unavailable resources if n_jobs > 1
        
        self.results = Parallel(n_jobs = 1)(delayed(self.train)
                                                      (chrom, pool_id, pool_bands,
                                                       self.nLatent, self.dataset.batchFlag,
                                                       self.gpuFlag) 
                                                      for (chrom, pool_id, pool_bands) in self.jobs_params)
       
        all_latents = [res[0] for res in self.results if res != None]
        self.concatenated_latent = np.hstack(all_latents)
        full_latent_path = os.path.join(self.latent_dir, 'full_latent')
        np.save(full_latent_path, self.concatenated_latent)
        
        self.store_imputed_cools()
        
    def train(self, chrom, pool_id, bands, nLatent, batchFlag, gpuFlag):
        
        if len(bands) > 1:
            pool_mat = np.hstack([self.dataset.all_band_mats[chrom][b] for b in bands])
            bins_id = np.hstack([self.dataset.get_band_lps(chrom, b) for b in bands])
        else:
            pool_mat = self.dataset.all_band_mats[chrom][bands[0]]
            bins_id = self.dataset.get_band_lps(chrom, bands[0])
            
        orig_num_cells = pool_mat.shape[0]
        valid_cells = list(np.where(pool_mat.sum(axis = 1) > 0)[0])
        if len(valid_cells) == 0:
            return 
        else:
            pool_mat = pool_mat[valid_cells,:]
            orig_num_lps = pool_mat.shape[1]
            # valid locus pairs
            valid_lps = list(np.where(pool_mat.sum(axis=0) > 0)[0])
            pool_mat = pool_mat[:, valid_lps]
            avg_poolDepth = pool_mat.sum(axis = 1).mean()
            adata = sc.AnnData(pool_mat)
            if batchFlag:
                adata.obs['batch'] = self.dataset.batch_names[valid_cells]
                adata.obs['cell_type'] = self.dataset.cell_types[valid_cells]
                scvi.model.SCVI.setup_anndata(adata, batch_key = 'batch')
            else: 
                scvi.model.SCVI.setup_anndata(adata)
                
            model = scvi.model.SCVI(adata, n_latent = nLatent)
            model.train(use_gpu = gpuFlag)
            
            pool_latent = np.zeros((orig_num_cells, nLatent))
            pool_latent[valid_cells, :] = model.get_latent_representation()
            
            
            imputed_pool = np.zeros((orig_num_cells, orig_num_lps))
            if batchFlag:
                for batch_name in list(set(adata.obs['batch'])):
                    sub_imputed_pool = model.get_normalized_expression(
                        library_size = avg_poolDepth,
                        transform_batch = batch_name)
                    imputed_pool[np.ix_(valid_cells, valid_lps)] += sub_imputed_pool
                imputed_pool /= len(self.dataset.batch_types)
            else:
                sub_imputed_pool = model.get_normalized_expression(
                    library_size = avg_poolDepth)
                imputed_pool[np.ix_(valid_cells, valid_lps)] = sub_imputed_pool
            
            imputed_bands = self.pool2bands(imputed_pool, chrom, bands)
            tmp_latent_path = os.path.join(self.latent_dir, '{}_pool{}_latent'.format(chrom, pool_id))
            np.save(tmp_latent_path, pool_latent)
            
            
            bin1_id = bins_id[0, valid_lps]
            bin2_id = bins_id[1, valid_lps]
            for c, ct1 in enumerate(list(set(adata.obs['cell_type']))):
                for ct2 in list(set(adata.obs['cell_type']))[(c+1):]:
                    diff = model.differential_expression(groupby = "cell_type",
                                                         group1 = ct1, group2 = ct2)
                    diff['chrom'] = chrom
                    diff['bin1_id'] = bin1_id 
                    diff['bin2_id'] = bin2_id
                    if not os.path.exists(self.lvm_DCC_filepaths[ct1,ct2]):
                        diff.to_csv(self.lvm_DCC_filepaths[ct1,ct2], sep = "\t", index = False,
                                    mode = 'w')
                    else:
                        diff.to_csv(self.lvm_DCC_filepaths[ct1,ct2], sep = "\t", index = False,
                                    header = None, mode = 'a')
            
            return pool_latent, imputed_bands, model
            
    
    def pool2bands(self, pool_mat, chrom, bands):
        
        band_lengths = [self.dataset.all_band_mats[chrom][b].shape[1] for b in bands]
        assert pool_mat.shape[1] == sum(band_lengths), "pool matrix dimension is different than sum of corresponding band lengths"
        band_ends = list(np.cumsum(band_lengths))
        band_starts = [0] + band_ends[:-1]
        band_mats = [pool_mat[:,s:e] for s, e in zip(band_starts, band_ends)]
        return band_mats
     
    
    def store_imputed_cools(self):
        
        print('start storing imputed cools.')
        imputed_cools_dir = os.path.join(self.out_dir, 'imputed_cools')
        os.mkdir(imputed_cools_dir)
        
        bins, offsets = make_bins(self.dataset.chrom_sizes_filepath,
                                  self.dataset.chromosomes,
                                  self.dataset.resolution)
        
        def store_imputed_cool(c):
            pixels = pd.DataFrame()
            for j, (chrom, pool_id, pool_bands) in enumerate(self.jobs_params):
                num_bins = self.dataset.num_bins[chrom]
                for b, band in enumerate(pool_bands):
                    if not self.results[j] is None:
                        bin1_id = np.arange(0, num_bins - band) + offsets[chrom]
                        bin2_id = np.arange(band, num_bins) + offsets[chrom]
                        pixels_tmp = pd.DataFrame({'bin1_id': bin1_id,
                                                'bin2_id': bin2_id,
                                                'count': self.results[j][1][b][c,:]})
                        pixels_tmp = pixels_tmp[pixels_tmp['count'] != 0]
                        if pixels.empty:
                            pixels = pixels_tmp
                        else:
                            pixels = pd.concat([pixels, pixels_tmp])
            coolname = '{}.cool'.format(os.path.splitext(os.path.basename(self.dataset.filepaths[c]))[0])
            coolpath = os.path.join(imputed_cools_dir, coolname)
            cooler.create_cooler(coolpath,
                        bins = bins,
                        pixels = pixels,
                        dtypes={'count': np.float64})
        Parallel(n_jobs = self.n_cpus)(delayed(store_imputed_cool)(c) 
                                       for c in tqdm(range(self.dataset.cell_num)))
                                                    
    

def get_pool_idx(max_diag):
    
    curr_diag = 1
    pool_idx = [[curr_diag]]
    curr_pool_size = 2
    while (curr_diag + 1) < max_diag:
        pool_idx.append([i for i in np.arange(curr_diag + 1, min(curr_diag + curr_pool_size + 1, max_diag))])
        curr_diag = curr_diag + curr_pool_size
        curr_pool_size += 1
    return pool_idx

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