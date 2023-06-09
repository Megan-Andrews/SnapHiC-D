library('rhdf5')
library(pryr)
library(purrr)
library(InteractionSet)
library(dplyr)
library(data.table)
library(Matrix)
library(sparsesvd)
library(umap)
library(csaw)
library(multiHiCcompare)
#TODO: add txt loader

### this function is from https://github.com/TaoYang-dev/hicrep/blob/master/R/cool2matrix.R
cool2matrix <- function(file, chr = 'chr1') {

  pixels <- h5read(file, c('pixels'));H5close()
  bins <- h5read(file, c('bins'));H5close()
  chrom.offset <- h5read(file, 'indexes/chrom_offset');H5close()

  chr2index <- function(chr) {
    if (substr(chr, 1, 3) == 'chr') {
      chr = substr(chr, 4, nchar(chr))
    }
    if (chr == 'X') {
      index = 23
    }
    else if (chr == 'Y') {
      index = 24
    }
    else if (chr == 'M') {
      index = 25
    }
    else {
      index = as.integer(chr)
    }
    return(index)
  }

  chrom.index <- chr2index(chr)
  chrom.range <- chrom.offset[chrom.index:(chrom.index+1)] + c(0, -1) 
  n.rows <- chrom.range[2] - chrom.range[1] + 1 
  bin.range <- which(pixels$bin1_id >= chrom.range[1] & 
                       pixels$bin1_id <= chrom.range[2] & 
                       pixels$bin2_id >= chrom.range[1] & 
                       pixels$bin2_id <= chrom.range[2])
  n.bins <- length(bin.range) 
  mat <- matrix(0, ncol = n.rows, nrow = n.rows)
  for (i in 1:n.bins) {
    mat[pixels$bin1_id[bin.range[i]] - chrom.range[1] + 1, 
        pixels$bin2_id[bin.range[i]] - chrom.range[1] + 1] <- 
      pixels$count[bin.range[i]]
    mat[pixels$bin2_id[bin.range[i]] - chrom.range[1] + 1, 
        pixels$bin1_id[bin.range[i]] - chrom.range[1] + 1] <- 
      pixels$count[bin.range[i]]
  }
  chrom.ranges <- (chrom.range[1]+1):(chrom.range[2]+1)
  # mat <- as.data.frame(mat)
  
  return(mat)
}


cool2sparse <- function(file, chr = 'chr1', resolution) {
  pixels <- h5read(file, c('pixels'));H5close()
  bins <- h5read(file, c('bins'));H5close()
  chrom.offset <- h5read(file, 'indexes/chrom_offset');H5close()
  chr2index <- function(chr) {
    if (substr(chr, 1, 3) == 'chr') {
      chr = substr(chr, 4, nchar(chr))
    }
    if (chr == 'X') {
      index = 23
    }
    else if (chr == 'Y') {
      index = 24
    }
    else if (chr == 'M') {
      index = 25
    }
    else {
      index = as.integer(chr)
    }
    return(index)
  }

  chrom.index <- chr2index(chr)
  chrom.range <- chrom.offset[chrom.index:(chrom.index+1)] + c(0, -1) 
  df = data.frame(bin1_id = pixels$bin1_id,
                  bin2_id = pixels$bin2_id, 
                  count = pixels$count)
  df = df[(df$bin1_id >= chrom.range[1]) & 
            (df$bin1_id <= chrom.range[2]) & 
            (df$bin2_id >= chrom.range[1]) & 
            (df$bin2_id <= chrom.range[2]),]
  df[,c('bin1_id','bin2_id')] = (df[,c('bin1_id','bin2_id')]-chrom.range[1])*resolution
  df$chr = chr
  df = df[,c('chr','bin1_id','bin2_id','count')]
  return(df)
}

### convert 3 columns data frame to contact map
df2mat <- function(df, chr_size, resolution){
  colnames(df) = c('bin1_id', 'bin2_id', 'count')
  chr_size = ceiling(chr_size/resolution)
  mat <- matrix(0, nrow=chr_size, ncol=chr_size)
  df[,c('bin1_id','bin2_id')] = df[,c('bin1_id','bin2_id')] / resolution
  mat[as.matrix(df[,c('bin1_id','bin2_id')]+1)] <- df[,'count']
  mat[as.matrix(df[,c('bin2_id','bin1_id')]+1)] <- df[,'count']
  mat
}

make_diffhic_object <- function(input_format, files, chr_name, 
                                chr_size_file, resolution){
  g = read.table(chr_size_file)
  chr_size = as.numeric(g[g[,1]==chr_name,2])
  if (input_format == 'txt'){
    dfs = lapply(files, read.table)
    mats = lapply(dfs, df2mat, chr_size = chr_size, resolution = resolution)
  }
  else if (input_format == 'cool'){
    mats = lapply(files, cool2matrix, chr = chr_name)
  }
  else{
    print('txt or cool formats are available.')
    return 
  }
  print(object_size(mats))
  chr_size = ceiling(chr_size/resolution)
  regions = GRanges(rep(chr_name,chr_size), IRanges(c(0:(chr_size-1)),c(1:chr_size)))
  cms = lapply(mats, ContactMatrix, anchor1 = c(1:chr_size),
                     anchor2 = c(1:chr_size), regions = regions)
  print(object_size(cms))
  to.keep = Reduce("|", lapply(cms, function(cm){as.matrix(cm)!=0}))
  isets = lapply(cms, deflate, extract = to.keep)
  data = Reduce(cbind, isets)
  interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
  assayNames(data) = 'counts'
  return(data)
}

# make_multiHiCcompare_object <- function(input_format, files, chr_name, 
#                                         chr_size_file, resolution, groups){
#   if (input_format == 'txt'){
#     dfs = lapply(files, function(f){
#       df = read.table(f)
#       colnames(df) = c('bin1_id', 'bin2_id', 'count')
#       df$chr = chr_name
#       df = df[,c('chr', 'bin1_id', 'bin2_id', 'count')]
#       return (df)
#     })
#   }
#   else if (input_format == 'cool'){
#     dfs = lapply(files, cool2sparse, chr = chr_name, resolution = resolution)
#   }
#   else{
#     print('txt or cool formats are available.')
#     return 
#   }
#   hicexp = make_hicexp(data_list=dfs, groups = groups, A.min = 5, zero.p = 1)
#   return(hicexp)
# }

bin_hic <- function(df, resolution){
  df$bin1_id = as.integer(df$bin1_id/resolution)
  df$bin2_id = as.integer(df$bin2_id/resolution)
  min_bin_id = pmin(df$bin1_id, df$bin2_id)
  max_bin_id = pmax(df$bin1_id, df$bin2_id)
  df$bin1_id = min_bin_id
  df$bin2_id = max_bin_id
  d_threshold = max(1000000/resolution,10)
  df = df[df$bin2_id-df$bin1_id <= d_threshold,]
  binned_df <- setDT(df)[,list(count=.N),names(df)]
  binned_df
}

make_regions_vecs <- function(regions_df, include, chrom_size_filepath, 
                                chromosomes, resolution){
  chrom_sizes = read.table(chrom_size_filepath)
  bool_vecs = list()
  for (chrom in chromosomes){
    region_ids = regions_df[regions_df$chr == chrom, 'bin_id']
    region_ids = as(region_ids, 'vector')$bin_id
    chrom_size = ceiling(chrom_sizes[chrom_sizes$V1 == chrom, 'V2']/resolution)
    bool_vecs[[chrom]] = vector(mode = "logical", length = chrom_size)
    bool_vecs[[chrom]][1:chrom_size] = !include
    bool_vecs[[chrom]][region_ids] = include
  }
  bool_vecs
}

### TODO: faster filtering? filtering before or after normalization?

read_hic_df <- function(path, chrom_columns, pos_columns, resolution,
                        filter_regions_path, TSS_regions_path, chrom_size_filepath){
  df = data.frame(fread(path, sep = "\t"))
  df = df[df[,chrom_columns[1]]==df[,chrom_columns[2]],]
  df = df[df[,chrom_columns[1]] %in% paste0('chr', c(1:22)),]
  df = df[,c(chrom_columns[1],pos_columns)]
  colnames(df) = c('chr', 'bin1_id', 'bin2_id')
  binned_df = bin_hic(df,resolution)
  
  # only including valid bins 
  # filter_regions = fread(filter_regions_path)
  # filter_regions = filter_regions[,c(1:2)]
  # colnames(filter_regions) = c('chr','start')
  # filter_regions$bin_id = as.integer(filter_regions$start/resolution)
  # filter_vecs = make_regions_vecs(filter_regions, FALSE, chrom_size_filepath, 
  #                                 paste0('chr', c(1:22)), resolution)
  # is_valid = function(x){filter_vecs[[x['chr']]][as.integer(x['bin1_id'])] & 
  #     filter_vecs[[x['chr']]][as.integer(x['bin2_id'])]}
  # valid = apply(binned_df, 1, is_valid)
  # valid = as.numeric(valid)
  # valid[is.na(valid)] = 0
  # 
  # # only including TSS regions 
  # TSS_regions = fread(TSS_regions_path)
  # TSS_regions$bin_id = as.integer((TSS_regions$start+TSS_regions$end)/(2*resolution))
  # TSS_vecs = make_regions_vecs(TSS_regions, TRUE, chrom_size_filepath, 
  #                                 paste0('chr', c(1:22)), resolution)
  # #return(list('TSS_vecs'=TSS_vecs,'df'=binned_df))
  # is_TSS = function(x){TSS_vecs[[x['chr']]][as.integer(x['bin1_id'])] |
  #     TSS_vecs[[x['chr']]][as.integer(x['bin2_id'])]}
  # TSS = apply(binned_df, 1, is_TSS)
  # TSS = as.numeric(TSS)
  # TSS[is.na(TSS)] = 0
  # binned_df = binned_df[valid&TSS,]
  binned_df
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# output: one hic dataframe including all cells 

rbind_hic_dfs <- function(hic_dfs){
  all_hic_df = lapply(c(1:length(hic_dfs)), function(id){
    df = hic_dfs[[id]]
    df = df[,c('chr','bin1_id','bin2_id','count')]
    df$cell = paste0('cell',id)
    df
  })
  all_hic_df = rbindlist(all_hic_df)
  all_hic_df
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# input: resolution
# input: vector of cell types corresponding to cells
# output: hicexp object for multiHiCcompare
# TODO: zero.p is used for single cell resolution, A.min for pseudobulk?

make_hicexp_obj <- function(hic_dfs, res, groups){
  hicexp_dfs = lapply(hic_dfs, function(df){
    df = df[df$bin1_id!=df$bin2_id,]
    df[,c(2,3)] = df[,c(2,3)]*res
    data.frame(df)
  })
  hicexp = make_hicexp(data_list=hicexp_dfs, groups = groups, A.min = 0, zero.p = 0.9)
}

# input: one hic dataframe including chr, bin1_id, bin2_id, count, and cell columns
# output: one hic dataframe with a shape of original one including BandNorm counts

bandnorm <- function(hic_df){
  hic_df$diag = hic_df$bin2_id - hic_df$bin1_id
  hic_df = hic_df[hic_df$diag != 0, ]
  band_info <- hic_df %>% group_by(diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(diag) %>% summarise(depth = mean(band_depth))
  hic_df <- hic_df %>% left_join(alpha_j, by = c("diag")) %>% 
    left_join(band_info, by = c("diag", "cell")) %>% 
    mutate(BandNorm = count/band_depth * depth) %>% 
    select(-c(band_depth, depth, count))
  #colnames(hic_df) = c('bin1_id', 'bin2_id', 'cell', 'diag', 'count')
  hic_df
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# output: one InteractionSet object 

make_iset <- function(dfs){
  dfs = lapply(c(1:length(dfs)), function(data_frames, id){
    df = data_frames[[id]]
    colnames(df) = c('chr', 'bin1_id', 'bin2_id', paste0('count',id))
    data.table(df, key = c('chr', 'bin1_id', 'bin2_id'))
  }, data_frames = dfs)
  merged_df =  Reduce(function(dt1,dt2){merge.data.table(dt1,dt2,all.x = TRUE)}, 
                      dfs)
  merged_df[is.na(merged_df)] = 0
  source_gr = GRanges(seqnames = merged_df$chr, 
                      ranges = IRanges(merged_df$bin1_id, end = merged_df$bin1_id+1))
  target_gr = GRanges(seqnames = merged_df$chr,
                      ranges = IRanges(merged_df$bin2_id, end = merged_df$bin2_id+1))
  GI = GInteractions(source_gr, target_gr)
  iset = InteractionSet(as(merged_df[,4:ncol(merged_df)], 'matrix'), GI)
  iset
}

### TODO: coarsen_iset?

# input: one hic dataframe including all cells with chr, bin1_id, bin2_id, \
# count, and cell columns or one hic matrix of size #cells \times #IFs
# input: format is df or mat
# input: vector of cell types corresponding to cells
# output: dataframe with a cell per row and columns correponding to its UMAP embeddings and cell type 

hic_embedding <- function(hic_data, format, groups){
  if (format == 'mat'){
    hic_mat = as(hic_data, 'sparseMatrix')
    hic_mat = hic_mat[,colSums(hic_mat!=0)>20]
  }
  if (format == 'df'){
    hic_df = hic_data
    hic_df$diag = hic_df$bin2_id - hic_df$bin1_id
    hic_df = hic_df[hic_df$diag != 0, ]
    names = unique(hic_df$cell)
    hic_df$cell = as.numeric(factor(hic_df$cell, level = names))
    hic_df = hic_df %>% 
      mutate(featureIndex = paste(bin1_id, bin2_id, sep = "_")) %>%
      select(-c(bin1_id, bin2_id, diag))
    hic_df$featureIndex = as.numeric(factor(hic_df$featureIndex))
    hic_mat = sparseMatrix(i = hic_df$cell, j = hic_df$featureIndex, x = hic_df$count, 
                           dims = c(max(hic_df$cell), max(hic_df$featureIndex)),
                           index1 = TRUE)
    hic_mat = hic_mat[,colSums(hic_mat!=0)>20]
  }
  
  pca_mat = sparsesvd(hic_mat, 20)
  pca_mat = pca_mat$u %*% diag(pca_mat$d)
  
  embedding = umap(pca_mat)$layout
  embedding = data.frame(embedding)
  colnames(embedding) = c("X1", "X2")
  embedding$cell_type = groups
  embedding
}


