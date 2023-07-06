library('rhdf5')
library(pryr)
library(InteractionSet)

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
   mats <- list()
    for (file in files) {
      print(file)
      tryCatch(
        {
          mat <- cool2mat(file, chr_name = chr_name)
          mats[[length(mats) + 1]] <- mat
        },
        error = function(err) {
          print("HDF5 error")
        }
      )
    }
    # mats = lapply(files, cool2matrix, chr = chr_name)
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

make_multiHiCcompare_object <- function(input_format, files, chr_name, 
                                        chr_size_file, resolution, groups){
  if (input_format == 'txt'){
    dfs = lapply(files, function(f){
      df = read.table(f)
      colnames(df) = c('bin1_id', 'bin2_id', 'count')
      df$chr = chr_name
      df = df[,c('chr', 'bin1_id', 'bin2_id', 'count')]
      return (df)
    })
  }
  else if (input_format == 'cool'){
    dfs = lapply(files, cool2sparse, chr = chr_name, resolution = resolution)
  }
  else{
    print('txt or cool formats are available.')
    return 
  }
  hicexp = make_hicexp(data_list=dfs, groups = groups, A.min = 5, zero.p = 1)
  return(hicexp)
}