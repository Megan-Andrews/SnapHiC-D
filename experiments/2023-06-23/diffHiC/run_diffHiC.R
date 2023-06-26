library(diffHic)
library(rhdf5)
library(GenomicRanges)
library(InteractionSet)
source("cool2matrix.R")

input_directory <- "/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/"  #should we run diffHiC on 10kb data 100kb data or the imputed data? 
astro_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-23/diffHiC/file_lists/Astro_100kb_file_list.txt"  
mg_file_list_path <- "/home/maa160/SnapHiC-D/experiments/2023-06-23/diffHiC/file_lists/MG_100kb_file_list.txt" 
astro_file_list <- readLines(astro_file_list_path)
mg_file_list <- readLines(mg_file_list_path)

# m <- cool2matrix("/project/compbio-lab/scHi-C/Lee2019/Human_single_cell_100kb_cool/181218_21yr_2_C1_AD002_Astro_100kb_contacts.cool
# ")
# print(dim(m)[1])

astro_count <- 0
mg_count <- 0

matrix_size <- 514
chromosome <- rep("chr22", matrix_size)
start <- c(0:(matrix_size-1))
end <- c(1:matrix_size)
regions <- GRanges( chromosome, IRanges(start, end))

print(paste0("matrix size ", matrix_size))
print(paste0("region length ", length(regions)))
# print(paste0("regions  ", regions))

print("creating matrices")
true_matrix <- as.data.frame(matrix(TRUE, nrow = matrix_size, ncol = matrix_size))
true_CM <- ContactMatrix(true_matrix, regions, regions)
false_matrix <- matrix(FALSE, nrow = matrix_size, ncol = matrix_size)
empty_matrix <- as.data.frame(matrix(0, nrow = matrix_size, ncol = matrix_size))
empty_CM <- ContactMatrix(empty_matrix, regions, regions)
# print(empty_matrix)
astro_data <- NULL #ContactMatrix(empty_matrix, regions, regions)
mg_data <- NULL #ContactMatrix(empty_matrix, regions, regions)


procces_into_matrix <- function(file_name) {
  file_path <- paste0(input_directory, file_name)
  matrix <- cool2matrix(file_path, chr="chr22")
  cm <- ContactMatrix(matrix, regions, regions)
  return(cm)
}

print("processing astro matrices")
cm1 <- procces_into_matrix(astro_file_list[1])
cm2 <- procces_into_matrix(astro_file_list[2])

to.keep <- as.matrix(cm1)!=0 | as.matrix(cm2)!=0
print(c("length", length(cm1), "dim", dim(cm1)))
print(c("length", length(cm2), "dim", dim(cm2)))
print(c("length", length(to.keep), "dim", dim(to.keep)))

iset1 <- deflate(cm1, extract=to.keep)
iset2 <- deflate(cm2, extract=to.keep)
data <- cbind(iset1, iset2)
# for (astro_file in astro_file_list) {
#   matrix <- procces_into_matrix(astro_file)
#   if (is.null(astro_data)) {
#     astro_data <- matrix
#   } else {
#     # print(paste("matrix result dim:", dim(matrix), "astro data:", dim(astro_data)))
#     to.keep <- false_matrix
#     to.keep[as.matrix(astro_data) != 0 | as.matrix(matrix) ] <- TRUE
#     print(to.keep)
#     print(as.matrix(astro_data))
#     print(as.matrix(matrix))
#     print(c("length", length(astro_data), "dim", dim(astro_data)))
#     print(c("length",length(matrix), "dim", dim(matrix)))
#     print(c("length",length(true_matrix), "dim", dim(true_matrix)))

#     print("deflate data")
#     iset1 <- deflate(astro_data, extract=to.keep)
#     print("deflate new matrix")
#     iset2 <- deflate(matrix, extract=to.keep)
#     astro_data <- cbind(iset1, iset2)
#     # print(astro_data)
#     # astro_data <- cbind(astro_data, matrix)
#   }
#   astro_count <- astro_count + 1
#   print(paste0("astro count: ", astro_count))
#   if(astro_count==30){
#     break
#   }
# }
# print(astro_data)

# print("processing mg matrices")
# for (mg_file in mg_file_list) {
#   matrix_result <- procces_into_matrix(mg_file)  
#   mg_data <- mergeCMs(matrix_result, mg_data)
#   mg_count <- mg_count + 1
#   print(paste0("astro count: ",mg_count))
# }
# print(mg_data)

# print("Create a diffHicSet object")
# hicset <- new("diffHicSet", matrices = list(astro_matrices, mg_matrices))

# print("Preprocess and normalize the data")
# hicset <- filterEmptyBins(hicset)
# hicset <- filterByRPM(hicset)
# hicset <- normalize(hicset, method = "loess", span = 0.3)

# print("Perform differential analysis")
# hicset <- testInteractions(hicset, method = "binom")

# print("Extract significant interactions") 
# significant_interactions <- getSignificantInteractions(hicset)

# print("Adjust p-values for multiple testing") 
# significant_interactions <- adjustSignificance(hicset, fdr.method = "BH")

# print("Extract differentially interacting regions")
# dcc_regions <- getDCC(significant_interactions)

# print("Print the results")
# print(dcc_regions)