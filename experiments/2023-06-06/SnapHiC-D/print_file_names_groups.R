
for(i in 1:22){
  CHR=sprintf('chr%s',i)
  
  dir.in.A  <- '/project/compbio-lab/scHi-C/SnapHiC-D_example_data/sox2_ESC/'
  dir.in.B  <- '/project/compbio-lab/scHi-C/SnapHiC-D_example_data/sox2_NPC/'
  dir.out <- '/project/compbio-lab/scHi-C/SnapHiC-D_example_data/sox2_outdir/'
  
  myFiles_A=list.files(dir.in.A,pattern =sprintf("\\%s.normalized.rwr.bedpe",CHR))
  myFiles_B=list.files(dir.in.B,pattern =sprintf("\\%s.normalized.rwr.bedpe",CHR))
  
  myFiles_A = cbind(myFiles_A,'A')
  myFiles_B = cbind(myFiles_B,'B')
  output=rbind(myFiles_A,myFiles_B)
  colnames(output)= c('name','group')
  
  write.table(output, file=paste(dir.out,sprintf('file_list.txt',CHR), sep=''),row.names=F, col.names=T, sep=' ', quote=F)
}


#### printing subset of files from an existing file list

# for(i in 1:22){
#   CHR=sprintf('chr%s',i)
#   dir.out <- '/dir/to/file/output/'
# ann <- read.table(sprintf('./file_lists/file_name_%s.txt',CHR), head=T)
# A=ann[which(ann$group=="A"),]
# B=ann[which(ann$group=="B"),]
# output=rbind(A[1:100,], B)
# 
# write.table(output, file=paste(dir.out,sprintf('name_file_100_%s.txt',CHR), sep=''),row.names=F, col.names=T, sep=' ', quote=F)
# 
# }