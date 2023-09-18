# read in individual coolers, aggregate them and convert them into a file with the format

# evaluatig the QC score or reproducibility of the coolers? which method? all?

import hicrep
import hicrep.utils
import cooler 

cell1 = 'GM12878' 
cell2 = 'H1Esc' #from Kim2020

Kim2020_cools = "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool"
rwr_cools = "/project/compbio-lab/scHi-C/Kim2020/kim2020_cool_rwr3"
scVI_cools = "/project/compbio-lab/scHi-C/Kim2020/results/2023-09-14/Kim2020_scVI_model/imputed_cools"
higashi_k0_cools = "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool_higashi_k0"
higashi_k5_cools =  "/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool_higashi_k5"





