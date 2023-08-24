#!/bin/bash
#SBATCH -J Kim_rwr_job
#SBATCH --gres=gpu:0
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=long
#SBATCH --output=Kim_rwr.out
#SBATCH --nodelist=cs-venus-03

############################################################################
###                            User Variables                            ###
############################################################################


source ~/miniconda3/etc/profile.d/conda.sh
conda activate SnapEnv

SnapHiC_D_dir="/home/maa160/SnapHiC-D"
file_list="/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool/file_list.txt"
input_dir="/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool"
output_dir="/project/compbio-lab/scHi-C/Kim2020/Kim2020_cool_rwr"
res=500000
chr_lens="/home/maa160/SnapHiC-D/ext/hg19.chrom.sizes.ordered"
format="cooler"
window_size=10000000
rp=0.05
extension=".cool"
upper_distance=40000000

############################################################################


mpirun -n 10 python ${SnapHiC_D_dir}/python/snapHiC_preprocessing_mpi.py -f $file_list --indir $input_dir -o $output_dir -r $res -l $chr_lens --format $format --window-size $window_size --rp $rp --extension $extension --upper-distance $upper_distance
