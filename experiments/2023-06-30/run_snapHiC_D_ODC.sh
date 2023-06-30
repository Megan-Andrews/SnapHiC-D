#!/bin/bash
#SBATCH -J Astro_MG_diff     # Name that will show up in squeue
#SBATCH --gres=gpu:0         # Request 4 GPU "generic resource"
#SBATCH --cpus-per-task=8
#SBATCH --time=0-03:00       # Max job time is 3 hours
#SBATCH --mem=16G            # Max memory (CPU) 16GB
#SBATCH --output=Astro_MG_diff.out   # Terminal output to file named (hostname)-(jobid).out
#SBATCH --partition=long     # long partition (allows up to 7 days runtime)
#SBATCH --nodelist=cs-venus-03   # if needed, set the node you want (similar to -w xyz)


# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate vae


############################################################################
###                            User Variables                            ###
############################################################################

SnapHiC_D_dir="/home/maa160/SnapHiC-D/experiments/2023-06-30/SnapHiC-D"
group_A_dir="/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/ODC_1_100kb_imputed_cool"
group_B_dir="/project/compbio-lab/scHi-C/Lee2019/100kb_imputed_cool/ODC_2_100kb_imputed_cool"
out_dir="/project/compbio-lab/scHi-C/Lee2019/results/2023-06-29/ODC_diff_batched"
chr="chr22"
genome="hg19"
num_CPUs=8
bin_size=100000
fdr_threshold=0.1
max_gap=101
min_gap=2



############################################################################


python ${SnapHiC_D_dir}/snapHiC_diff2.py -s $SnapHiC_D_dir -i $group_A_dir -j $group_B_dir -o $out_dir -c $chr -g $genome -n $num_CPUs --binsize $bin_size --fdr_threshold $fdr_threshold --maxi_gap $max_gap --mini_gap $min_gap
