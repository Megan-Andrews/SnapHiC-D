#!/bin/bash
#SBATCH -J Lee_MG_job
#SBATCH --gres=gpu:0
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 20
#SBATCH --output=Lee_MG.out
#SBATCH --partition=long
#SBATCH --nodelist=cs-venus-03


source ~/miniconda3/etc/profile.d/conda.sh
conda activate snapHiC

mpirun -n 10  python ../../snapHiC_preprocessing_mpi.py  -f ../../ext/Lee_MG_samples.txt -r 10000 -l ../../ext/hg19.chrom.sizes.sub --format tab-sep --window-size 10000000 --rp 0.05 -o /project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_rwr_cool --chrom-columns 1 3 --pos-columns 2 4 -o temp --extension .txt.gz --upper-distance 2000000 --log-dir log/Lee_MG_log
