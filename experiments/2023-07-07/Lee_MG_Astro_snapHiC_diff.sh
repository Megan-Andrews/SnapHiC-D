#!/bin/bash
#SBATCH -J snapHiC_diff
#SBATCH --gres=gpu:0
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 11
#SBATCH --output=snapHiC_diff.out
#SBATCH --partition=long
#SBATCH --nodelist=cs-venus-03


source ~/miniconda3/etc/profile.d/conda.sh
conda activate snapHiC

mpirun -n 11  python ../../snapHiC_diff_mpi.py  --indir-A /project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_rwr_cool/Astro --indir-B /project/compbio-lab/scHi-C/Lee2019/Human_single_cell_10kb_rwr_cool/MG -r 10000 -o /project/compbio-lab/scHi-C/Lee2019/results/2023-07-10/snapHiC_diff_MG_Astro_10kb --min-distance 10000 --max-distance 2000000 --snapHiC-dir ~/projects/SnapHiC-D -g hg19
