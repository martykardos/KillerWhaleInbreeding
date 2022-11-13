#!/bin/sh
#SBATCH --job-name=r_interleave
#SBATCH --output=r_interleave.log
#SBATCH --mem=90GB
#SBATCH -t 2000:00:00
#SBATCH -c 10
#SBATCH -D /scratch/mkardos/orca/cetaceanGenomes/cactus
module load R
Rscript rCode_interleave.R
