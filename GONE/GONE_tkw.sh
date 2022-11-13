#!/bin/bash
#SBATCH --job-name= GONE_tkw
#SBATCH --output=GONE_twk.log
#SBATCH --mail-user=martin.kardos@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH --mem=10GB
#SBATCH -t 500:00:00
#SBATCH -D /scratch/mkardos/orca/GONE/
bash script_GONE.sh kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25_tkwOnly_autosomes &

