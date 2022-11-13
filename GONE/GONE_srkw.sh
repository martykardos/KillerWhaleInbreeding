#!/bin/bash
#SBATCH --job-name= GONE_srkw
#SBATCH --output=GONE_srwk.log
#SBATCH --mail-user=martin.kardos@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH --mem=10GB
#SBATCH -t 500:00:00
#SBATCH -D /scratch/mkardos/orca/GONE/
bash script_GONE.sh kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1.hetDepthFiltered.maxMiss0.25.minGQ30.srkw.maxMiss25Perc.mac1.filtered

