#!/bin/bash
#SBATCH --job-name= GONE_arkw
#SBATCH --output=GONE_arwk.log
#SBATCH -c 20
#SBATCH --mem=10GB
#SBATCH -t 72:00:00
#SBATCH -D /scratch/mkardos/orca/GONE/
exec script_GONE_2.sh kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1.hetDepthFiltered.maxMiss0.25.minGQ30.arkw.maxMiss25Perc.mac1.filtered

