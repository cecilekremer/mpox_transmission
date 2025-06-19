#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch_sapphirerapids
#SBATCH --account=lp_h_simpact
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96

module load R/4.1.0-foss-2021a 

Rscript $VSC_SCRATCH/simNet/03a_network_vsc.R 
