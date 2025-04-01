#!/bin/bash -l

#SBATCH --clusters=wice
#SBATCH --partition=batch
#SBATCH --account=lp_h_simpact
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load R/4.1.0-foss-2021a 

Rscript $VSC_SCRATCH/simNet/03c_serial_interval_reg_age_hh_VSC.R 
