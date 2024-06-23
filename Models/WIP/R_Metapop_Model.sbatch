#!/bin/bash

#SBATCH --job-name=R-metapop-model
#SBATCH --account=pi-cbrook
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4G

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."

module load R

Rscript ./Scripts/6_ConnectMat_pt2_WIP.R
