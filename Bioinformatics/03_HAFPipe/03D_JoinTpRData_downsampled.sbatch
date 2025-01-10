#!/bin/bash
#SBATCH --job-name=JoinTpRData_downsampled
#SBATCH --time=08:00:00
#SBATCH --mem=70GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=hns
#SBATCH --out=out.JoinTpRData_downsampled 

ml R

Rscript JoinTpRData_downsampled.R
