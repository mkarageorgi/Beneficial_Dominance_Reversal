#!/bin/bash
#SBATCH --job-name=Get_tp_RData
#SBATCH --time=10:00:00
#SBATCH --mem=70GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=hns
#SBATCH --out=out.Get_tp_RData

ml R

Rscript Get_tp_RData_Downsampled.R
