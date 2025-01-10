#!/bin/bash
#SBATCH --job-name=01.orch.2021.alignment_RD1
#SBATCH --time=48:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=10
#SBATCH --partition=hns
#SBATCH --out=out.01.orch.2021.alignment_RD1

module load biology
module load bwa

Dir=${1}

cd ${Dir}
###file containing all sequencing files
awk 'NR>0' orch2021.samps.csv > orch2021.samps.txt

#loop througha ll sequencing files and align to reference genome (v. 5) using bwa
while read samp
do
    bwa mem -t 20 ../../../../Dmel_RefGenome/all_dmel.fasta "$samp"_R1_001.fastq.gz "$samp"_R2_001.fastq.gz > "$samp".sam
done < orch2021.samps.txt

