#!/bin/bash
#SBATCH --job-name=02B_SortBams
#SBATCH --time=40:00:00
#SBATCH --mem=70GB
#SBATCH --cpus-per-task=8
#SBATCH --partition=hns
#SBATCH --out=out.02B_SortBams_3Files

module load biology
module load samtools

Dir=${1}

cd ${Dir}

for i in ./*.filtered.bam;do
    samtools sort -@ 16 $i > ${i/.filtered.bam/.filtered.sorted.bam}
done
