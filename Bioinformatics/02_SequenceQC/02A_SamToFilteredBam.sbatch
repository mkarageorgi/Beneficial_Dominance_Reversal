#!/bin/bash
#SBATCH --job-name=02A.SamToBam
#SBATCH --time=72:00:00
#SBATCH --mem=80GB
#SBATCH --cpus-per-task=12
#SBATCH --partition=hns
#SBATCH --out=out.02A.SamToBam

module load biology
module load samtools
module load bamtools

Dir=${1}

cd ${Dir}

for i in ./*.sam;do
    samtools view -@ 24 -S -bq 20 $i > ${i/.sam/.filtered.bam}
done
