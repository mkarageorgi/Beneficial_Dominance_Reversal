#!/bin/bash
#SBATCH --job-name=02D_IndexBams
#SBATCH --time=08:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=hns
#SBATCH --out=out.02D_IndexBams

module load biology
module load samtools

#cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/
Dir=${1}

cd ${Dir}

for i in ./*.bam;do
    samtools index -b $i
done
