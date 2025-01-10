#!/bin/bash
#SBATCH --job-name=02E_ComputeDepth
#SBATCH --time=16:00:00
#SBATCH --mem=300GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.02E.ComputeDepth

module load biology
module load samtools

#cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/
Dir=${1}

cd ${Dir}

for i in ./*.bam;do
    samtools depth -a $i > ${i/.bam/.dp.txt}
done

