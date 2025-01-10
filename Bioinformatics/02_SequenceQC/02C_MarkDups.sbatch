#!/bin/bash
#SBATCH --job-name=02C_MarkDups
#SBATCH --time=60:00:00
#SBATCH --mem=80GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=hns
#SBATCH --out=out.02C_MarkDups

ml biology
ml java

Dir=${1}

cd ${Dir}


for i in ./*.filtered.sorted.bam;do
    java -jar /home/users/mcbitter/dpetrov/SOFTWARE/picard.jar MarkDuplicates \
	REMOVE_DUPLICATES=TRUE \
	I=./$i \
        O=./${i/.bam/.deduped.bam} \
        M=/scratch/groups/dpetrov/MarkB/SequencingPractice/${i/.bam/marked_dup_metrics.txt}

done
