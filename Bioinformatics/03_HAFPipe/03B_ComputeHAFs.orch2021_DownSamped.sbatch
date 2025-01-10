#!/bin/bash
#SBATCH --job-name=ComputeHAFs
#SBATCH --time=50:00:00
#SBATCH --mem=70GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=hns
#SBATCH --out=out.ComputeHAFs


cd /home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/03_bams/DownSamp_F1_bams/

module load R
ml py-numpy/1.14.3_py27
ml biology
ml samtools
ml system
ml curl


Tpt=${1}
Gens=${2}


for Chrom in 2L 2R 3L 3R X; do
    for i in ./tp${Tpt}_*.bam;do
	/home/users/mcbitter/OrchardProject/Code/HAFpipe-line/HAFpipe_wrapper.sh -t 3,4 \
	    -c ${Chrom} \
	    -s ../../founders/snptables/Orchard2021/inbredv2_withHets.orch2021.${Chrom}.snpTable \
	    -i npute \
	    -b $i \
	    -r ../../../Dmel_RefGenome/all_dmel.fasta \
	    -g ${Gens} \
	    -o ../../04_HAFs_V2/
    done
done                                                                                                                                                                                     
