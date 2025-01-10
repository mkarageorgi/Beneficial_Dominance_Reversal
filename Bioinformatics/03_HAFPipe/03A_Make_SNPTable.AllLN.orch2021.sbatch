#!/bin/bash                                                                                                          
#SBATCH --job-name=make_snptable                                                                                     
#SBATCH --time=08:00:00                                                                                              
#SBATCH --mem=40G                                                                                                    
#SBATCH --cpus-per-task=1                                                                                            
#SBATCH --partition=hns                                                                                              
#SBATCH --out=make_snptable.log                                                                                      

cd /home/users/mcbitter/dpetrov/MarkB/

ml R
ml py-numpy/1.14.3_py27
ml biology
ml samtools

chroms=(2L 2R 3L 3R X)

for chrom in ${chroms[*]};do \
    /home/users/mcbitter/OrchardProject/Code/HAFpipe-line/HAFpipe_wrapper.sh -t 1,2 \
        -l HAF-pipe.make.snptable.log \
        -v ./Orchard2021Data/founders/vcfs/inbredv2.filtered.Orchard2021.vcf \
        -c ${chrom} \
        -s ./Orchard2021Data/founders/inbredv2_withHets.orch2021.${chrom}.snpTable \
        -k --keephets \
        --mincalls 40 \
        -i npute

done
