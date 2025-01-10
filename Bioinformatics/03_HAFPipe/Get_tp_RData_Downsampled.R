#Doing each time point separately

library(data.table)
library(tidyverse)



tpts = c(1:11, 13)
for (tpt in tpts){
    setwd("/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/04_HAFs_V2/")
    data = data.frame()
    files <- list.files(pattern=paste0("tp",tpt,"_"), full.names=TRUE, recursive=FALSE)
    for (file in files){
    #Step 1 - eactract full sample name from each HAFpipe output file 
    clean.file <- strsplit(file, "./")[[1]][2]
    full.sample.name <- strsplit(clean.file, ".bam")[[1]][1]
    chr <- strsplit(clean.file, ".afSite")[[1]][1] #extract chrom info
    chr <- strsplit(chr, ".bam.")[[1]][2]
    freqs <- read.csv(file, header = TRUE)
    freqs = freqs %>%
        mutate(sample.name = full.sample.name) %>%
        mutate(chrom = chr)
    freqs = freqs %>%
        select(sample.name, chrom, pos, af)
    data = rbind(data, freqs)        
    }
    sp.data = data %>%
        spread(sample.name, af)

#####sp.data = na.omit(sp.data)#####not doing this until after treatments segregated

    sites = as.data.frame(sp.data %>%
    dplyr::select(chrom, pos))
    colnames(sites) = c('chrom', 'pos')

    samps = as.data.frame(colnames(sp.data)[3:ncol(sp.data)])
    colnames(samps) = c('sample.name')

    afmat = as.matrix(sp.data[,3:ncol(sp.data)])

    save(sites, samps, afmat, file = paste0("../RData_V2/orch2021_Downsampled_tp",tpt,".RData"))
}
