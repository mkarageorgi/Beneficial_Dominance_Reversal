##joining TP .RData
library(tidyverse)
library(data.table)

setwd("/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/RData_V2/")
files <- list.files(pattern="*_tp*", full.names=TRUE, recursive=FALSE)


chrom = c("2L", "2R", "3L", "3R", "X")
pos = NA
af.meta = as.data.frame(cbind(chrom, pos))
af.meta$pos = as.integer(af.meta$pos)
samps.meta = data.frame()

for (file in files){
    load(file)
    df.samps = samps
    df.af = cbind(sites, afmat)
    af.meta = right_join(af.meta, df.af)
    samps.meta = rbind(samps.meta, df.samps)
}


afmat = na.omit(af.meta)
sites = afmat[,c(1:2)]
afmat = afmat[,-c(1:2)]
samps = samps.meta

save(afmat, sites, samps, file = "./orch2021_Downsampled_META_RAW_V2.RData")