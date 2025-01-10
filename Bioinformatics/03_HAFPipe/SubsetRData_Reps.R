setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/')
library(tidyverse)
load('./orch2021_Downsampled_META_Filtered.RData', verbose = TRUE)
df = cbind(samps, t(afmat))
df.eec = cbind(samps, t(eec))

breps.df = df %>% ##get list of samples with biological reps
    filter(biol.rep == "Yes")
breps = breps.df$sample    

df = df %>%
    filter(sample %in% breps) %>%
    filter(tech.rep == "No")

df.eec = df.eec %>%
    filter(sample %in% breps) %>%
    filter(tech.rep == "No")

samps = df[,1:ncol(samps)]

df = df[,-c(1:ncol(samps))]
df.eec = df.eec[,-c(1:ncol(samps))]

afmat = t(df)
eec = t(df.eec)


save(samps, sites,afmat, eec, file = "./orch2021_Downsampled_BiolReps.RData" )


##Getting All Cages technical reps RData - Run on cluster in SubsetRData.R
setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/RData/Downsampled/')
library(tidyverse)


load('./orch2021_Downsampled_META_Filtered.RData', verbose = TRUE)

dim(afmat)
dim(eec)
df = cbind(samps, t(afmat))
df.eec = cbind(samps, t(eec))

techreps.df = df %>% ##get list of samples with biological reps
    filter(tech.rep == "Yes")
techreps = techreps.df$sample


df = df %>%
    filter(sample %in% techreps) %>%
    filter(!full.sample.name %in% c("tp11_F1_E2_downsamped", 'tp11_F1_E3_downsamped', 'tp11_F1_E4_downsamped', 'tp11_F1_E6_downsamped','tp11_F1_E6_downsamped', 'tp11_F1_E9_downsamped', "tp3_F1_P7_downsamped" )) 

df.eec = df.eec %>%
    filter(sample %in% techreps) %>%
    filter(!full.sample.name %in% c("tp11_F1_E2_downsamped", 'tp11_F1_E3_downsamped', 'tp11_F1_E4_downsamped', 'tp11_F1_E6_downsamped','tp11_F1_E6_downsamped', 'tp11_F1_E9_downsamped', "tp3_F1_P7_downsamped" )) 



samps = df[,1:ncol(samps)]

df = df[,-c(1:ncol(samps))]
df.eec = df.eec[,-c(1:ncol(samps))]

afmat = t(df)
eec = t(df.eec)

save(samps, sites,afmat, eec, file = "./orch2021_Downsampled_TechReps.RData")