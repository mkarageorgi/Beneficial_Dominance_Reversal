library(tidyverse)
library(purrr)
library(data.table)
chr2L<- fread("./glm_malathion_2L.csv")
chr2R<- fread("./glm_malathion_2R.csv")
chr3L<- fread("./glm_malathion_3L.csv")
chr3R<- fread("./glm_malathion_3R.csv")
chrX<- fread("./glm_malathion_X.csv")

chromall<- rbind(chr2L,chr2R,chr3L,chr3R,chrX) %>% rename(pos = column)

  sites.cage_comparison = chromall%>%
  
  separate_wider_delim(treat_range, delim = ".", names = c("cage", "comparison")) %>%

  #filter tpt term - otherwise FDR will also account the p.values from the intercept
  filter(term =="tpt") %>% 
  
  # group by treat_range to get all chromosomes for each time range for E or P cages
  group_by(cage,comparison) %>% 
  
  #fdr correction for all chromosomes per treat_range
  mutate(p.value.adjusted = p.adjust(p.value, method="BH"), .after = p.value) %>%
  
  # score SNPs based on their significance - be aware that the ordering matters
  # condition are ordered from the most to the least stringest
  mutate(sigLevel = case_when(
    p.value.adjusted < 0.01 & abs(effect_size) > 0.02 ~ 3,  
    p.value.adjusted < 0.05 & abs(effect_size) > 0.02 ~ 2,  
    p.value.adjusted < 0.2  ~ 1,                           
    TRUE ~ 0
    ), .after = p.value.adjusted) %>%
  # convert sigLevel to a factor
  mutate(sigLevel = factor(sigLevel, levels = c(0, 1, 2, 3))) %>%
  ungroup()

  write_csv(sites.cage_comparison, "../data/raw/sigsite_malathion.csv")
