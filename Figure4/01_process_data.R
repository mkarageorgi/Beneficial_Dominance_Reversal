library(tidyverse)
library(reticulate)
# loads two matrices:
#  - afmat with allele frequencies
#  - eec
# and two dataframes:
#  - samps
#  - sites
load("data/raw/orch2021_Downsampled_META_Filtered.RData")

write_csv(samps, "data/raw/samps.csv")
write_csv(sites, "data/raw/sites.csv")

# np <- import("numpy")
#
# np$save("./afmat.npy", r_to_py(afmat))
# np$save("./eec.npy", r_to_py(eec))

samps <- as_tibble(samps)
sites <- as_tibble(sites)
freqs <- as.matrix(afmat)

sites <- sites %>%
  mutate(site_idx = seq_len(nrow(sites))) # %>%
  # drop rows for chrom X
  # filter(chrom != "X")

samps <- samps %>%
  mutate(freq_idx = seq_len(nrow(samps))) %>%
  # make the columns numeric
  mutate(cage = as.numeric(str_extract(cage, "\\d+")), tpt = as.numeric(tpt)) %>%
  # drop the I treatment
  filter(treatment != "I") %>%
  # remove E11 and E12 to have balanced data between E and P cages
  filter(cage <= 10) %>%
  # drop biol.rep and tech.rep
  filter(biol.rep != "Yes", tech.rep != "Yes") %>%
  # select matching timepoints
  filter(tpt %in% c(1, 3, 5, 7, 9, 10, 11, 12)) %>%
  # relabel timepoints similar to marianna
  mutate(tpt = case_when(
    tpt == 1 ~ 1,
    tpt == 3 ~ 2,
    tpt == 5 ~ 3,
    tpt == 7 ~ 4,
    tpt == 9 ~ 5,
    tpt == 10 ~ 6,
    tpt == 11 ~ 7,
    tpt == 12 ~ 8, TRUE ~ as.numeric(tpt)
  )) %>%
  # sort by treatment, timepoint, cage
  arrange(treatment, tpt, cage)

missing <- samps %>%
  mutate(ttpt = paste0(treatment, "t", as.character(tpt))) %>%
  group_by(ttpt, cage) %>%
  dplyr::summarise(n = n()) %>%
  pivot_wider(names_from = cage, values_from = n) %>%
  arrange(ttpt)

freqs <- freqs[sites$site_idx, samps$freq_idx]

# update site_idx and freq_idx in sites and samps again
sites <- sites %>%
  mutate(site_idx = seq_len(nrow(sites)))
samps <- samps %>%
  mutate(freq_idx = seq_len(nrow(samps)))

write_csv(samps, "data/processed/samps.csv")
write_csv(sites, "data/processed/sites.csv")

np <- import("numpy")
#reticulate::use_virtualenv()
np$save("data/processed/afmat.npy", r_to_py(freqs))
