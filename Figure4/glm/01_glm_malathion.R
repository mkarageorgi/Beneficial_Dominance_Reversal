###########################
### Code glm - all cages ##########
#The script subsets the data based on a specific chromosome provided as an argument, formats the data, 
#fits quasi-binomial generalized linear models (GLMs) with a formula count~time for different treatment and time ranges
#applies the repeat masker,
#and stores the results in a CSV file. 
###########################

### SETUP - Load libraries
##########################
library(tidyverse)  # data manipulation and visualization
library(future)     #for parallel computation
library(purrr)      # functional programming tools for working with data
library(furrr)      # functional parallel programming
library(broom)      # tidy and summarize model output
library(parallel)   # embarrasingly parallel processing
library(biomartr)   # tool designed for genomic data retrieval and functional annotation (incl. repeat masker)
library(data.table) # faster data manipulation and improved performance.
##########################


######################################
### READ IN ARGUMENTS FROM BASH SCRIPT
######################################

# Check if an argument was provided
if (length(args) == 0) {
  stop("No chromosome argument provided.")
}
# Get command-line argument for chromosome
chromosome <- commandArgs(trailingOnly = TRUE)[1]  # Read the argument from the Bash script, eg "2L", "2R", "3L", "3R","X"

#####################################

##########################
### MAIN ANALYSIS
##########################

# 1. Load data and subset by chromosome of interest, and format data for subsetting treatment and time-range

load("../data/raw/orch2021_Downsampled_META_Filtered.RData", verbose = TRUE)

# Select specific chromosome
sites.afmat <- cbind(sites, afmat) %>% # Combine "sites" and "afmat" data frames
  filter(chrom == chromosome) # Filter chromosome based on the argument


# Transpose the dataframe and set the "pos" values as column names
sites.afmat <- sites.afmat[, -1]  # remove the "chrom" column
t_sites.afmat <- t(sites.afmat)   # transpose the site.afmat
colnames(t_sites.afmat) <- sites.afmat$pos #set the "pos" values as column names
t.sites.afmat<- t_sites.afmat[-1,] #remove row pos
t.sites.afmat<- as.data.frame(t.sites.afmat, stringsAsFactors = FALSE) #convert to dataframe


# Format a new dataframe which allows selecting a specific treatment & tpt 
samps.t.sites.afmat <- cbind(samps, t.sites.afmat) %>%  # combine "samps" and transposed data frame
  filter(tpt %in% c("1","3","5","7", "9", "10","11","12")) %>% # keep only matched tpts between E and P cages
  mutate(tpt = case_when( # turn the tpt column into numeric and re-name tpt to have continuous numbers
    tpt == "1" ~ 1, tpt == "3" ~ 2, tpt == "5" ~ 3, tpt == "7" ~ 4, tpt == "9" ~ 5, tpt == "10" ~ 6, 
    tpt == "11" ~ 7, tpt == "12" ~ 8, TRUE ~ as.numeric(tpt))) %>%
  mutate(generation = case_when( # add a column for generation time - necessary for rd.site
    tpt == 1 ~ 6, tpt == 2 ~ 7, tpt == 3 ~ 8, tpt == 4 ~ 9, tpt == 5 ~ 12, tpt == 6 ~ 13, 
    tpt == 7 ~ 14, tpt == 8 ~ 15), .after = tpt) %>%
  filter(treatment %in% c("E", "P")) %>% # keep only E and P cages 
  filter(biol.rep == "No" & tech.rep == "No") %>%  # remove biol and tech reps
  filter(!(cage %in% c("E11", "E12")))  # exclude cages E11 and E12


# 2. Fit quasi-binomial GLM (model: counts ~ time) for each treatment/time-range 

print("Fitting GLMs...")

source("./quasi_GLM_mdl1.R")

# Create a vector of treatments and a vector of time point ranges
treatments <- unique(samps.t.sites.afmat$treatment)
timepoint_ranges <- list(c(1,2), c(2,6), c(6,8)) #run code in time ranges based on malathion application

# Create an empty list to store the results
results_list <- list()

# Create a dataframe with chromosome-specific arguments for rd.site
chrom.args <- data.frame(
  chromosome = c("2L", "2R", "3L", "3R", "X"),
  pct_missing = c(6.81, 6.15, 6.85, 7.06, 7.96),
  nof_snps = c(620700, 482594, 576993, 596326, 335095),
  chrom_length = c(23011544, 21146708, 24543557, 27905053, 22422827),
  recomb_rate = c(0.0000000239, 0.0000000266, 0.0000000179, 0.0000000196, 0.0000000295),
  stringsAsFactors = FALSE
)

# Function to retrieve chromosome-specific arguments
get_chrom_args <- function(chromosome) {
  args <- chrom.args[chrom.args$chromosome == chromosome, ]
  return(args)
}

# Loop over the treatments
results_list <-lapply(treatments, function(treatment) {
  # Loop over the time point ranges
  res <- lapply(timepoint_ranges, function(tp_range) {
    # Subset the data based on the time point range and treatment
    subset_data <- samps.t.sites.afmat %>% filter(treatment == !!treatment & tpt >= tp_range[1] & tpt <= tp_range[2])
    
    # Print the current treatment and time range
    print(paste("Processing:", treatment, tp_range[1], tp_range[2]))
    
    # Fit the GLM model with chromosome-specific arguments
    args <- get_chrom_args(chromosome)
    
    # Fit the GLM model with chromosome-specific arguments
    result <- fit_GLM(
      af.site = subset_data %>% dplyr::select(9:ncol(.)),
      rd.site = calc_expected_ec(
        rd = 8,
        gen = subset_data$generation,
        pct_missing = args$pct_missing, 
        nof_snps = args$nof_snps,         
        chrom_length =args$chrom_length, 
        recomb_rate = args$recomb_rate    
      ),
      sampData = subset_data %>% dplyr::select(sample, tpt, treatment),
      formulaString = "cts ~ tpt",
      poolSize = 100,
      numCores = 19
    )

    # Add the result to the list
    treat_range <- paste(treatment, paste(tp_range, collapse = "_"), sep = ".")
    result$treat_range <- treat_range
    results_list[[treat_range]] <- result
  })
  
  res
})

print("GLMs fitted...")


# 3. Combine summary model data in one dataframe 
result_df <- bind_rows(results_list, .id = "treatment") %>%
  dplyr::mutate(chrom = chromosome)  # Add a column with the chromosome value 

# Reset row names
rownames(result_df) <- NULL

# 4. Apply repeat masker

# Use the read_rm function which reads an organism specific Repeat Masker output file

repeat.masker<- read_rm("/scratch/users/mkarag/dm3.fa.out") %>%
  mutate(qry_id = str_replace(qry_id, "^chr", ""))

setDT(repeat.masker) #turn to data.table
setDT(result_df) # convert result_df to a data.table

result_df[, column := as.integer(column)]
repeat.masker[, qry_start := as.integer(qry_start)]
repeat.masker[, qry_end := as.integer(qry_end)]

result_df <- result_df[!repeat.masker, on = .(chrom = qry_id, column >= qry_start, column <= qry_end)] #column is for position

# 5. Save the data

# Define the filename with the chromosome argument
filename <- paste0("./glm_malathion_", chromosome, ".csv")

# Write the dataframe to a CSV file with the chromosome argument
write.csv(result_df, filename, row.names = FALSE)
