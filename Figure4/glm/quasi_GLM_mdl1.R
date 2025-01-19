## The calc_expected_ec function estimates the expected effective coverage based on a number input parameters.
calc_expected_ec = function(rd, gen, pct_missing, nof_snps, chrom_length, recomb_rate) {
  # Define parameters
  q = 18
  mycoeffs = data.frame(a = 0.5199118, b = -0.6909052, c = 0.3553630)
  
  # Calculate window size based on recombination rate and chromosome length
  winSize = round(qexp(q / 100, 1 / ((chrom_length) / ((recomb_rate) * (chrom_length) * (gen) + 1))) / 1000)
  
  # Calculate number of reads per window based on read depth, number of SNPs, and window size
  nReadsPerWin = (rd) * (nof_snps) * (winSize) * 1000 / (chrom_length)
  
  # Calculate expected effective coverage (ec) based on actual read depth, percent missing, and mycoeffs
  ec = 10^(mycoeffs$a * log10(nReadsPerWin) + mycoeffs$b * log10(1 + pct_missing) + mycoeffs$c)
  
  # Return expected effective coverage (ec)
  return(ec)
}

# calculate effect size where b0 is the intercept b1 is the time coefficient
effect.size.mdl1 <- function(b0, b1, x1, x2) {
  p1 <- exp(b0 + b1 * x1) / (1 + exp(b0 + b1 * x1))
  p2 <- exp(b0 + b1 * x2) / (1 + exp(b0 + b1 * x2))
  dp <- p2 - p1
  
  return(list(freq1 = p1, freq2 = p2, effect_size = dp))
}

## fit quasi-binomial GLM 
fit_GLM <- function(af.site, rd.site, sampData, formulaString, poolSize, cmpAll = NULL, dontReport = NULL, numCores) {
  
  # Calculate the effective sample size based on the total number of alleles (poolSize *2) and coverage per site based on HARP (rd.site)
  Neff <- ((poolSize * 2 * rd.site) - 1) / (poolSize * 2 + rd.site)
  
  # Function to fit a quasibinomial GLM to a column and extract coefficient and p-value information using broom::tidy()
  fit_GLM_one <- function(col) {
    # Calculate the allele counts based on the estimated effective sample size and observed allele frequencies. This is equivalent to cbind(A_Cnt,Tot_Cnt-A_Cnt)
    cts <- cbind(round(Neff * col), round(Neff * (1 - col)))
    
    # Fit a quasibinomial GLM to the data using the specified formula
    model <- glm(as.formula(formulaString), family = "quasibinomial", data = sampData) 
    
    # Extract coefficient and p-value information using broom::tidy()
    tidy_stats <- tidy(model) 
    
    # Calculate effect size
    x1 <- min(sampData$tpt)  
    x2 <- max(sampData$tpt)  
    
    b0 <- tidy_stats$estimate[tidy_stats$term == "(Intercept)"]
    b1 <- tidy_stats$estimate[tidy_stats$term == "tpt"]
    
    effect_size <- effect.size.mdl1(b0, b1, x1, x2)
    
    tidy_stats <- tidy_stats %>%
      mutate(freq1 = effect_size$freq1,
             freq2 = effect_size$freq2,
             effect_size = effect_size$effect_size)    
    
    # Return the tidy statistics data frame
    return(tidy_stats)
  }
  
  # Set the number of cores for parallel processing
  plan(future::multicore, workers = numCores)
  
  # Apply the fit_GLM_one function to each column of af.site using furrr::future_map_dfr()
  result <- furrr::future_map_dfr(af.site, fit_GLM_one, .id = "column")
  
  # Reset the plan to the default (sequential computation)
  plan(future::sequential)
  
  # Return the resulting data frame
  return(result)
  }

##extract coefficients and p values
extract_coef_pval = function(model, cmpAll=NULL, dontReport=NULL) {
  
  # If cmpAll is not NULL, generate pairwise comparisons of the coefficients using the Tukey method
  if (!is.null(cmpAll)) {
    
    # Use the glht() function to generate the comparisons
    model.multcomp = summary(eval(parse(text=paste0("glht(model, mcp(", cmpAll, "='Tukey'))"))))
    
    # Extract the coefficient and p-value information for the comparisons, and format it into a matrix called cp
    cp = cbind(coefficients(model.multcomp), model.multcomp$test$pvalues)
    row.names(cp) = gsub(" - ", "_", row.names(cp))
    
  } else {
    
    # If cmpAll is NULL, extract the coefficient and p-value information for all coefficients in the model using summary()
    cp = summary(model)$coefficients[-1, c(1, 4), drop=FALSE]
    
  }
  
  # If dontReport is not NULL, remove any rows from cp that contain the specified character strings
  if (!is.null(dontReport)) {
    cp = cp[grep(paste0("(", paste0(dontReport, collapse="|"), ")"), row.names(cp), invert=TRUE), , drop=FALSE]
  }
  
  # Return the resulting coefficient and p-value matrix cp
  return(cp)
  
}
