## Figure 4

This code reproduces Figure 4.

How to generate the figure:

* First, obtain the following files that were generated earlier: the allele frequency data (`orch2021_Downsampled_META_Filtered.RData`) and inbred line SNP tables for each chromosome (`Orchard2021/inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric`). Place the first file under `data/raw/` (`data/raw/orch2021_Downsampled_META_Filtered.RData`). Place the SNP tables under `data/snptables/Orchard2021/` (`inbredv2_withHets.orch2021.{chromosome}.snpTable.numeric`).

* Run the code in the `glm` subfolder to generate GLM effect size estimates. A repeat masker file `dm3.fa.out` is required for the GLM. Obtain this file from and put it under `data/raw/` as `data/raw/dm3.fa.out`. First, run `01_glm_malation.R`, passing each chromosome arm as an argument:

```
Rscript --vanilla 01_glm_malation.R 2L
Rscript --vanilla 01_glm_malation.R 2R
Rscript --vanilla 01_glm_malation.R 3L
Rscript --vanilla 01_glm_malation.R 3R
Rscript --vanilla 01_glm_malation.R X
```

These scripts are highly parallelized and computationally intensive, so a cluster environment may be neccessary. Then, run `02_summarize_glm.R`. This script generates a file `data/raw/sigsite_malation.csv` that would be used in downstream computations.

* Run all numbered R and Python scripts in this directory in order. These scripts will generate small tables in the `plot_data` folder that are used for plotting.

* Run `plot.Rmd` to generate the figure panels.


Code and analysis by Marianna Karageorgi and Egor Lappo

### Programming environment

The R scripts require the following packages: `c("tidyverse", "reticulate", "broom", "egg", "ggnewscale", "purrr", "furrr", "biomartr")`.

The Python programming environment should be created with `conda` using the `env.yml` file in this directory.