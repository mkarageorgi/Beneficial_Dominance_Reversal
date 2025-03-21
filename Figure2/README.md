## Figure 2
This code reproduces Figure 2.

### How to generate the figure
* First, obtain the relevant file within the data folder. Place the file under `data/`
* Run each R markdown document to generate the figure plot

### Programming environment
The R scripts require the following packages: `c("tidyverse", "data.table", "broom", "boot", "boot.pval", "multcomp", "emmeans", "msm", "drc", "modelsummary", "ggpubr", "ggthemes", "ggnewscale", "kableExtra", "htmltools", "webshot")`.

These packages serve different purposes:
* Data manipulation: `tidyverse`, `data.table` 
* Statistical analysis: `broom`, `boot`, `boot.pval`, `multcomp`, `emmeans`, `msm`, `drc`, `modelsummary`  
* Visualization: `ggpubr`, `ggthemes`,`ggnewscale`
* Output formatting: `kableExtra`, `htmltools`, `webshot`

Data collected by Zach Mouza, Anastasia Andreeva, and Marianna Karageorgi.  
Code and analysis by Marianna Karageorgi

StatTable.R written by Danae Papadopetraki (https://github.com/danae1968)
