# Extended Data Figure 2
This code reproduces Extended Data Figure 2, which shows dose-response curves for each Ace genotype (a) and estimates LD50 values for each Ace genotype derived from four-parameter logistic dose-response modeling (b).

### How to generate the figure
* First, obtain the relevant file within the data folder. Place the file under `data/`
* Run the R markdown document to generate the figure plot

### Programming environment
The R scripts require the following packages: c("tidyverse", "data.table", "ggpubr", "ggthemes", "drc", "broom", "kableExtra", "htmltools", "webshot").
These packages serve different purposes:

Data manipulation: tidyverse, data.table
Statistical analysis: drc, broom
Visualization: ggpubr, ggthemes
Output formatting: kableExtra, htmltools, webshot

Data collected by Zach Mouza, Anastasia Andreeva, and Marianna Karageorgi.  
Code and analysis by Marianna Karageorgi
