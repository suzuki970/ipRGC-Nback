## Experimental data for *"Selective Activation of ipRGC Modulates Working Memory Performance"*
Copyright 2025 Yuta Suzuki


### Article information
Suzuki, Y., Nakauchi, S., Liao, H. Selective Activation of ipRGC Modulates Working Memory Performance


## Requirements
Python
- pre-peocessing (**'/[Python]PreProcessing/PupilAnalysisToolbox'**)
- numpy
- scipy
- os
- json

R
Run the following code first to install the required packages.
```
pkg_list <- read.csv("loaded_packages_versions.csv", stringsAsFactors = FALSE)

for (i in 1:nrow(pkg_list)) {
  if (!pkg_list[i,"Package"] %in% installed.packages()[,"Package"]){
    install.packages(pkg_list[i,"Package"])
  }
}
```
## Raw data
Raw data can be found at **'[Python]PreProcessing/results'**

## Pre-processing
- Raw data (.asc) are pre-processed by **'[Python]PreProcessing/parseData.py'**

	- Pre- processed data is saved as **‘s[subject ID]_trial.json’**
	
- Artifact rejection and data epoch are performed by **'[Python]PreProcessing/dataAnalysis.py'**

## Figure and statistics
- *‘[Rmd]Results/figure.Rmd’* is to generate figures and statistical results.


### Article information
Suzuki, Y., Nakauchi, S. & Liao, H.-I. Selective activation of ipRGC modulates working memory performance. PLOS One 20, e0327349 (2025).
  
