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
- library(rjson)
- library(ggplot2)
- library(ggpubr)
- library(Cairo)
- library(gridExtra)
- library(effsize)
- library(BayesFactor)
- library(rjson)
- library(reshape)
- library(lme4)
- library(permutes)

## Raw data
raw data can be found at **'[Python]PreProcessing/results'**

## Pre-processing
- Raw data (.asc) are pre-processed by **'[Python]PreProcessing/parseData.py'**

	- Pre- processed data is saved as **‘s[subject ID]_trial.json’**
	
- Artifact rejection and data epoch are performed by **'[Python]PreProcessing/dataAnalysis.py'**

## Figure and statistics
- *‘[Rmd]Results/figure.Rmd’* is to generate figures and statistical results.


### Article information

