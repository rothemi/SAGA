# FIGURE 5 - R SAGA FINALSET 152

This R-Script repeats the feature selection routines developed before for the complete SAGA dataset (152 samples with known IVIM properties). 

## Dataset preparation
The starting point was the quantile normalized and batch corrected expression matrix (36,226 annotated probes, 152 samples with known IVIM properties).  

## Feature Selection
The data set was subjected to unsupervised filtering (genefilter), yielding 1243 probes at an IQR = 0.8, which were supplied to SVM-RFE. SVM-RFE was performed with the same parameters as described in the “FIGURE 3_R_SAGA_Development Phase” R-Script. SVM-RFE retained 20 predictors, which were used as input for SVM-GA. After 14 iterations the genetic algorithm found an optimal combination of 11 predictors for the complete dataset. Principal component analysis of the SAGA dataset reduced to these 11 probes was performed with the base function “prcomp”. For comparison, PCA was performed on 11 randomly selected probes. 

## Availability of R script and workspace file

*	[FIGURE_5_R_SAGA_FINALSET_152.RData](https://www.dropbox.com/s/k1gj33layaufthg/.RData?dl=0)
*	[20200826_FeatureSelection_rfe_GA40_FinalSet_152.R](./20200826_FeatureSelection_rfe_GA40_FinalSet_152.R)
