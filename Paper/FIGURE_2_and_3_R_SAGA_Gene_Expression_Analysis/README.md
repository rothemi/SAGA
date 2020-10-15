# FIGURE 2 and 3 - R SAGA Gene Expression Analysis

This R-script contains the generic gene expression analysis of the SAGA dataset. 

## Dataset preparation
First, all raw files (.txt)  were read in separately for each array design with the function “read.maimages” (limma). A merged dataset was created by extracting all probes from the original Agilent Mouse Genome Oligo Microarray 4x44K v2 array from the four array platforms and combining them using the function “cbind.EList” (limma). The Raw data was quantile-normalized and log2-transformed followed by averaging of the four within-array replicates of each probe resulting in a dataset with 39,428 unique probes. Batch correction between different SAGA assays was performed on quantile-normalized log2-values using ComBat as implemented in the R package “sva”. 

## t-SNE
Gene expression profiles were visualized by t-distributed stochastic neighbor embedding (t-SNE) using the  Barnes-Hut implementation of t-SNE from the “Rtsne”-package. For each t-SNE plot, Barnes-Hut t-SNE was run 1000 times with different random seeds and the iteration with the lowest Kullback-Leibler divergence was selected for visualization as a 2D plot.

## Differential expression analysis
Differentially expressed probes between the subgroups were computed on the quantile normalized and batch corrected gene expression matrix (36,226 annotated probes) using ”limma” with Benjamini-Hochberg multiple testing correction. We computed the top differentially expressed genes for the following contrasts: ”transforming – mock”, “safe – mock”, “transforming – safe” and “transforming – (mock+safe)/2”  for 152 SAGA samples with known IVIM properties (65 transforming, 32 mock and 55 safe). 

## Gene set enrichment analysis
The quantile normalized and batch corrected gene expression matrix (36,226 annotated probes, 152 samples) was filtered for the 24,664 gene symbols that appear at least once in the interrogated MSigDB.v6.2 gene set collections. In cases with multiple probes per gene, the probe with the highest standard deviation across the samples was selected, resulting in a gene expression matrix consisting of 15,376 probes. From this matrix .gct files were generated and used as input for the Broad GSEA software together with a .chip file containing the annotation for the 15,376 probes. The enrichment results were visualized by plotting the normalized enrichment score (NES) against the FDR. We additionally performed GSEA with ROAST (rotation gene set tests for complex microarray experiments) and CAMERA (competitive gene set test accounting for inter-gene correlation) from the limma package by applying both functions to the matrix of 15,376 probes and computing the same contrasts as with the Broad GSEA tool. 


## Availability of R script and workspace file

*	[FIGURE_2_and_3_R_SAGA_Gene_Expression_Analysis.RData](https://www.dropbox.com/sh/2rlpjnhece4tl8p/AAAVibkBQVRXtV2DRBeeOeIia?dl=0)
*	[R_SAGA_Figure 2 and 3_GeneExpressionAnalysis.R](./R_SAGA_Figure%202%20and%203_GeneExpressionAnalysis.R)
