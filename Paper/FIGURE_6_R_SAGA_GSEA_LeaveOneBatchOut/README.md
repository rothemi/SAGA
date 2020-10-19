# FIGURE 6 - R SAGA GSEA LeaveOneBatchOut

This R-script contains the code used for leave-one-batch-out validation phase, extended by GSEA for the left-out test batches. 

## Dataset preparation
Raw intensities of 169 arrays from 19 experimental batches were read in and combined into an “EListRaw” object without further modification. 15 samples with unknown ground truth were subsequently removed from the dataset, resulting in 154 assays. For iteration 1, the raw data of all microarrays belonging to batch 1 (IVIM #120411) was set aside as an independent test set, all other batches (2-19) were used as training set and were quantile normalized, averaged and batch corrected using limma and ComBat, respectively. For iteration 2, the raw data of all microarrays belonging to batch 2 (IVIM #150128) were set aside and batches 1,3-19 were used as training set. This splitting scheme was repeated for all 18 batches that contained mock control samples (IVIM #171102 was excluded from SAGA-GSEA since it had no mock samples available). 

## SAGA_GSEA
During the leave-one-batch-out procedure, the optimal predictors found by SVM-RFE and SVM-GA for the training set of each iteration were used as gene set for GSEA. The raw data of the left-out test batch was read in, quantile-normalized, averaged and log2-transformed using the R package “limma”. The preprocessed gene expression matrix of the left-out test batch was converted into an “epheno” object using the function “ExpressionPhenoTest” from the package “phenoTest” with the phenotype variable (“Group”) set to 1 for all mock samples in each assay and a unique value {2,3,…,n} for each of the samples to be tested against the mock controls. The normalized enrichment score, p-values and fdr were calculated for every sample against the mock samples using the function “gsea” from “phenoTest”. The procedure was repeated for each of the 18 test batches with mock controls available. The GSEA results were aggregated over the 18 test batches. The ROC curve for SAGA-GSEA and the best normalized enrichment score (NES) cutoff were computed using the function “roc” on the normalized enrichment scores and the true class labels. A vector was assigned to the class “transforming” when its NES was greater than the optimal ROC-cutoff computed on the dataset after exclusion of the strongly transforming LTR.SFFV.eGFP samples (NES>1.3).

## Availability of R script and workspace file

*	[FIGURE_6_R_SAGA_GSEA_LeaveOneBatchOut.RData](https://www.dropbox.com/s/72ykrdzfxghofuf/.RData?dl=0)
*	[R_FeatureSelection and GSEA_LeaveOneBatchOut_154.R](./R_FeatureSelection%20and%20GSEA_LeaveOneBatchOut_154.R)
