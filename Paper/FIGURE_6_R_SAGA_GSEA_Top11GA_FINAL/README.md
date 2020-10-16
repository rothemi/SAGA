# FIGURE 6 - R SAGA GSEA Top11GA FINAL

For the final implementation of SAGA-GSEA to be used in the R-package, the 11 optimal predictors determined on the complete dataset (FIGURE 5_R_SAGA_FINALSET_152) for the final SAGA classifier are used as gene set. 

The raw data of each SAGA assay was read in separately, quantile-normalized, averaged and log2-transformed using the R package “limma”. The preprocessed gene expression matrix was converted into an “epheno” object using the function “ExpressionPhenoTest” from the package “phenoTest”.  The phenotype variable (“Group”) was set to “1” for all mock samples in each assay and to a unique value {2,3,…,n} for each of the samples to be tested against the mock controls. The normalized enrichment score, p-values and fdr were calculated for every sample against the mock samples using the function “gsea” from “phenoTest” using the 11 optimal predictors as gene set. The GSEA results were aggregated over the 18 batches with mock controls available. The ROC curve for SAGA-GSEA and the best normalized enrichment score (NES) cutoff were computed using the function “roc” on the normalized enrichment scores and the true class labels. A vector was assigned to the class “transforming” when its NES was greater than the optimal ROC-cutoff computed on the dataset after exclusion of the strongly transforming LTR.SFFV.eGFP samples (NES=1.0).

## Availability of R script and workspace file

*	[FIGURE_6_R_SAGA_GSEA_Top11GA_FINAL.RData](https://www.dropbox.com/s/oqrw52r29ek4agu/.RData?dl=0)
*	[R_SAGA_GSEA_Top11GA_FINAL.R](./R_SAGA_GSEA_Top11GA_FINAL.R)
