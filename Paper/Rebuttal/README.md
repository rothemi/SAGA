# Rebuttal - Phase 3

Here are two scripts and one Supplementary Table generate for the current rebuttal phase of the manuscript.

## REBUTTAL - R SAGA FeatureSelection LeaveOneBatchOut Reviewer

In this R-Script the leave-one-batch-out scenario is implemented as suggested by the reviewer.

### Dataset preparation
Raw intensities of 169 arrays from 19 experimental batches were read in and combined into an “EListRaw” object without further modification. 15 samples with unknown ground truth were subsequently removed from the dataset, resulting in 154 assays including two mock duplicates (X6374.1, X6379.1 from batch 17). For iteration 1, the raw data of all micorarrays belonging to batch 1 (IVIM #120411) was set aside as an independent test set, all other batches (2-19) were used as training set and were quantile normalized, averaged and batch corrected using limma and ComBat, respectively. 

### Feature Selection
Each preprocessed training set was subjected to unsupervised filtering (genefilter), SVM-RFE and SVM-GA (Caret) as performed for the development phase. The numbers of subset sizes to assess during SVM-RFE encompassed 1,2,3…,40,45,50, all predictors = 43 predictor subsets in total. 

### Prediction of independent test batches
After having determined the optimal predictors in the training set, an SVM with radial kernel was trained on the training set reduced to the 8 optimal predictors found by the feature selection routines. The hyperparameters of the SVM (sigma and cost) were determined by 20 times repeated 10-fold cross-validation. The microarrays of test set 1 (IVIM #120411) were quantile normalized, and averaged using limma. No batch correction was performed. The test set was reduced to the 8 optimal predictors determined on the training set and predicted using the SVM trained before and the caret function “predict”. The unadjusted test set was predicted using random forest, generalized linear models, linear discriminant analysis, KNN and quadratic discriminant analysis (qda) after training the different classifiers on the 8-optimal- predictor training set

*	[REBUTTAL_R_SAGA_FeatureSelection_LeaveOneBatchOut_Reviewer.RData](https://www.dropbox.com/s/pyktuolm5k6gdyu/.RData?dl=0)
*	[Rebuttal_FeatureSelection_batchwise_154_20200724.R](./Rebuttal_FeatureSelection_batchwise_154_20200724.R)

## REBUTTAL - R SAGA FeatureSelection LeaveOneBatchOut bapred upfront
This R-Script is similar to “FIGURE 4_R_SAGA_LeaveOneBatchOut Phase”, except that the training set is quantile-normalized and batch corrected using the R-package “bapred”. Details are given in the rebuttal letter.

*	[REBUTTAL_R_SAGA_FeatureSelection_LeaveOneBatchOut_bapred_upfront.RData](https://www.dropbox.com/s/v3amnhl8oav3ekh/.RData?dl=0)
*	[20200807_154_bapred upfront_batchwise_FINAL.R](./20200807_154_bapred%20upfront_batchwise_FINAL.R)

## REBUTTAL - Supplementary Data 1
Details are given in the rebuttal letter.

*	[Rebuttal_Supplementary_Data_1.xlsx](./Rebuttal_Supplementary_Data_1.xlsx)
