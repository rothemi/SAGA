# FIGURE 4 - R SAGA LeaveOneBatchOut Phase

This R-script contains the code used for leave-one-batch-out validation phase. 

## Dataset preparation
Raw intensities of 169 arrays from 19 experimental batches were read in and combined into an “EListRaw” object without further modification. 15 samples with unknown ground truth were subsequently removed from the dataset, resulting in 154 assays including two mock duplicates (X6374.1, X6379.1 from batch 17, Supplementary data 8). For iteration 1, the raw data of batch 1 (IVIM #120411) was set aside as an independent test set, all other batches (2-19) were used as training set and were quantile normalized, averaged and batch corrected using limma and ComBat, respectively. 

## Feature Selection
Each preprocessed training set was subjected to unsupervised filtering (genefilter), SVM-RFE and SVM-GA (Caret) as performed for the development phase. The numbers of subset sizes to assess during SVM-RFE encompassed 1,2,3…,40,45,50, all predictors = 43 predictor subsets in total. 

## Prediction of independent test batches
After having determined the optimal predictors in the training set, the raw training set was again quantile normalized and batch-corrected by the R package “bapred” in order to estimate and store the parameters necessary for the later add-on correction of the test set. An SVM with radial kernel was trained on the bapred-adjusted training set reduced to the optimal predictors found by the feature selection routines. The hyperparameters of the SVM (sigma and cost) were determined by 20 times repeated 10-fold cross-validation as before. At this point, the optimal features had been determined and the classifier had been trained and fixed using the training set only, whereas the test set had not been used. This was followed by add-on quantile normalization and add-on batch correction of the raw-test set using the bapred functions “qunormaddon” and “combatbaaddon”, respectively. The add-on adjusted test set was reduced to the optimal predictors determined on the training set and predicted using the SVM trained before and the caret function “predict”. The complete procedure was repeated 18 additional times with every available batch to be used one time as independent test set. The results from the 19 iterations of building SAGA and predicting the independent test batches were aggregated and analyzed by the functions “confusionMatrix” (Caret), “roc” (pROC) and “pr.curve” (PRROC).

## Availability of R script and workspace file

*	[FIGURE_4_R_SAGA_LeaveOneBatchOut_Phase.RData](https://www.dropbox.com/s/7o7wq3e1eifdmzi/Feature%20Selection%20Leave%20One%20Batch%20Out.RData?dl=0)
*	[R_FeatureSelection_LeaveOneBatchOut_154_FINAL_20200312.R](./R_FeatureSelection_LeaveOneBatchOut_154_FINAL_20200312.R)
