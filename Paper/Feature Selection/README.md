# Feature Selection
## SAGA classifier development and implementation
We implemented the classifier development steps using a training-evaluation scheme, where the whole training cohort (68 samples; 34 transforming, 17 non-transforming; 17 mock samples) is split into multiple training and external evaluation sets by repeated cross-validation. For high dimensional datasets it is key to determine the optimal features for the classification model because they have a direct impact on the performance of many machine learning algorithms, including support vector machines. They can be substantially impaired by large numbers of irrelevant predictors45. Furthermore, models using less predictors are quicker to compute, less prone to overfitting and generally better interpretable than models based on thousands of features27. We first implemented a variance-based filter method to exclude probes showing very little variation in the dataset. We used the R-package “genefilter” to select probes showing an interquartile range of log2-expression values greater than 0.8 in the quantile-normalized and batch corrected training cohort, which selected 1188 out of 35492 annotated features. Next, we performed recursive feature elimination28,46 (SVM-RFE) as implemented in the R-package “caret” version 6.0-7747. To control for overfitting the search routine of SVM-RFE was embedded in an outer resampling layer using five times repeated 10-fold cross-validation. SVM-RFE initially ranks all predictors (n=1188) in each resample according to their individual receiver operating characteristic (ROC) in the training data. In each iteration in each subset, less important predictors are removed and the prediction performance on the hold-out samples is calculated (external validation). The metric to be maximized by the SVM-RFE search algorithm was set to “ROC”. For SVM-RFE we trained an SVM with a radial basis function kernel and tested subset sizes of 1-20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400 and 500 predictors on the mean-centered and scaled training set. The value for the hyperparameter c (cost of violating the margin of the hyperplane) of the SVM was determined in each fold by nested 10-fold cross-validation using the parameter “tune.length” set to 10. The hyperparameter sigma, which determines the width of the Gaussian kernel, is estimated by caret-svm using the function “sigest” of the R-package “kernlab”. After completion of the calculation of the resamples caret aggregates the results of the hold-out predictions to determine the best subset size and repeats the process using the complete training data in order to find the best performing probe subset. The initial results suggested 10-20 probes performed best in conjunction with the radial SVM classifier (Fig. 2i). SVM-RFE is a greedy algorithm designed to quickly discard large numbers of features, therefore it can erroneously reject some informative probes as well as retain less informative ones. Therefore, we repeated SVM-RFE twenty times using different random seeds tabulating the probes most often selected by the algorithm, resulting in 29 probes. Eleven of them were selected in each of the different SVM-RFE runs (Supplementary Table 4). As input for the genetic algorithm feature search we chose all probes that were selected more than twice out of twenty runs (“Top19rfe”; Supplementary Table 4). The genetic algorithm optimization procedure was implemented using the function “gafs” from the R package “caret”. Computations for the genetic algorithm driven feature search were performed on a c4.8xlarge Amazon Web Service EC2 instance with 36 cores running an RStudio Server Amazon Machine Image (AMI) with RStudio 1.1.383 and R 3.4.1. The algorithm was set to maximize internal accuracy as internal and external fitness measure. Parameters used for the genetic algorithm were: generations: 100, population size: 50, crossover probability: 0.8, mutation probability: 0.1, elitism: 3. Similar to SVM-RFE the feature search was conducted repeatedly in an outer layer of resampling iterations which was set to five times 10-fold cross-validation, using one tenth of the data to measure external prediction accuracy in each fold. The computations for the outer resampling were parallelized over the 36 cores of the EC2 instance. For reasons of computational efficacy a support vector machine with linear kernel was used for class- prediction during the genetic algorithm search, whose hyperparameter c (cost) was estimated by nested three times cross-validation with the parameter “tune.length” set to 5. During the search of the feature space the genetic algorithm found three different probe combinations that achieved an external prediction accuracy of more than 97%. Two subsets consisted of eight probes and differed only by one probe. For reasons of stability we chose the third solution (iteration 56; “Top9GA”; Supplementary Table 4) comprised of 9 probes, that was identical to the second solution, but included the probe “A_52_P663904“, interrogating Lhfpl1, the only gene always selected by SVM-RFE that was downregulated in samples transduced with transforming vectors. Using the reduced 9- probe list from the genetic algorithm the final model was built using the package “e1071” to train a SVM with radial basis function kernel on the complete set of training samples. Hyperparameters gamma (g) and cost (c) were optimized by grid search using the function “tune” from “e1071” with an initial range of 10-4:100 for g and 10-1:102 for c and five times 10-fold crossvalidation. We found that ranges for gamma= 0.05…0.15 and cost = 1…10 resulted in stable prediction performance. The SVM was subsequently trained and implemented in the SAGA- package with the optimal combination of hyperparameters using the function “svm” with the parameter “probability” set to “TRUE”. In the R-package SAGA, unknown samples are add-on quantile normalized and add-on batch corrected using the functions “qunormaddon” and “combatbaaddon” of the package “bapred”48 to prevent changes of the training data set by the addition of unknown samples to be predicted by the SAGA classifier. The probability of unknown samples belonging to the class “transforming” or “non-transforming” are computed with the function “predict” of the R package e1071. A vector was considered belonging to the class “transforming” when the probability computed by the SVM model was greater than 0.5. 
For the implementation of SAGA-GSEA assays were read in batch-wise, quantile-normalized, averaged and log2-transformed within each assay using the R package “limma”. The preprocessed and unfiltered expression matrix with the Agilent ProbeIDs as row names was directly converted into an “epheno” object using the function “ExpressionPhenoTest” from the package “phenoTest” 49 with the phenotype variable (“Group”) set to 1 for all mock samples in each assay  and a unique value {2,3,…,n} for each of the samples to be tested against the mock samples. As probeset to be tested for enrichment we took the ProbeIDs that were selected by SVM-RFE as described above. A_52_P663904 / Lhfpl1 was excluded from 19-probeset as found by SVM-RFE since it was downregulated in transformed samples, leaving 18 upregulated probes as SAGA-Core probeset. Of note, we did not use the 9 probeset selected by SVM-GA, since very small probesets (<15) can lead to instability of the GSEA-algorithm22. The normalized enrichment score and its p-values were calculated for every sample against the mock samples using the function “gsea” from “phenoTest”. A vector was assigned to the class “transforming” when its NES was greater than 0.

## Performance of SAGA classifiers
In order to assess the reproducibility and stability of the SAGA pipeline the final classifier was tested with an independent validation set (23 samples; 11 transforming, 8 non-transforming; 4 mock) not involved in any of the classifier construction procedures. Prediction of the unknown samples was performed as described above after the classifier had been trained and the optimal hyperparameters had been found using the complete training cohort. Performance was assessed using the function “confusionMatrix” from the package “caret” to tabulate observed and predicted classes and compute sensitivity, specificity, accuracy and Cohen’s kappa. The receiver operating characteristic was calculated with the function “roc” from the package “pROC”. To assess performance of SAGA we used 5 times repeated 10-fold cross-validation of the combined training and validation set. For each of the 50-folds a radial SVM was trained und tuned as described above for the complete SAGA classifier. We predicted the hold-out samples as described above and estimated the performance using the function “confusionMatrix” and the AUC of the receiver operating characteristic. The procedure was repeated with the 19- probeset as found by SVM-RFE and the 1188 probes used as input into the feature selection routines yielding the ROC-curves depicted in Figure 2m. 

## Availability of raw data and R workspace file

The following files can be downloaded [here](https://owncloud.gwdg.de/index.php/s/nHjGbn5sAnU93SV):
*	Feature_Selection.RData
*	[20190306_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R](./20190306_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R)
*	20190220_rfe_GA_152_10TestSets and FINALSET_Download.zip
*	20190306_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R
*	AllResults_Accuracy_full_TestSet vs Resampling_IQR0.8.pdf
*	AllResults_Accuracy_full_TestSet vs Resampling_IQR0.8_Median.pdf
*	AllResults_Accuracy_GA_TestSet vs Resampling_IQR0.8.pdf
*	AllResults_Accuracy_GA_TestSet vs Resampling_IQR0.8_median.pdf
*	AllResults_Accuracy_rfe_TestSet vs Resampling_IQR0.8.pdf
*	AllResults_Accuracy_rfe_TestSet vs Resampling_IQR0.8_median.pdf
*	AllResults_rfe_GA40_10TestSets_FinalSet152_IQR0.8.txt
*	AllTestSets_ConfusionMatrix_SAGA.txt
*	AllTestSets_LTR.RV_SAGA_vs_IVIM.pdf
*	AllTestSets_Other_SAGA_vs_IVIM.pdf
*	AllTestSets_ROC_LTR.SF vs other.R
*	AllTestSets_ROC_SAGA.pdf
*	AllTestSets_ROC_SAGA_RFA_GA_IVIM.pdf
*	AllTestSets_ROC_SAGA_vs_IVIM.pdf
*	AllTestSets_ROC_SVM.Full.pdf
*	AllTestSets_ROC_SVM.GA.pdf
*	AllTestSets_ROC_SVM.rfe.pdf
*	AllTestSets_SAGA_vs_IVIM152.pdf
*	AllTestSets_SAGA_vs_IVIM152_LTR.SF.pdf
*	AllTestSets_SAGA_vs_IVIM152_noLTR.SF.pdf
*	AllTestSetsConfusionMatrix_FULL.txt
*	AllTestSetsConfusionMatrix_rfe.txt
*	Annotation.matrix.train.FINAL.txt
*	Annotation_SAGA_FINAL_KNOWN_20181128.txt
*	DeLongTest_Compound vs IVIM.txt
*	DeLongTest_GA vs IVIM.txt
*	DeLongTest_RFE vs IVIM.txt
*	DeLongTest_SAGA vs IVIM.txt
*	ESET_RMA_COMBAT_KNOWN_FullSagaSet152_FINAL.txt
*	FINAL_VariableImportance_ROC_IQR0.8.txt
*	FinalSet152_Accuracy over FeatureNumber.pdf
*	FinalSet152_Accuracy over FeatureNumber_Zoom.pdf
*	FinalSet152_Accuracy over Iteration_GA.pdf
*	FinalSet152_Accuracy over Iteration_GA_external.pdf
*	FinalSet152_ModelDifferences_rfe vs full.txt
*	FinalSet152_optVars_GeneticAlgorithm.txt
*	FinalSet152_optVars_rfe.txt
*	FinalSet152_PCA_allVars.pdf
*	FinalSet152_PCA_allVars_V2.pdf
*	FinalSet152_PCA_optVars_GA.pdf
*	FinalSet152_PCA_optVars_GA_V2.pdf
*	FinalSet152_Resamples_rfe vs full.txt
*	FinalSet152_Results_rfe.txt
*	IVIM_MTT_ROC.txt
*	IVIM_Total_ConfusionMatrix.txt
*	SAGA_Targets_FINAL_152.txt
*	TestSet01_ConfusionMatrix_allVars.txt
*	TestSet01_ConfusionMatrix_optVars.txt
*	TestSet01_PCA_allVars.pdf
*	TestSet01_PCA_optVars_rfe.pdf
*	TestSet01_Predictions_allVars.txt
*	TestSet01_Predictions_rfe.txt
*	TestSet01_ROC_allVars.pdf
*	TestSet01_ROC_optVars_rfe.pdf
*	TestSet010_ROC_allVars.pdf
*	TestSet010_ROC_optVars_GA.pdf
*	TestSet010_ROC_optVars_rfe.pdf
*	TestSet02_ConfusionMatrix_FULL.txt
*	TestSet02_ConfusionMatrix_GA.txt
*	TestSet02_ConfusionMatrix_rfe.txt
*	TestSet02_Predictions_allVars.txt
*	TestSet02_Predictions_optVars_GA.txt
*	TestSet02_Predictions_optVars_rfe.txt
*	TestSet02_ROC_allVars.pdf
*	TestSet02_ROC_optVars_GA.pdf
*	TestSet02_ROC_optVars_rfe.pdf
*	TestSet03_ConfusionMatrix_FULL.txt
*	TestSet03_ConfusionMatrix_GA.txt
*	TestSet03_ConfusionMatrix_rfe.txt
*	TestSet03_Predictions_allVars.txt
*	TestSet03_Predictions_optVars_GA.txt
*	TestSet03_Predictions_optVars_rfe.txt
*	TestSet03_ROC_allVars.pdf
*	TestSet03_ROC_optVars_GA.pdf
*	TestSet03_ROC_optVars_rfe.pdf
*	TestSet04_ConfusionMatrix_FULL.txt
*	TestSet04_ConfusionMatrix_optVars_rfe.txt
*	TestSet04_Predictions_allVars.txt
*	TestSet04_Predictions_optVars_rfe.txt
*	TestSet04_ROC_allVars.pdf
*	TestSet04_ROC_optVars_rfe.pdf
*	TestSet05_ConfusionMatrix_FULL.txt
*	TestSet05_ConfusionMatrix_optVars_rfe.txt
*	TestSet05_Predictions_allVars.txt
*	TestSet05_Predictions_optVars_rfe.txt
*	TestSet05_ROC_allVars.pdf
*	TestSet05_ROC_optVars_rfe.pdf
*	TestSet06_ConfusionMatrix_FULL.txt
*	TestSet06_ConfusionMatrix_GA.txt
*	TestSet06_ConfusionMatrix_rfe.txt
*	TestSet06_Predictions_allVars.txt
*	TestSet06_Predictions_optVars_GA.txt
*	TestSet06_Predictions_optVars_rfe.txt
*	TestSet06_ROC_allVars.pdf
*	TestSet06_ROC_optVars_GA.pdf
*	TestSet06_ROC_optVars_rfe.pdf
*	TestSet07_ConfusionMatrix_FULL.txt
*	TestSet07_ConfusionMatrix_GA.txt
*	TestSet07_ConfusionMatrix_rfe.txt
*	TestSet07_Predictions_allVars.txt
*	TestSet07_Predictions_optVars_GA.txt
*	TestSet07_Predictions_optVars_rfe.txt
*	TestSet07_ROC_allVars.pdf
*	TestSet07_ROC_optVars_GA.pdf
*	TestSet07_ROC_optVars_rfe.pdf
*	TestSet08_ConfusionMatrix_FULL.txt
*	TestSet08_ConfusionMatrix_GA.txt
*	TestSet08_ConfusionMatrix_rfe.txt
*	TestSet08_Predictions_allVars.txt
*	TestSet08_Predictions_optVars_GA.txt
*	TestSet08_Predictions_optVars_rfe.txt
*	TestSet08_ROC_allVars.pdf
*	TestSet08_ROC_optVars_GA.pdf
*	TestSet08_ROC_optVars_rfe.pdf
*	TestSet09_ConfusionMatrix_FULL.txt
*	TestSet09_ConfusionMatrix_optVars_rfe.txt
*	TestSet09_Predictions_allVars.txt
*	TestSet09_Predictions_optVars_rfe.txt
*	TestSet09_ROC_allVars.pdf
*	TestSet09_ROC_optVars_rfe.pdf
*	TestSet10_ConfusionMatrix_FULL.txt
*	TestSet10_ConfusionMatrix_GA.txt
*	TestSet10_ConfusionMatrix_rfe.txt
*	TestSet10_Predictions_allVars.txt
*	TestSet10_Predictions_optVars_GA.txt
*	TestSet10_Predictions_optVars_rfe.txt
*	TestSets2367810_ConfusionMatrix_Full.txt
*	TestSets2367810_ConfusionMatrix_GA.txt
*	TestSets2367810_ConfusionMatrix_RFE.txt
*	TestSets2367810_ROC_full_rfe_GA.pdf
*	TrainingSet01_Accuracy over FeatureNumber.pdf
*	TrainingSet01_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet01_ModelDifferences_rfe vs full.txt
*	TrainingSet01_OptVars_rfe.txt
*	TrainingSet01_Resamples_rfe vs full.txt
*	TrainingSet01_Results_rfe.txt
*	TrainingSet010_Accuracy over FeatureNumber.pdf
*	TrainingSet010_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet010_Accuracy over Iteration_GA.pdf
*	TrainingSet010_Accuracy over Iteration_GA_external.pdf
*	TrainingSet02_Accuracy over FeatureNumber.pdf
*	TrainingSet02_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet02_Accuracy over Iteration_GA.pdf
*	TrainingSet02_Accuracy over Iteration_GA_external.pdf
*	TrainingSet02_ModelDifferences_rfe vs full.txt
*	TrainingSet02_optVars_GeneticAlgorithm.txt
*	TrainingSet02_OptVars_rfe.txt
*	TrainingSet02_Resamples_rfe vs full.txt
*	TrainingSet02_Results_rfe.txt
*	TrainingSet03_Accuracy over FeatureNumber.pdf
*	TrainingSet03_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet03_Accuracy over Iteration_GA.pdf
*	TrainingSet03_Accuracy over Iteration_GA_external.pdf
*	TrainingSet03_ModelDifferences_rfe vs full.txt
*	TrainingSet03_optVars_GeneticAlgorithm.txt
*	TrainingSet03_optVars_rfe.txt
*	TrainingSet03_Resamples_rfe vs full.txt
*	TrainingSet03_Results_rfe.txt
*	TrainingSet04_Accuracy over FeatureNumber.pdf
*	TrainingSet04_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet04_ModelDifferences_rfe_vs_full.txt
*	TrainingSet04_optVars_rfe.txt
*	TrainingSet04_Resamples_rfe_vs_full.txt
*	TrainingSet04_Results_rfe.txt
*	TrainingSet05_Accuracy over FeatureNumber.pdf
*	TrainingSet05_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet05_ModelDifferences_rfe_vs_full.txt
*	TrainingSet05_optVars_rfe.txt
*	TrainingSet05_resamples_rfe_vs_full.txt
*	TrainingSet05_Results_rfe.txt
*	TrainingSet06_Accuracy over FeatureNumber.pdf
*	TrainingSet06_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet06_Accuracy over Iteration_GA.pdf
*	TrainingSet06_Accuracy over Iteration_GA_external.pdf
*	TrainingSet06_ModelDifferences_rfe_vs_full.txt
*	TrainingSet06_optVars_GeneticAlgorithm.txt
*	TrainingSet06_optVars_rfe.txt
*	TrainingSet06_resamples_rfe_vs_full.txt
*	TrainingSet06_Results_rfe.txt
*	TrainingSet07_Accuracy over FeatureNumber.pdf
*	TrainingSet07_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet07_Accuracy over Iteration_GA.pdf
*	TrainingSet07_Accuracy over Iteration_GA_external.pdf
*	TrainingSet07_ModelDifferences_rfe vs full.txt
*	TrainingSet07_optVars_GeneticAlgorithm.txt
*	TrainingSet07_optVars_rfe.txt
*	TrainingSet07_resamples_rfe vs full.txt
*	TrainingSet07_Results_rfe.txt
*	TrainingSet08_Accuracy over FeatureNumber.pdf
*	TrainingSet08_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet08_Accuracy over Iteration_GA.pdf
*	TrainingSet08_Accuracy over Iteration_GA_external.pdf
*	TrainingSet08_Modeldifferences_rfe vs full.txt
*	TrainingSet08_optVars_GeneticAlgorithm.txt
*	TrainingSet08_optVars_rfe.txt
*	TrainingSet08_resamples_rfe vs full.txt
*	TrainingSet08_Results_rfe.txt
*	TrainingSet09_Accuracy over FeatureNumber.pdf
*	TrainingSet09_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet09_ModelDifferences_rfe_vs_full.txt
*	TrainingSet09_optVars_rfe.txt
*	TrainingSet09_Resamples_rfe_vs_full.txt
*	TrainingSet09_Results_rfe.txt
*	TrainingSet10_ModelDifferences_rfe vs full.txt
*	TrainingSet10_optVars_GeneticAlgorithm.txt
*	TrainingSet10_optVars_rfe.txt
*	TrainingSet10_Resamples_rfe vs full.txt
*	TrainingSet10_Results_rfe.txt
