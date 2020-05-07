# Feature Selection
## SAGA classifier development
The development of the predictive model of SAGA was implemented using the R package “caret”52 based on a support vector machine with a radial basis function kernel (method = "svmRadial"). Unless otherwise specified, all calls to functions mentioned in this paragraph belong to the “caret” package with key parameters specified in parentheses after the name of the function or directly discussed in the text. All computations were run on a c5.18xlarge Amazon Web Service EC2 instance with 72 cores and 144 Gb RAM running RStudio 1.1.456 and R 3.5.1. The data splitting and resampling scheme to assess the performance of the models and control for overfitting is outlined in Supplementary Figure 3. First, the quantile normalized and batch corrected expression matrix (36,226 annotated probes, 152 samples with known IVIM properties) was partitioned into a training set comprised of 70% of the samples (107 samples) and an independent test set of 30% of the samples (45 samples). The test set was not used at any point for feature selection or model tuning. To allocate samples to the test or training set the caret function “createDataPartition” (p=0.7) was used, which performs stratified sampling based on the class labels to keep the distribution of transforming and nontransforming samples equal between the training and test sets. Since a single training / test set split can lead to a biased assessment of model building and feature selection28 ten stratified random training / test set splits of the dataset were created and the complete model building pipeline was run for ten times for a more unbiased and reliable assessment of the predictive modeling process. Predictive performance of many models, especially support vector machines, can be significantly affected by large numbers of irrelevant predictors53. Furthermore, models using fewer predictors are quicker to compute, less prone to overfitting and generally better interpretable than models based on thousands of predictors28. Therefore, a combination of feature selection steps was performed to reduce the number of predictors as far as possible while maintaining or increasing predictive power. First, we applied an unsupervised filter to each training set to exclude probes interrogating genes that were not expressed at all or show only little variation in the dataset. This step helped to reduce computation time and avoided the selection of features by the subsequent SVM-RFE step that have a good discriminatory power between the classes based on their AUROC, but display only a small absolute fold-change between the different classes. The R-package “genefilter” was used to discard probes with an interquartile range (IQR) of log2-expression values less than 0.8 in the quantile-normalized and batch corrected training cohort, which retained a median of 1,195 out of 36,226 annotated probes (Supplementary Data 4 tab 1). IQR = 0.8 was chosen empirically, since it consistently selected around 1,000 features in all test/training set splits. Setting the IQR lower (e.g. IQR= 0.5) retained too many features (median around 4,500), leading to a substantial increase in overall computation time as well as a failure to reduce the number of features in the subsequent SVM-RFE step in 3 out of 10 training/test splits. In contrast, setting IQR=1.2 selected on average around 250 features, which could be efficiently handled by SVM-RFE. However, at IQR=1.2 important predictors, such as A_55_P2077048/Itih5 (AUROC= 0.98) were already discarded before the actual feature selection step. The implementations using IQR 0.5 and IQR = 1.2 are available at GITHUB.  Next, we performed recursive feature elimination (SVM-RFE) on the training set using the function “rfe”. Since feature selection is part of the model building process, it needs to be conducted inside of a resampling layer (“external resampling layer”, Supplementary Figure 3) to assess the impact of the selection process on the model performance and to prevent overfitting of the model to the predictors. To establish the external resampling layer, 200 resamples of the training set were created by twenty times repeated 10-fold cross-validation using the function “createMultiFolds” (Parameters: k=10, times = 20). The function divides the entire training set (107 samples) into 10 subsets (folds) of equal size and the first fold (11 samples, “external holdouts”) is predicted by a model fit to the remaining 9 folds (96 samples, “external training”) of the data. This is repeated with the second fold after the first one has been returned to the training set and so on, resulting in 10 resamples for each of the twenty repeats of 10-fold CV. Importantly, the 200 identical resamples were used to fit the full models using all predictors, to allow a direct comparison of the SVM-RFE model and the full model using the resampling accuracies. The 200 resamples were submitted to the helper function “rfeControl”, which controls the details of the external resampling process of the function “rfe”. The feature selection process itself was carried out for each of the 200 resamples separately and computed in parallel by setting the “rfeControl” parameter: “allowParallel = TRUE”. To ensure reproducibility of the analysis, a fixed set of random seeds that “rfe” uses at each resampling iteration was created and submitted to “rfeControl” via the “seeds” parameter. Within each resample, SVM-RFE ranks all predictors according to their individual receiver operating characteristic (ROC) on the 96 training samples. In each iteration less important predictors are removed, the model is fitted to the 96 training samples and the 11 holdout samples are predicted.  The metric to be maximized by “rfe” was set to “Accuracy”.  After initial inspection of the resampling profiles, we noted that accuracy peaked most often between 5-30 predictors. For maximum resolution within these ranges, all subset sizes from 1-40 predictors were tested. Outside of this range, wider intervals were used (45, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500 predictors), resulting in 52 subset sizes in total. For each tested subset within each resample of the external layer an additional “inner layer” of resampling had to be established to determine the tuning parameters of the SVM-model. The details of the inner resampling layer were specified by the helper function “trainControl” and set to three times repeated 10-fold cross-validation (30 resamples). To be precise, each training set from the external layer (96 samples) was partitioned further into 30 internal resamples comprised of 86 “internal training” and 10 “internal holdout” samples, respectively (Supplementary Figure 3). For each value of tuning parameters and each internal resample, the SVMrad model was fit to the 86 internal training samples and the remaining 10 internal holdout samples were predicted. The prediction accuracy from the 30 internal resamples over the different tested hyperparameter values was used to determine the optimal value for the tuning parameters and these parameters were passed to the external layer to fit the model and predict the external hold-outs. SVM-RFE with a radial basis function kernel has two tuning parameters: cost (penalty parameter) and sigma (inverse width of the gaussian kernel). For the cost parameter, the parameter “tuneLength” of the “rfe”- function was set to 20, resulting in cost values ranging from 2-2 - 217. For the sigma parameter an analytical estimate was used which is calculated by “rfe” internally by calling the function “sigest” from the R-package kernlab54. “sigest” uses the methodology proposed by Caputo et al55 to estimate a value for sigma which results in a good prediction performance when used with a radial kernel SVM. We validated this approach initially by a manual search for sigma over a wide range of values (2-15 - 20), but could not find substantially better solutions for our dataset than suggested by “sigest” (data not shown). Hence, using a fixed value for sigma estimated with “sigest” and tuning the SVM over the cost parameter only resulted in a substantially smaller hyperparameter space and reduced computation time for SVM-RFE. To find the best subset size for the entire training set, the prediction accuracy of the external holdout samples for each subset size and each resample was averaged into a resampling profile (Fig. 3b, Fig. 4a), which allowed to determine the best average subset size across all resamples. To generate the final set of predictors, “rfe” repeated the process on the complete training set with the optimal subset size determined from the resampling profile. The performance of the SVM-RFE model was compared to the full model using all predictors using the caret functions “resamples” and “diff”, which compare resampling results of different models on a common data set comprised of identical resamples using a paired t-test56. The resampling-based results for the ten training / test set splits and the final model are tabulated in Supplementary Data 4 tab 1 (P-value_Resampling_full_vs_rfe). The GA procedure was implemented using the function “gafs” and its helper function “gafsControl” from the R package “caret”. The gene expression matrix reduced to the probes found by the preceding SVM-RFE step was used as input into GA. Similarly to the SVM-RFE implementation, SVM-GA was conducted inside an external resampling layer to assess the performance of the GA-model over the generations (external resampling accuracy). 50 external resamples of the training sets or the final dataset were created with the function “createMultiFolds” (k=10, times = 5) and passed to the function “gafsControl”, which controls the outer resampling process of the GA. The computational burden of SVM-GA is higher than for SVM-RFE, so only 50 external resamples were used to complete the analysis in a reasonable amount of time. The prediction performances on the external hold-out samples at each generation across all external resamples were averaged into the external resampling profile (Fig. 3 c, Fig. 4b), which was used to determine the optimal number of iterations the algorithm should proceed (Supplementary Figure 3). To determine the final feature set, “gafs” applied the GA to the entire training set for the optimal number of generations from the resampling process. Further parameters of “gafsControl” were set to enable parallel computing for the external layer, to maximize the test statistic (accuracy) and to use fixed random seeds for reproducibility. In initial runs using the default settings of “gaf” feature reduction was quite inefficient, leading to the removal of only 3-5 predictors on average. For a more effective reduction of feature numbers the size of the initial predictor subsets (chromosomes) in the starting population was reduced. Therefore, the helper function of GA (caretGA$initial) that creates the initial population was modified to produce chromosomes comprised of a random 40% of predictors, instead of creating initial subsets ranging from 10% to 90% of predictors. The GA procedure itself was run for 40 generations, with a population size of 40, a crossover probability of 0.7, a mutation probability of 0.1. Elitism was set to 3, meaning that the best three solutions survive to the next generation. The metric to optimize was set to “accuracy”, the classification method to “svmRadial”. Similarly to the SVM-RFE process, the GA had an additional inner layer of resampling conducted at each generation within each resample and for each chromosome to tune the SVM. The inner resampling layer of GA was set to two times repeated 10-fold cross-validation (20 resamples) by the helper function “trainControl”. For the cost parameter of the SVM the parameter “tuneLength” was set to 12, for cost values between 2-2 – 29. The reduced tune length was chosen to save computation time after it had been determined from the preceding steps that the optimal cost parameter for the SVM was in the range of 2-2 - 27.  For the sigma parameter the estimate computed by “sigest” function from “kernlab” was used as described above. 

## Availability of raw data and R workspace file

The following files can be downloaded [here](https://owncloud.gwdg.de/index.php/s/nHjGbn5sAnU93SV):

*	Feature_Selection.RData
*	[20200322_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median](./20200322_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median)
* AllResults_Accuracy_full_TestSet vs Resampling_IQR0.8.pdf
* AllResults_Accuracy_full_TestSet vs Resampling_IQR0.8_Median.pdf
* AllResults_Accuracy_GA_TestSet vs Resampling_IQR0.8.pdf
* AllResults_Accuracy_GA_TestSet vs Resampling_IQR0.8_median.pdf
* AllResults_Accuracy_rfe_TestSet vs Resampling_IQR0.8.pdf
* AllResults_Accuracy_rfe_TestSet vs Resampling_IQR0.8_median.pdf
* AllResults_rfe_GA40_10TestSets_FinalSet152_IQR0.8.txt
* AllTestSets_AUPRC_SAGA_vs_IVIM_allVectors.pdf
* AllTestSets_AUPRC_SAGA_vs_IVIM_LTR.SF.pdf
* AllTestSets_AUPRC_SAGA_vs_IVIM_nonLTR.SF.pdf
* AllTestSets_AUPRC_SAGA_vs_IVIM152.pdf
* AllTestSets_AUPRC_SAGA_vs_IVIM152_LTR.SF.pdf
* AllTestSets_AUPRC_SAGA_vs_IVIM152_nonLTR.SF.pdf
* AllTestSets_AUROC_SAGA_vs_IVIM152.pdf
* AllTestSets_AUROC_SAGA_vs_IVIM152_LTR.SF.pdf
* AllTestSets_AUROC_SAGA_vs_IVIM152_noLTR.SF.pdf
* AllTestSets_ConfusionMatrix_SAGA.txt
* AllTestSets_LTR.RV_SAGA_vs_IVIM.pdf
* AllTestSets_Other_SAGA_vs_IVIM.pdf
* AllTestSets_ROC_SAGA.pdf
* AllTestSets_ROC_SAGA_RFA_GA_IVIM.pdf
* AllTestSets_ROC_SAGA_vs_IVIM.pdf
* AllTestSets_ROC_SVM.Full.pdf
* AllTestSets_ROC_SVM.GA.pdf
* AllTestSets_ROC_SVM.rfe.pdf
* AllTestSetsConfusionMatrix_FULL.txt
* AllTestSetsConfusionMatrix_rfe.txt
* Annotation.matrix.train.FINAL.txt
* Annotation_SAGA_FINAL_KNOWN_20181128.txt
* ConfusionMatrix_IVIM_all 502 samples_.txt
* DeLongTest_Compound vs IVIM.txt
* DeLongTest_GA vs IVIM.txt
* DeLongTest_RFE vs IVIM.txt
* DeLongTest_SAGA vs IVIM.txt
* ESET_RMA_COMBAT_KNOWN_FullSagaSet152_FINAL.txt
* FINAL_VariableImportance_ROC_IQR0.8.txt
* FinalSet152_Accuracy over FeatureNumber.pdf
* FinalSet152_Accuracy over FeatureNumber_Zoom.pdf
* FinalSet152_Accuracy over Iteration_GA.pdf
* FinalSet152_Accuracy over Iteration_GA_external.pdf
* FinalSet152_ModelDifferences_rfe vs full.txt
* FinalSet152_optVars_GeneticAlgorithm.txt
* FinalSet152_optVars_rfe.txt
* FinalSet152_PCA_11randomVars.pdf
* FinalSet152_PCA_allVars.pdf
* FinalSet152_PCA_allVars_V2.pdf
* FinalSet152_PCA_optVars_GA.pdf
* FinalSet152_PCA_optVars_GA_V2.pdf
* FinalSet152_Resamples_rfe vs full.txt
* FinalSet152_Results_rfe.txt
* IVIM_MTT_ROC.txt
* IVIM_Total_ConfusionMatrix.txt
* SAGA_Targets_FINAL_152.txt
* TestSet01_ConfusionMatrix_allVars.txt
* TestSet01_ConfusionMatrix_optVars.txt
* TestSet01_PCA_allVars.pdf
* TestSet01_PCA_optVars_rfe.pdf
* TestSet01_Predictions_allVars.txt
* TestSet01_Predictions_rfe.txt
* TestSet01_PRROC_optVars_full.pdf
* TestSet01_PRROC_optVars_rfe.pdf
* TestSet01_ROC_allVars.pdf
* TestSet01_ROC_optVars_rfe.pdf
* TestSet010_PRROC_allVars.pdf
* TestSet010_PRROC_optVars_GA.pdf
* TestSet010_PRROC_optVars_rfe.pdf
* TestSet010_ROC_allVars.pdf
* TestSet010_ROC_optVars_GA.pdf
* TestSet010_ROC_optVars_rfe.pdf
* TestSet02_ConfusionMatrix_FULL.txt
* TestSet02_ConfusionMatrix_GA.txt
* TestSet02_ConfusionMatrix_rfe.txt
* TestSet02_Predictions_allVars.txt
* TestSet02_Predictions_optVars_GA.txt
* TestSet02_Predictions_optVars_rfe.txt
* TestSet02_PRROC_optVars_allVars.pdf
* TestSet02_PRROC_optVars_GA.pdf
* TestSet02_PRROC_optVars_rfe.pdf
* TestSet02_ROC_allVars.pdf
* TestSet02_ROC_optVars_GA.pdf
* TestSet02_ROC_optVars_rfe.pdf
* TestSet03_ConfusionMatrix_FULL.txt
* TestSet03_ConfusionMatrix_GA.txt
* TestSet03_ConfusionMatrix_rfe.txt
* TestSet03_Predictions_allVars.txt
* TestSet03_Predictions_optVars_GA.txt
* TestSet03_Predictions_optVars_rfe.txt
* TestSet03_PRROC_optVars_allVars.pdf
* TestSet03_PRROC_optVars_GA.pdf
* TestSet03_PRROC_optVars_rfe.pdf
* TestSet03_ROC_allVars.pdf
* TestSet03_ROC_optVars_GA.pdf
* TestSet03_ROC_optVars_rfe.pdf
* TestSet04_ConfusionMatrix_FULL.txt
* TestSet04_ConfusionMatrix_optVars_rfe.txt
* TestSet04_Predictions_allVars.txt
* TestSet04_Predictions_optVars_rfe.txt
* TestSet04_PRROC_allVars.pdf
* TestSet04_PRROC_optVars_rfe.pdf
* TestSet04_ROC_allVars.pdf
* TestSet04_ROC_optVars_rfe.pdf
* TestSet05_ConfusionMatrix_FULL.txt
* TestSet05_ConfusionMatrix_optVars_rfe.txt
* TestSet05_Predictions_allVars.txt
* TestSet05_Predictions_optVars_rfe.txt
* TestSet05_PRROC_allVars.pdf
* TestSet05_PRROC_optVars_rfe.pdf
* TestSet05_ROC_allVars.pdf
* TestSet05_ROC_optVars_rfe.pdf
* TestSet06_ConfusionMatrix_FULL.txt
* TestSet06_ConfusionMatrix_GA.txt
* TestSet06_ConfusionMatrix_rfe.txt
* TestSet06_Predictions_allVars.txt
* TestSet06_Predictions_optVars_GA.txt
* TestSet06_Predictions_optVars_rfe.txt
* TestSet06_PRROC_allVars.pdf
* TestSet06_PRROC_optVars_GA.pdf
* TestSet06_PRROC_optVars_rfe.pdf
* TestSet06_ROC_allVars.pdf
* TestSet06_ROC_optVars_GA.pdf
* TestSet06_ROC_optVars_rfe.pdf
* TestSet07_ConfusionMatrix_FULL.txt
* TestSet07_ConfusionMatrix_GA.txt
* TestSet07_ConfusionMatrix_rfe.txt
* TestSet07_Predictions_allVars.txt
* TestSet07_Predictions_optVars_GA.txt
* TestSet07_Predictions_optVars_rfe.txt
* TestSet07_PRROC_allVars.pdf
* TestSet07_PRROC_optVars_GA.pdf
* TestSet07_PRROC_optVars_rfe.pdf
* TestSet07_ROC_allVars.pdf
* TestSet07_ROC_optVars_GA.pdf
* TestSet07_ROC_optVars_rfe.pdf
* TestSet08_ConfusionMatrix_FULL.txt
* TestSet08_ConfusionMatrix_GA.txt
* TestSet08_ConfusionMatrix_rfe.txt
* TestSet08_Predictions_allVars.txt
* TestSet08_Predictions_optVars_GA.txt
* TestSet08_Predictions_optVars_rfe.txt
* TestSet08_PRROC_optVars_allVars.pdf
* TestSet08_PRROC_optVars_GA.pdf
* TestSet08_PRROC_optVars_rfe.pdf
* TestSet08_ROC_allVars.pdf
* TestSet08_ROC_optVars_GA.pdf
* TestSet08_ROC_optVars_rfe.pdf
* TestSet09_ConfusionMatrix_FULL.txt
* TestSet09_ConfusionMatrix_optVars_rfe.txt
* TestSet09_Predictions_allVars.txt
* TestSet09_Predictions_optVars_rfe.txt
* TestSet09_PRROC_allVars.pdf
* TestSet09_PRROC_optVars_rfe.pdf
* TestSet09_ROC_allVars.pdf
* TestSet09_ROC_optVars_rfe.pdf
* TestSet10_ConfusionMatrix_FULL.txt
* TestSet10_ConfusionMatrix_GA.txt
* TestSet10_ConfusionMatrix_rfe.txt
* TestSet10_Predictions_allVars.txt
* TestSet10_Predictions_optVars_GA.txt
* TestSet10_Predictions_optVars_rfe.txt
* TestSets2367810_ConfusionMatrix_Full.txt
* TestSets2367810_ConfusionMatrix_GA.txt
* TestSets2367810_ConfusionMatrix_RFE.txt
* TestSets2367810_ROC_full_rfe_GA.pdf
* Toplist_optimal predictors_10testsets.txt
* TrainingSet01_Accuracy over FeatureNumber.pdf
* TrainingSet01_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet01_ModelDifferences_rfe vs full.txt
* TrainingSet01_OptVars_rfe.txt
* TrainingSet01_Resamples_rfe vs full.txt
* TrainingSet01_Results_rfe.txt
* TrainingSet010_Accuracy over FeatureNumber.pdf
* TrainingSet010_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet010_Accuracy over Iteration_GA.pdf
* TrainingSet010_Accuracy over Iteration_GA_external.pdf
* TrainingSet02_Accuracy over FeatureNumber.pdf
* TrainingSet02_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet02_Accuracy over Iteration_GA.pdf
* TrainingSet02_Accuracy over Iteration_GA_external.pdf
* TrainingSet02_ModelDifferences_rfe vs full.txt
* TrainingSet02_optVars_GeneticAlgorithm.txt
* TrainingSet02_OptVars_rfe.txt
* TrainingSet02_Resamples_rfe vs full.txt
* TrainingSet02_Results_rfe.txt
* TrainingSet03_Accuracy over FeatureNumber.pdf
* TrainingSet03_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet03_Accuracy over Iteration_GA.pdf
* TrainingSet03_Accuracy over Iteration_GA_external.pdf
* TrainingSet03_ModelDifferences_rfe vs full.txt
* TrainingSet03_optVars_GeneticAlgorithm.txt
* TrainingSet03_optVars_rfe.txt
* TrainingSet03_Resamples_rfe vs full.txt
* TrainingSet03_Results_rfe.txt
* TrainingSet04_Accuracy over FeatureNumber.pdf
* TrainingSet04_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet04_ModelDifferences_rfe_vs_full.txt
* TrainingSet04_optVars_rfe.txt
* TrainingSet04_Resamples_rfe_vs_full.txt
* TrainingSet04_Results_rfe.txt
* TrainingSet05_Accuracy over FeatureNumber.pdf
* TrainingSet05_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet05_ModelDifferences_rfe_vs_full.txt
* TrainingSet05_optVars_rfe.txt
* TrainingSet05_resamples_rfe_vs_full.txt
* TrainingSet05_Results_rfe.txt
* TrainingSet06_Accuracy over FeatureNumber.pdf
* TrainingSet06_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet06_Accuracy over Iteration_GA.pdf
* TrainingSet06_Accuracy over Iteration_GA_external.pdf
* TrainingSet06_ModelDifferences_rfe_vs_full.txt
* TrainingSet06_optVars_GeneticAlgorithm.txt
* TrainingSet06_optVars_rfe.txt
* TrainingSet06_resamples_rfe_vs_full.txt
* TrainingSet06_Results_rfe.txt
* TrainingSet07_Accuracy over FeatureNumber.pdf
* TrainingSet07_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet07_Accuracy over Iteration_GA.pdf
* TrainingSet07_Accuracy over Iteration_GA_external.pdf
* TrainingSet07_ModelDifferences_rfe vs full.txt
* TrainingSet07_optVars_GeneticAlgorithm.txt
* TrainingSet07_optVars_rfe.txt
* TrainingSet07_resamples_rfe vs full.txt
* TrainingSet07_Results_rfe.txt
* TrainingSet08_Accuracy over FeatureNumber.pdf
* TrainingSet08_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet08_Accuracy over Iteration_GA.pdf
* TrainingSet08_Accuracy over Iteration_GA_external.pdf
* TrainingSet08_Modeldifferences_rfe vs full.txt
* TrainingSet08_optVars_GeneticAlgorithm.txt
* TrainingSet08_optVars_rfe.txt
* TrainingSet08_resamples_rfe vs full.txt
* TrainingSet08_Results_rfe.txt
* TrainingSet09_Accuracy over FeatureNumber.pdf
* TrainingSet09_Accuracy over FeatureNumber_Zoom.pdf
* TrainingSet09_ModelDifferences_rfe_vs_full.txt
* TrainingSet09_optVars_rfe.txt
* TrainingSet09_Resamples_rfe_vs_full.txt
* TrainingSet09_Results_rfe.txt
* TrainingSet10_ModelDifferences_rfe vs full.txt
* TrainingSet10_optVars_GeneticAlgorithm.txt
* TrainingSet10_optVars_rfe.txt
* TrainingSet10_Resamples_rfe vs full.txt
* TrainingSet10_Results_rfe.txt
