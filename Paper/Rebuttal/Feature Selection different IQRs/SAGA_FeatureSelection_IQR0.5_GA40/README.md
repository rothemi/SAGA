# Feature Selection with IQR 0.5
In addition to the main [Feature Selection script](../../../../Paper/Feature%20Selection/20190306_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R) where an IQR of 0.8 was used for nonspecific feature prefiltering, we also constructed the SAGA-SVM classifier with an IQR of 0.5 using the script [20190321_rfe_GA40_152_10TestSets_IQR0.5_FINAL.R](./20190321_rfe_GA40_152_10TestSets_IQR0.5_FINAL.R).

## Availability of raw data and R workspace file

The following files can be downloaded [here](https://owncloud.gwdg.de/index.php/s/K7Rr2ZJlMdyRGgd):
*	AllResults_Accuracy_full_TestSet vs Resampling_IQR0.8.pdf
*	AllResults_Accuracy_GA_TestSet vs Resampling_IQR0.5.pdf
*	AllResults_Accuracy_rfe_TestSet vs Resampling_IQR0.5.pdf
*	AllResults_rfe_GA40_10TestSets_FinalSet152_IQR0.5.txt
*	AllResults_rfe_GA40_10TestSets_FinalSet152_IQR0.5.xlsx
*	Annotation_SAGA_FINAL_KNOWN_20181128.txt
*	ESET_RMA_COMBAT_KNOWN_FullSagaSet152_FINAL.txt
*	FINAL_VariableImportance_ROC.txt
*	FinalSet152_Accuracy over FeatureNumber.pdf
*	FinalSet152_ModelDifferences_rfe vs full.txt
*	FinalSet152_PCA_allVars.pdf
*	FinalSet152_Resamples_rfe vs full.txt
*	SAGA_Targets_FINAL_152.txt
*	TestSet01_ConfusionMatrix_allVars.txt
*	TestSet01_ConfusionMatrix_FULL.txt
*	TestSet01_ConfusionMatrix_GA.txt
*	TestSet01_ConfusionMatrix_optVars.txt
*	TestSet01_ConfusionMatrix_rfe.txt
*	TestSet01_Predictions_allVars.txt
*	TestSet01_Predictions_optVars_GA.txt
*	TestSet01_Predictions_optVars_rfe.txt
*	TestSet01_Predictions_rfe.txt
*	TestSet01_ROC_allVars.pdf
*	TestSet01_ROC_optVars_GA.pdf
*	TestSet01_ROC_optVars_rfe.pdf
*	TestSet02_ConfusionMatrix_FULL.txt
*	TestSet02_ConfusionMatrix_optVars_rfe.txt
*	TestSet02_Predictions_allVars.txt
*	TestSet02_Predictions_optVars_rfe.txt
*	TestSet02_ROC_allVars.pdf
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
*	TestSet07_ConfusionMatrix_rfe.txt
*	TestSet07_Predictions_allVars.txt
*	TestSet07_Predictions_optVars_rfe.txt
*	TestSet07_ROC_allVars.pdf
*	TestSet07_ROC_optVars_rfe.pdf
*	TestSet08_ConfusionMatrix_FULL.txt
*	TestSet08_ConfusionMatrix_rfe.txt
*	TestSet08_Predictions_allVars.txt
*	TestSet08_Predictions_optVars_rfe.txt
*	TestSet08_ROC_allVars.pdf
*	TestSet08_ROC_optVars_rfe.pdf
*	TestSet09_ConfusionMatrix_FULL.txt
*	TestSet09_ConfusionMatrix_optVars_rfe.txt
*	TestSet09_Predictions_allVars.txt
*	TestSet09_Predictions_optVars_rfe.txt
*	TestSet09_ROC_allVars.pdf
*	TestSet09_ROC_optVars_rfe.pdf
*	TestSet10_ConfusionMatrix_FULL.txt
*	TestSet10_ConfusionMatrix_rfe.txt
*	TestSet10_Predictions_allVars.txt
*	TestSet10_Predictions_optVars_rfe.txt
*	TestSet10_ROC_allVars.pdf
*	TestSet10_ROC_optVars_rfe.pdf
*	TrainingSet01_Accuracy over FeatureNumber.pdf
*	TrainingSet01_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet01_Accuracy over Iteration_GA.pdf
*	TrainingSet01_ModelDifferences_rfe vs full.txt
*	TrainingSet01_optVars_GeneticAlgorithm.txt
*	TrainingSet01_OptVars_rfe.txt
*	TrainingSet01_Resamples_rfe vs full.txt
*	TrainingSet01_Results_rfe.txt
*	TrainingSet02_Accuracy over FeatureNumber.pdf
*	TrainingSet02_ModelDifferences_rfe vs full.txt
*	TrainingSet02_OptVars_rfe.txt
*	TrainingSet02_Resamples_rfe vs full.txt
*	TrainingSet03_Accuracy over FeatureNumber.pdf
*	TrainingSet03_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet03_Accuracy over Iteration_GA.pdf
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
*	TrainingSet06_ModelDifferences_rfe_vs_full.txt
*	TrainingSet06_optVars_GeneticAlgorithm.txt
*	TrainingSet06_optVars_rfe.txt
*	TrainingSet06_resamples_rfe_vs_full.txt
*	TrainingSet06_Results_rfe.txt
*	TrainingSet07_Accuracy over FeatureNumber.pdf
*	TrainingSet07_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet07_ModelDifferences_rfe vs full.txt
*	TrainingSet07_optVars_rfe.txt
*	TrainingSet07_resamples_rfe vs full.txt
*	TrainingSet07_Results_rfe.txt
*	TrainingSet08_Accuracy over FeatureNumber.pdf
*	TrainingSet08_Modeldifferences_rfe vs full.txt
*	TrainingSet08_optVars_rfe.txt
*	TrainingSet08_resamples_rfe vs full.txt
*	TrainingSet08_Results_rfe.txt
*	TrainingSet09_Accuracy over FeatureNumber.pdf
*	TrainingSet09_Accuracy over FeatureNumber_Zoom.pdf
*	TrainingSet09_ModelDifferences_rfe_vs_full.txt
*	TrainingSet09_optVars_rfe.txt
*	TrainingSet09_Resamples_rfe_vs_full.txt
*	TrainingSet09_Results_rfe.txt
*	TrainingSet10_Accuracy over FeatureNumber.pdf
*	TrainingSet10_ModelDifferences_rfe vs full.txt
*	TrainingSet10_optVars_rfe.txt
*	TrainingSet10_Resamples_rfe vs full.txt
*	TrainingSet10_Results_rfe.txt
