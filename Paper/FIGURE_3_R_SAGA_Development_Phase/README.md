# FIGURE 3 - R SAGA Development Phase

This R-script contains the code used for the development phase of SAGA. 

## Dataset preparation
The starting point was the quantile normalized and batch corrected expression matrix (36,226 annotated probes, 152 samples with known IVIM properties).  Machine learning was implemented using the R package “caret” based on a support vector machine with a radial basis function kernel (method = "svmRadial"). The dataset was split into 10 different training and test sets, with the training sets comprised of 70% of the samples and test sets comprised of 30% of the samples. To allocate samples to the test or training set the caret function “createDataPartition” (p=0.7) was used, which performs stratified sampling based on the class labels to keep the distribution of transforming and nontransforming samples equal between the training and test sets.

## Unsupervised filtering
First, the R-package “genefilter” was used to discard probes with an interquartile range (IQR) of log2-expression values less than 0.8 to exclude probes interrogating genes that were not expressed at all or show only little variation in the dataset. 

## Recursive feature elimination (SVM-RFE)
SVM-RFE was performed on the training sets by using the function “rfe”. For the external resampling layer, which assesses the impact of the selection process on the model performance, twenty times repeated 10-fold cross-validation was used, implemented by the function “createMultiFolds” (Parameters: k=10, times = 20, 200 resamples in total). The 200 resamples were submitted to the helper function “rfeControl”, which controls the details of the external resampling process of the function “rfe”. The feature selection process itself was carried out for each of the 200 resamples separately and computed in parallel by setting the “rfeControl” parameter: “allowParallel = TRUE”. A fixed set of random seeds that “rfe” uses at each resampling iteration was created and submitted to “rfeControl” via the “seeds” parameter. The metric to be maximized by “rfe” was set to “Accuracy”.  52 predictors subset sizes in total were tested (1…40, 45, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500 predictors).  For each tested subset within each resample of the external layer an additional “inner layer” of resampling had to be established to determine the tuning parameters of the SVM. The details of the inner resampling layer were specified by the helper function “trainControl” and set to three times repeated 10-fold cross-validation (30 resamples). SVM-RFE with a radial basis function kernel has two tuning parameters: cost (penalty parameter) and sigma (inverse width of the gaussian kernel). For the cost parameter, the parameter “tuneLength” of the “rfe”- function was set to 20, resulting in cost values ranging from 2-2 - 217. For the sigma parameter an analytical estimate was used which is calculated by “rfe” internally by calling the function “sigest” from the R-package kernlab. 

## Genetic Algorithm
The GA procedure was implemented using the function “gafs” and its helper function “gafsControl” from the R package “caret”. The gene expression matrix reduced to the probes found by the preceding SVM-RFE step was used as input into GA. 50 external resamples of the training sets were created with the function “createMultiFolds” (k=10, times = 5) and passed to the function “gafsControl”, which controls the outer resampling process of the GA. Further parameters of “gafsControl” were set to enable parallel computing for the external layer, to maximize the test statistic (accuracy) and to use fixed random seeds for reproducibility. For a more effective reduction of feature numbers, the size of the initial predictor subsets (chromosomes) in the starting population was reduced. Therefore, the helper function of GA (caretGA$initial) that creates the initial population was modified to produce chromosomes comprised of a random 40% of predictors. The GA procedure was run for 40 generations, with a population size of 40, a crossover probability of 0.7, a mutation probability of 0.1. Elitism was set to 3, meaning that the best three solutions survive to the next generation. The metric to optimize was set to “accuracy”, the classification method to “svmRadial”. Similarly to the SVM-RFE process, the GA had an additional inner layer of resampling conducted at each generation within each resample and for each chromosome to tune the SVM. The inner resampling layer of GA was set to two times repeated 10-fold cross-validation (20 resamples) by the helper function “trainControl”. 

## Prediction of test sets
Samples in the hold-out test sets were predicted after training a support vector machine with radial kernel on the training set using all predictors (full model) or reduced to the optimal predictors found by SVM-RFE and SVM-GA (reduced models) by using the caret functions “train” and “predict”, respectively. For training the full and the reduced SVM-models, identical parameters and resamples were specified in the “train” function (method = "svmRadial", metric = "Accuracy", tuneLength = 20, twenty repeats of 10-fold cross-validation). The function “predict” was used with the parameter “type” set to “prob”, which computes the probability that a sample belongs to a given class. An unknown sample was considered belonging to the class “transforming” when the probability for class “transforming” was greater than 0.5. Performance estimates (sensitivity, specificity, accuracy, kappa) for the predicted test sets were computed using the function “confusionMatrix” on the predicted and the true class labels, respectively. The “pROC” R-package63 (v1.15.3) was used to compute and visualize the ROC curves for the test sets using the function “roc” on the probability for class “transforming” as output by the “predict” function. P values to compare the difference between the AUROC of two unpaired ROC curves were performed with the “roc.test” function using the “delong” method and the alternative hypothesis set to “greater”. Precision recall curves were generated using the R-package “PRROC” (v.1.3.1). 


## Availability of R script and workspace file

*	[FIGURE_3_R_SAGA_Development_Phase.RData](https://www.dropbox.com/s/eauuqpgr2bfnl9v/Feature_Selection.RData?dl=0)
*	[20200322_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R](./20200322_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R)