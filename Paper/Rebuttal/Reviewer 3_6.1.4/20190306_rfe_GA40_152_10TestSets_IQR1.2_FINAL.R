#############################################################################################################################################
#################################### Construction of SAGA-SVM Classifier IQR 1.2 ############################################################
#############################################################################################################################################
library(limma)
library(genefilter)
library(Rtsne)
library(RColorBrewer)
library(caret)
library(kernlab)
library(ggplot2)
library(dplyr)
library(pROC)
library(doMC)
registerDoMC(cores = 10)    

#############################################################################################################################################
#### Read in Data / Normalized and Batch corrected log2-ExpressionSet of 152 SAGA Samples ###################################################
#############################################################################################################################################

pData      <- read.delim("SAGA_Targets_FINAL_152.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation <- read.delim("Annotation_SAGA_FINAL_KNOWN_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
eset.batch <- as.matrix(read.delim("ESET_RMA_COMBAT_KNOWN_FullSagaSet152_FINAL.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE,row.names = 1))

#### visualize ESET
set.seed(476)  
tsne_out <- Rtsne(t(eset.batch),dims = 2, initial_dims = 20, perplexity = 11,
                  theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
                  verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData$Design_Color, pch=16, cex=1.3)

#############################################################################################################################################
####  Set global Parameters for GENETIC ALGORITHM ###########################################################################################
#############################################################################################################################################
svmGA <- caretGA  # predefined helper functions for the genetic algorithm 

# define function to create an initial population with individuals that consist of 40% of the features on average
initial40 <- function (vars, popSize, ...) {x <- matrix(NA, nrow = popSize, ncol = vars)
probs <- rep(0.6,length = popSize)   
for (i in 1:popSize) {
  x[i, ] <- sample(0:1, replace = TRUE, size = vars, prob = c(probs[i],1 - probs[i]))
}
var_count <- apply(x, 1, sum)
if (any(var_count == 0)) {
  for (i in which(var_count == 0)) {
    x[i, ] <- sample(0:1, replace = TRUE, size = vars)
  }
}
x
}
environment(initial40) <- asNamespace('caret')
svmGA$initial <- initial40

# set all seeds for running the genetic algorithm in parallel over the 50 different resamples ##
set.seed(123)
seeds.GA <- vector(mode = "integer", length = length(index.GA.2)+1)    # B+1 elements where B is the number of resamples = 51
for(i in 1:length(index.GA.2)+1) seeds.GA[[i]] <- sample.int(10000, 1)

#############################################################################################################################################
#############################################################################################################################################
#### I. Split 01 / TestSet 01 ###############################################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 01 and TestSet 01 ############################################################
##############################################################################################################
set.seed(123)
split.1        <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.1 <- eset.batch[, split.1] 
matrix.test.1  <- eset.batch[,-split.1] 
pData.train.1  <- pData[split.1,]
pData.test.1   <- pData[-split.1,]

table(pData.train.1$Design)
table(pData.test.1$Design)

##############################################################################################################
#### 2. nonspecific feature prefiltering  ####################################################################
##############################################################################################################
f1       <- function(x) (IQR(x) > 0.8)    
fselect.1  <- genefilter(matrix.train.1, filterfun(f1))
summary(fselect.1)                            # 1151 genes selected with interquartile range of log2-int >0.8
matrix.train.1 <-matrix.train.1[fselect.1,]

##############################################################################################################
#### 3. SVM: FULL MODEL ######################################################################################
##############################################################################################################
matrix.train.1 <- (t(matrix.train.1))
labels.train.1 <- as.factor(pData.train.1$Class)

#### 3.1 SetUp SVM with radial kernel ########################################################################
##############################################################################################################
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))   # calculates Accuracy, Sens, Spec, ROC, kappa of external resamples

## create 200 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for model comparison
set.seed(123)
index.1 <- createMultiFolds(labels.train.1, k=10, times = 20)  

fullCtrl.1 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.1,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.1 <- train(matrix.train.1,labels.train.1,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.1)

svmFull.1  

# Resampling results for the FullModel with 107 samples 1151 predictors: 
# final tuning parameters: sigma = 0.0006816859, C = 0.5 
#    ROC       Sens   Spec       Accuracy   Kappa
# 0.9508452  0.84575  0.9377381  0.8978712  0.7885724


##############################################################################################################
#### 4. SVM-RFE on TrainingSet 01 ############################################################################
##############################################################################################################

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

## set number of features to test (subsetSizes)
FeatureNumbers <- c(seq(1,40,by=1),45,50,60,70,80,90,100,200,300,400,500)                 # 52 subset sizes to test incl. the full set 

## set all seeds for reproducibilty when running models in parallel 
set.seed(123)
seeds.rfe <- vector(mode = "list", length = length(index.1)+1)                            # 52 seeds for each of the 200 resamples + 1 for the complete dataSet
for(i in 1:length(index.1)) seeds.rfe[[i]] <- sample.int(10000, length(FeatureNumbers)+1) 
seeds.rfe[[length(index.1)+1]] <- sample.int(10000, 1)

## set details of outer resampling process
outerctrl.1      <- rfeControl(  method = "repeatedcv", repeats = 20,  # 10CVn20 = 200 resamples 
                                 index = index.1,                      # use same 200 resamples as in the full model for comparability
                                 seeds = seeds.rfe,                    # use defined seeds for nested CV within each resample for reproducibility
                                 saveDetails = TRUE,
                                 returnResamp="final", 
                                 verbose = TRUE, 
                                 allowParallel = TRUE)                 # parallel on 65 cores (AWS EC2 c5.18xlarge	72 vCPUs	144 Gb RAM) 

outerctrl.1$functions         <- caretFuncs                            
outerctrl.1$functions$summary <- fiveStats                             # calculate Accuracy, ROC, kappa, Sens + Spec for each resample                    

#### 4.2 Parameters for inner (nested) cross-validation loop for hyperparameter tuning within each resample and for each subset size #######
##############################################################################################################
innerctrl <- trainControl(method = "repeatedcv",repeats = 3,           # 10CVn3
                          verboseIter = FALSE,
                          classProbs = TRUE,
                          allowParallel = FALSE)                      

#### 4.3.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning x 20 values for cost parameter = 6,240,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.1  <- rfe(matrix.train.1, labels.train.1, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.1,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.1 # 9 optVars found
write.table(rfe.1$results, file = "TrainingSet01_Results_rfe.txt", sep="\t",col.names=NA)

optFeatures.rfe.1 <- cbind(rfe.1$optVariables, Annotation[rfe.1$optVariables,])
write.table(optFeatures.rfe.1, file = "TrainingSet01_OptVars_rfe.txt", sep="\t",col.names=NA)

# Plot Accuracy over FeatureNumber
trellis.par.set(caretTheme())
plot(rfe.1, type = c("g", "o"))
plot(rfe.1, type = c("g", "o"), xlim = c(0,61))

# Plot ROC over FeatureNumber
plot(rfe.1$results$Variables,rfe.1$results$ROC, xlim = c(0,201),type = "o",panel.first = grid())

                         
#### 4.4.compare resampling performances between full model and rfe ##########################################
##############################################################################################################

rfeResamples.1 <- resamples(list("SVM_full.1" = svmFull.1,"SVM_RFE.1" = rfe.1))
sink("TrainingSet01_Resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.1)
sink()

modelDifferences.1 <- diff(rfeResamples.1)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet01_ModelDifferences_Resamples_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.1)
sink()


##### 5  PCA on optimal variables ##############################################################
################################################################################################
matrix.train.rfe.1     <- t(matrix.train.1[,rfe.1$optVariables])         # subset to the optVars found by rfe
matrix.test.rfe.1      <- matrix.test.1[rfe.1$optVariables,]             # subset to the optVars found by rfe
matrix.opt.1           <- cbind(matrix.train.rfe.1, matrix.test.rfe.1)   # combine training and test set
pData.opt.1            <- rbind(pData.train.1,pData.test.1)
pData.opt.1$Set        <- ifelse(row.names(pData.opt.1)%in% row.names(pData.train.1),"train","test")
pData.opt.1$Design_Color <- ifelse(pData.opt.1$Set == "train", pData.opt.1$Design_Color,"#000000")   # TestSet samples are black

pca.1         <- prcomp(t(matrix.opt.1))           
plot(pca.1$x, pch=16, col=pData.opt.1$Design_Color, cex=1.5, asp=1)
legend(2,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.1$Design_Color), pch=16, bty="n", cex=0.8)

##### compared to PCA on all variables ########################################################
matrix.full.1 <- rbind(matrix.train.1, t(matrix.test.1[colnames(matrix.train.1),]))
pca.full.1    <- prcomp(matrix.full.1)           
plot(pca.full.1$x, pch=16, col=pData.opt.1$Design_Color, cex=1, asp=1)
legend(30,-20, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.1$Design_Color), pch=16, bty="n", cex=0.8)


##### 6  Train SVM on optVar FeatureSet and predict independent TestSet.1 ######################
################################################################################################

set.seed(721)
svmOpt.1  <- train(t(matrix.train.rfe.1),labels.train.1,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.1)   # use same parameters as for the SVMfull model (10CVn20/index.1)
svmOpt.1  

# Resampling results for the optVars Model with 107 samples and 7 predictors 
# final tuning parameters: sigma = 0.3736986, C = 32 

#    ROC       Sens   Spec       Accuracy   Kappa
# 0.9797976  0.95525  0.9582143  0.9571288  0.9124974 

# note: these values are better than the external resampling results (above) due to 
# positive bias since the optimal features were determined for this TestSet

##### Predict TestSet 1 with optimal predictors
Prediction_rfe.1 <- predict(svmOpt.1,t(matrix.test.rfe.1), type = "prob")
Prediction_rfe.1$Prediction_rfe.1 <- ifelse(Prediction_rfe.1$transforming>0.50,"transforming","untransforming")
Prediction_rfe.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_rfe.1)
write.table(Prediction_rfe.1, file = paste("TestSet01_Predictions_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet 1 with all 1151 predictors
matrix.test.1.full <- matrix.test.1[row.names(t(matrix.train.1)),]
Prediction_SVM_full.1 <- predict(svmFull.1,t(matrix.test.1.full), type = "prob")
Prediction_SVM_full.1$Prediction_SVM_full <- ifelse(Prediction_SVM_full.1$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_SVM_full.1)
write.table(Prediction_SVM_full.1, file = paste("TestSet01_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 7 Performance of opVars Classifier on TestSet 1 ########################################################
#############################################################################################################

#### 7.1 Confusion matrix  ##################################################################################  
sink("TestSet01_ConfusionMatrix_optVars_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.1$Prediction_rfe.1), as.factor(Prediction_rfe.1$TrueLabel))
sink()

sink("TestSet01_ConfusionMatrix_allVars.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.1$Prediction_SVM_full), as.factor(Prediction_SVM_full.1$TrueLabel))
sink()


#### 7.2.ROC on probability "transforming" on TestSet 1 #####################################################

Prediction_rfe.1$Class <- as.factor(ifelse(Prediction_rfe.1$TrueLabel == "transforming","transforming","nontransforming"))

roc1 <- roc(Prediction_rfe.1$Class,                    # response vector (factor or character)
            Prediction_rfe.1$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=0.5)


Prediction_SVM_full.1$Class <- as.factor(ifelse(Prediction_SVM_full.1$TrueLabel == "transforming","transforming","nontransforming"))
roc2 <- roc(Prediction_SVM_full.1$Class,                    # response vector (factor or character)
            Prediction_SVM_full.1$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### II. Split 02 / TestSet 02 ##############################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(1440)
split.2 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.2 <- eset.batch[, split.2] 
matrix.test.2  <- eset.batch[,-split.2] 
pData.train.2  <- pData[split.2,]
pData.test.2   <- pData[-split.2,]

table(pData.train.2$Design)
table(pData.test.2$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering TrainingSet 02 ######################################################
##############################################################################################################

fselect.2  <- genefilter(matrix.train.2, filterfun(f1))
summary(fselect.2)
matrix.train.2 <-matrix.train.2[fselect.2,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 02 #######################################################################
##############################################################################################################

matrix.train.2 <- (t(matrix.train.2))
labels.train.2 <- as.factor(pData.train.2$Class)


#### 2.1 SetUp SVMrad for full model #########################################################################
##############################################################################################################

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.2 <- createMultiFolds(labels.train.2, k=10, times = 20)  

fullCtrl.2 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.2,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.2 <- train(matrix.train.2,labels.train.2,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.2)

svmFull.2  

# Resampling results for the FullModel with 107 samples 1201 predictors: 
# final tuning parameters: sigma = 0.000713516, C = 64 
#    ROC       Sens   Spec       Accuracy   Kappa
# 0.9688571  0.86125  0.9110714  0.8902045  0.7739778


##############################################################################################################
#### 3. SVM-RFE TrainingSet 02 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.2      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               index = index.2,
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.2$functions         <- caretFuncs
outerctrl.2$functions$summary <- fiveStats


#### 3.2.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.2 <- rfe(matrix.train.2, labels.train.2, 
                         sizes=FeatureNumbers,
                         rfeControl=outerctrl.2,
                         metric = "Accuracy",
                         ## Options to train()
                         method="svmRadial",
                         tuneLength = 20,
                         trControl = innerctrl))

rfe.2 # 6 optVars found
write.table(rfe.2$results, file = "TrainingSet02_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.2, type = c("g", "o"))
plot(rfe.2, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.2 <- cbind(rfe.2$optVariables, Annotation[rfe.2$optVariables,])
write.table(optFeatures.rfe.2, file = "TrainingSet02_OptVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3.compare external resampling performances ############################################################
##############################################################################################################

rfeResamples.2 <- resamples(list("SVM_full.2" = svmFull.2,"SVM_RFE.2" = rfe.2))
sink("TrainingSet02_Resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.2)
sink()

modelDifferences.2 <- diff(rfeResamples.2)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet02_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.2)
sink()

##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.2     <- t(matrix.train.2[,rfe.2$optVariables])
matrix.test.rfe.2      <- matrix.test.2[rfe.2$optVariables,]
matrix.opt.2           <- cbind(matrix.train.rfe.2, matrix.test.rfe.2)
pData.opt.2            <- rbind(pData.train.2,pData.test.2)
pData.opt.2$Set        <- ifelse(row.names(pData.opt.2)%in% row.names(pData.train.2),"train","test")
pData.opt.2$Design_Color <- ifelse(pData.opt.2$Set == "train", pData.opt.2$Design_Color,"#000000")

pca.2         <- prcomp(t(matrix.opt.2))           
plot(pca.2$x, pch=16, col=pData.opt.2$Design_Color, cex=1, asp=1)
legend(4,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.2$Design_Color), pch=16, bty="n", cex=0.8)


##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and predict independent TestSet.2 ##############################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.train.rfe.2  <- matrix.train.2[,rfe.2$optVariables]
matrix.test.rfe.2  <- t(matrix.test.2[rfe.2$optVariables,])

# train SVM and predict TestSet02 on 35 optVars found by rfe #################################################
set.seed(721)
svmOpt.rfe.2  <- train(matrix.train.rfe.2,labels.train.2,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.2)
svmOpt.rfe.2  

Prediction_rfe.2 <- predict(svmOpt.rfe.2,matrix.test.rfe.2, type = "prob")
Prediction_rfe.2$Prediction_rfe.2 <- ifelse(Prediction_rfe.2$transforming>0.50,"transforming","untransforming")
Prediction_rfe.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_rfe.2)
write.table(Prediction_rfe.2, file = paste("TestSet02_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet02 with all predictors ##################################################################
matrix.test.2.full <- matrix.test.2[row.names(t(matrix.train.2)),]
Prediction_SVM_full.2 <- predict(svmFull.2,t(matrix.test.2.full), type = "prob")
Prediction_SVM_full.2$Prediction_SVM_full <- ifelse(Prediction_SVM_full.2$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_SVM_full.2)
write.table(Prediction_SVM_full.2, file = paste("TestSet02_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 02 #####################################################
#############################################################################################################

#### 7.1 Confusion matrices  ################################################################################  
#############################################################################################################
sink("TestSet02_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.2$Prediction_rfe.2), as.factor(Prediction_rfe.2$TrueLabel))
sink()

sink("TestSet02_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.2$Prediction_SVM_full), as.factor(Prediction_SVM_full.2$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 02 #######################################################
#############################################################################################################
Prediction_rfe.2$Class <- as.factor(ifelse(Prediction_rfe.2$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.2$Class,                    
               Prediction_rfe.2$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_SVM_full.2$Class <- as.factor(ifelse(Prediction_SVM_full.2$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.2$Class,                    # response vector (factor or character)
               Prediction_SVM_full.2$transforming,             # predictor vector (numeric)
              percent=TRUE, levels=c("nontransforming","transforming"),
              plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
              print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### III. Split 03 / TestSet 03 #############################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(876)
split.3 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.3 <- eset.batch[, split.3] 
matrix.test.3  <- eset.batch[,-split.3] 
pData.train.3  <- pData[split.3,]
pData.test.3   <- pData[-split.3,]

table(pData.train.3$Design)
table(pData.test.3$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering in TrainingSet 03 ###################################################
##############################################################################################################
fselect.3  <- genefilter(matrix.train.3, filterfun(f1))
summary(fselect.3)
matrix.train.3 <-matrix.train.3[fselect.3,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 03 #######################################################################
##############################################################################################################

matrix.train.3 <- (t(matrix.train.3))
labels.train.3 <- as.factor(pData.train.3$Class)

#### 2.1 SetUp SVMrad for full model #########################################################################
##############################################################################################################

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(12345)
index.3 <- createMultiFolds(labels.train.3, k=10, times = 20)  

fullCtrl.3 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.3,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.3 <- train(matrix.train.3,labels.train.3,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.3)

svmFull.3  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 03 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.3      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               index = index.3,
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.3$functions         <- caretFuncs
outerctrl.3$functions$summary <- fiveStats


#### 3.2.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.3  <- rfe(matrix.train.3, labels.train.3, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.3,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.3    # 29 predictors found 
write.table(rfe.3$results, file = "TrainingSet03_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.3, type = c("g", "o"))
plot(rfe.3, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.3 <- cbind(rfe.3$optVariables, Annotation[rfe.3$optVariables,])
write.table(optFeatures.rfe.3, file = "TrainingSet03_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3.compare resampling performances #####################################################################
##############################################################################################################

rfeResamples.3 <- resamples(list("SVM_full.3" = svmFull.3,"SVM_RFE.3" = rfe.3))
sink("TrainingSet03_Resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.3)
sink()

modelDifferences.3 <- diff(rfeResamples.3)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet03_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.3)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION ########################################
################################################################################################
matrix.train.rfe.3 <- matrix.train.3[,rfe.3$optVariables]   # subset fot the 36 optVars

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2210)
index.GA.3 <- createMultiFolds(labels.train.3, k=10, times = 5)  

outerctrl.GA.3 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.3,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                         external = TRUE),
                              allowParallel = TRUE)                                  

#### 4.2 run GA  ###############################################################################
################################################################################################

system.time(GA.3<- gafs(matrix.train.rfe.3, labels.train.3, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.3,
                        ## Now we pass options to `train` via "svmGA":               
                        metric = "Accuracy",
                        method = "svmRadial",
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.3 analyze results from Genetic Algorithm #################################################
################################################################################################

GA.3
optVars.GA.3 <- Annotation[GA.3$optVariables,]
write.table(optVars.GA.3, file = "TrainingSet03_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.3 <- GA.3$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.3 <- arrange(performance.external.3, Iter)
performance.3          <- GA.3$ga$internal               # export internal accuracy (within the resample)
performance.3$AccuracyExternal <- aggregate(performance.external.3$Accuracy, by=list(Iter=performance.external.3$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.3 <- data.frame(Iter = c(performance.3$Iter,performance.3$Iter), Accuracy = c(performance.3$AccuracyExternal,performance.3$Accuracy), Group=c(rep("external",40),rep("internal",40)))

## plot internal and external accuracy over the iterations 
ggplot(performance.long.3, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.5,se = T) +
  theme_bw()

##############################################################################################################
##### 5  PCA on optimal variables TestSet 03 #################################################################
##############################################################################################################
matrix.train.rfe.3     <- t(matrix.train.3[,GA.3$optVariables])
matrix.test.rfe.3      <- matrix.test.3[GA.3$optVariables,]
matrix.opt.3           <- cbind(matrix.train.rfe.3, matrix.test.rfe.3)
pData.opt.3            <- rbind(pData.train.3,pData.test.3)
pData.opt.3$Set        <- ifelse(row.names(pData.opt.3)%in% row.names(pData.train.3),"train","test")
pData.opt.3$Design_Color <- ifelse(pData.opt.3$Set == "train", pData.opt.3$Design_Color,"#000000")

pca.3         <- prcomp(t(matrix.opt.3))           
plot(pca.3$x, pch=16, col=pData.opt.3$Design_Color, cex=1, asp=1)
legend(4,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.3$Design_Color), pch=16, bty="n", cex=0.8)

# clean up
rm(matrix.train.rfe.3,matrix.test.rfe.3,matrix.opt.3,pData.opt.3)

##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet03 #######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.test.rfe.3  <- t(matrix.test.3[rfe.3$optVariables,])
matrix.test.GA.3   <- t(matrix.test.3[GA.3$optVariables,])
matrix.train.GA.3  <- matrix.train.3[,GA.3$optVariables]

# train SVM and predict TestSet03 on 36 optVars found by rfe
set.seed(721)
svmOpt.rfe.3  <- train(matrix.train.rfe.3,labels.train.3,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.3)
svmOpt.rfe.3  

Prediction_rfe.3 <- predict(svmOpt.rfe.3,matrix.test.rfe.3, type = "prob")
Prediction_rfe.3$Prediction_rfe.3 <- ifelse(Prediction_rfe.3$transforming>0.50,"transforming","untransforming")
Prediction_rfe.3 <- cbind(pData.test.3[,c(1:3)],TrueLabel=pData.test.3$Class,Prediction_rfe.3)
write.table(Prediction_rfe.3, file = paste("TestSet03_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

# train SVM and predict TestSet03 on 16 optVars found by rfe-GA
set.seed(721)
svmOpt.GA.3  <- train(matrix.train.GA.3,labels.train.3,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.3)
svmOpt.GA.3  

Prediction_GA.3 <- predict(svmOpt.GA.3,matrix.test.GA.3, type = "prob")
Prediction_GA.3$Prediction_GA.3 <- ifelse(Prediction_GA.3$transforming>0.50,"transforming","untransforming")
Prediction_GA.3 <- cbind(pData.test.3[,c(1:3)],TrueLabel=pData.test.3$Class,Prediction_GA.3)
write.table(Prediction_GA.3, file = paste("TestSet03_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


##### Predict TestSet03 with all predictors
matrix.test.3.full <- matrix.test.3[row.names(t(matrix.train.3)),]
Prediction_SVM_full.3 <- predict(svmFull.3,t(matrix.test.3.full), type = "prob")
Prediction_SVM_full.3$Prediction_SVM_full <- ifelse(Prediction_SVM_full.3$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.3 <- cbind(pData.test.3[,c(1:3)],TrueLabel=pData.test.3$Class,Prediction_SVM_full.3)
write.table(Prediction_SVM_full.3, file = paste("TestSet03_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 03 #####################################################
#############################################################################################################

#### 7.1 Confusion matrices  ################################################################################  
#############################################################################################################
sink("TestSet03_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.3$Prediction_rfe.3), as.factor(Prediction_rfe.3$TrueLabel))
sink()

sink("TestSet03_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.3$Prediction_GA.3), as.factor(Prediction_GA.3$TrueLabel))
sink()

sink("TestSet03_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.3$Prediction_SVM_full), as.factor(Prediction_SVM_full.3$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 03 ########################################################
##############################################################################################################
Prediction_rfe.3$Class <- as.factor(ifelse(Prediction_rfe.3$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.3$Class,                    
               Prediction_rfe.3$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_GA.3$Class <- as.factor(ifelse(Prediction_GA.3$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA <- roc(Prediction_GA.3$Class,                    # response vector (factor or character)
              Prediction_GA.3$transforming,             # predictor vector (numeric)
              percent=TRUE, levels=c("nontransforming","transforming"),
              plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
              print.auc=T,print.thres=0.5)

Prediction_SVM_full.3$Class <- as.factor(ifelse(Prediction_SVM_full.3$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.3$Class,                    # response vector (factor or character)
               Prediction_SVM_full.3$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### IV. Split 04 / TestSet 04 ##############################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(345)
split.4 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.4 <- eset.batch[, split.4] 
matrix.test.4  <- eset.batch[,-split.4] 
pData.train.4  <- pData[split.4,]
pData.test.4   <- pData[-split.4,]

table(pData.train.4$Design)
table(pData.test.4$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering on TrainingSet ######################################################
##############################################################################################################
fselect.4  <- genefilter(matrix.train.4, filterfun(f1))
summary(fselect.4)
matrix.train.4 <-matrix.train.4[fselect.4,]

##############################################################################################################
#### 2. SVM: FULL MODEL on TrainingSet #######################################################################
##############################################################################################################

matrix.train.4 <- (t(matrix.train.4))
labels.train.4 <- as.factor(pData.train.4$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(999)
index.4 <- createMultiFolds(labels.train.4, k=10, times = 20)  

fullCtrl.4 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.4,
                           allowParallel = TRUE)

set.seed(721)
svmFull.4 <- train(matrix.train.4,labels.train.4,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.4)

svmFull.4  

##############################################################################################################
#### 3. SVM-RFE on TrainingSet 04 ############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.4      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.4,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.4$functions         <- caretFuncs
outerctrl.4$functions$summary <- fiveStats


#### 3.2.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.4  <- rfe(matrix.train.4, labels.train.4, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.4,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.4    # 6 predictors chosen
write.table(rfe.4$results, file = "TrainingSet04_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.4, type = c("g", "o"))
plot(rfe.4, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.4 <- cbind(rfe.4$optVariables, Annotation[rfe.4$optVariables,])
write.table(optFeatures.rfe.4, file = "TrainingSet04_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances #####################################################################
##############################################################################################################

rfeResamples.4 <- resamples(list("SVM_full.4" = svmFull.4,"SVM_RFE.4" = rfe.4))
sink("TrainingSet04_Resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.4)
sink()

modelDifferences.4 <- diff(rfeResamples.4)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet04_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.4)
sink()


##### 4  PCA on optimal variables TestSet 04 #################################################################
##############################################################################################################
matrix.train.rfe.4     <- t(matrix.train.4[,rfe.4$optVariables]) # subset train matrix on the optimal predictors 
matrix.test.rfe.4      <- matrix.test.4[rfe.4$optVariables,]     # subset test matrix on the optimal predictors 
matrix.opt.4           <- cbind(matrix.train.rfe.4, matrix.test.rfe.4)
pData.opt.4            <- rbind(pData.train.4,pData.test.4)
pData.opt.4$Set        <- ifelse(row.names(pData.opt.4)%in% row.names(pData.train.4),"train","test")
pData.opt.4$Design_Color <- ifelse(pData.opt.4$Set == "train", pData.opt.4$Design_Color,"#000000")

pca.4         <- prcomp(t(matrix.opt.4))           
plot(pca.4$x, pch=16, col=pData.opt.4$Design_Color, cex=1.5, asp=1)
legend(2,-1, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.4$Design_Color), pch=16, bty="n", cex=0.8)


##### 5  Train SVM on optVar FeatureSet.4 and predict independent TestSet 04 #################################
##############################################################################################################

set.seed(721)
svmOpt.4  <- train(t(matrix.train.rfe.4),labels.train.4,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.4)
svmOpt.4  

Prediction_rfe.4 <- predict(svmOpt.4,t(matrix.test.rfe.4), type = "prob")
Prediction_rfe.4$Prediction_rfe.4 <- ifelse(Prediction_rfe.4$transforming>0.50,"transforming","untransforming")
Prediction_rfe.4 <- cbind(pData.test.4[,c(1:3)],TrueLabel=pData.test.4$Class,Prediction_rfe.4)
write.table(Prediction_rfe.4, file = paste("TestSet04_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet04 with all predictors #################################################################
matrix.test.4.full <- matrix.test.4[row.names(t(matrix.train.4)),]
Prediction_SVM_full.4 <- predict(svmFull.4,t(matrix.test.4.full), type = "prob")
Prediction_SVM_full.4$Prediction_SVM_full <- ifelse(Prediction_SVM_full.4$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.4 <- cbind(pData.test.4[,c(1:3)],TrueLabel=pData.test.4$Class,Prediction_SVM_full.4)
write.table(Prediction_SVM_full.4, file = paste("TestSet04_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 6 Performance of opVars Classifier on test samples TestSet 04 ##########################################
#############################################################################################################

#### 6.1 Confusion matrix  ##################################################################################  
sink("TestSet04_ConfusionMatrix_optVars_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.4$Prediction_rfe.4), as.factor(Prediction_rfe.4$TrueLabel))
sink()

sink("TestSet04_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.4$Prediction_SVM_full), as.factor(Prediction_SVM_full.4$TrueLabel))
sink()


#### 6.2 ROC of TestSet 04 #################################################################################

Prediction_rfe.4$Class <- as.factor(ifelse(Prediction_rfe.4$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.4$Class,                    
               Prediction_rfe.4$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_SVM_full.4$Class <- as.factor(ifelse(Prediction_SVM_full.4$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.4$Class,                    # response vector (factor or character)
               Prediction_SVM_full.4$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### Split 05 / TestSet 05 ##################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(1431)
split.5 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.5 <- eset.batch[, split.5] 
matrix.test.5  <- eset.batch[,-split.5] 
pData.train.5  <- pData[split.5,]
pData.test.5   <- pData[-split.5,]

table(pData.train.5$Design)
table(pData.test.5$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering on TrainingSet 05 ###################################################
##############################################################################################################
fselect.5  <- genefilter(matrix.train.5, filterfun(f1))
summary(fselect.5)
matrix.train.5 <-matrix.train.5[fselect.5,]

##############################################################################################################
#### 2. SVM: FULL MODEL on TrainingSet 05 ####################################################################
##############################################################################################################

matrix.train.5 <- (t(matrix.train.5))
labels.train.5 <- as.factor(pData.train.5$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.5 <- createMultiFolds(labels.train.5, k=10, times = 20)  

fullCtrl.5 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.5,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.5 <- train(matrix.train.5,labels.train.5,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.5)

svmFull.5  

##############################################################################################################
#### 3. SVM-RFE on TrainingSet 05 ############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.5      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.5,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.5$functions         <- caretFuncs
outerctrl.5$functions$summary <- fiveStats


#### 3.2 SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.5  <- rfe(matrix.train.5, labels.train.5, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.5,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.5    # 9 optimal predictors found 
write.table(rfe.5$results, file = "TrainingSet05_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.5, type = c("g", "o"))
plot(rfe.5, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.5 <- cbind(rfe.5$optVariables, Annotation[rfe.5$optVariables,])
write.table(optFeatures.rfe.5, file = "TrainingSet05_optVars_rfe_Split7S1431_IQR1.2_CVn20S1234.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances on TrainingSet 05 ###################################################
##############################################################################################################

rfeResamples.5 <- resamples(list("SVM_full.5" = svmFull.5,"SVM_RFE.5" = rfe.5))
sink("TrainingSet05_resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.5)
sink()

modelDifferences.5 <- diff(rfeResamples.5)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet05_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.5)
sink()

##############################################################################################################
##### 4  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.5     <- t(matrix.train.5[,rfe.5$optVariables])
matrix.test.rfe.5      <- matrix.test.5[rfe.5$optVariables,]
matrix.opt.5           <- cbind(matrix.train.rfe.5, matrix.test.rfe.5)
pData.opt.5            <- rbind(pData.train.5,pData.test.5)
pData.opt.5$Set        <- ifelse(row.names(pData.opt.5)%in% row.names(pData.train.5),"train","test")
pData.opt.5$Design_Color <- ifelse(pData.opt.5$Set == "train", pData.opt.5$Design_Color,"#000000")

pca.5         <- prcomp(t(matrix.opt.5))           
plot(pca.5$x, pch=16, col=pData.opt.5$Design_Color, cex=1.5, asp=1)
legend(-3.8,1.5, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.5$Design_Color), pch=16, bty="n", cex=1)

##############################################################################################################
##### 5 Train SVM on optVars from rfe and predict independent TestSet 05 #####################################
##############################################################################################################

set.seed(721)
svmOpt.5  <- train(t(matrix.train.rfe.5),labels.train.5,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.5)
svmOpt.5  

Prediction_rfe.5 <- predict(svmOpt.5,t(matrix.test.rfe.5), type = "prob")
Prediction_rfe.5$Prediction_rfe.5 <- ifelse(Prediction_rfe.5$transforming>0.50,"transforming","untransforming")
Prediction_rfe.5 <- cbind(pData.test.5[,c(1:3)],TrueLabel=pData.test.5$Class,Prediction_rfe.5)
write.table(Prediction_rfe.5, file = paste("TestSet05_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet05 with all predictors
matrix.test.5.full <- matrix.test.5[row.names(t(matrix.train.5)),]
Prediction_SVM_full.5 <- predict(svmFull.5,t(matrix.test.5.full), type = "prob")
Prediction_SVM_full.5$Prediction_SVM_full <- ifelse(Prediction_SVM_full.5$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.5 <- cbind(pData.test.5[,c(1:3)],TrueLabel=pData.test.5$Class,Prediction_SVM_full.5)
write.table(Prediction_SVM_full.5, file = paste("TestSet05_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 6 Performance of opVars Classifier on test samples TestSet 05 ##########################################
#############################################################################################################

#### 6.1 Confusion matrix  ##################################################################################  
sink("TestSet05_ConfusionMatrix_optVars_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.5$Prediction_rfe.5), as.factor(Prediction_rfe.5$TrueLabel))
sink()

sink("TestSet05_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.5$Prediction_SVM_full), as.factor(Prediction_SVM_full.5$TrueLabel))
sink()


#### 6.2 ROC on probability "transforming" TestSet 05 ########################################################

Prediction_rfe.5$Class <- as.factor(ifelse(Prediction_rfe.5$TrueLabel == "transforming","transforming","nontransforming"))

roc1 <- roc(Prediction_rfe.5$Class,                    # response vector (factor or character)
            Prediction_rfe.5$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=0.5)

Prediction_SVM_full.5$Class <- as.factor(ifelse(Prediction_SVM_full.5$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.5$Class,                    # response vector (factor or character)
               Prediction_SVM_full.5$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### Split 06 /  TestSet 06 #################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(1152)
split.6 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.6 <- eset.batch[, split.6] 
matrix.test.6  <- eset.batch[,-split.6] 
pData.train.6  <- pData[split.6,]
pData.test.6   <- pData[-split.6,]

table(pData.train.6$Design)
table(pData.test.6$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering TrainingSet 06 ######################################################
##############################################################################################################
fselect.6  <- genefilter(matrix.train.6, filterfun(f1))
summary(fselect.6)
matrix.train.6 <-matrix.train.6[fselect.6,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 06 #######################################################################
##############################################################################################################

matrix.train.6 <- (t(matrix.train.6))
labels.train.6 <- as.factor(pData.train.6$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(12345)
index.6 <- createMultiFolds(labels.train.6, k=10, times = 20)  

fullCtrl.6 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.6,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.6 <- train(matrix.train.6,labels.train.6,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.6)

svmFull.6  

##############################################################################################################
#### 3. SVM-RFE on TrainingSet 06 ############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.6      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               index = index.6,
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.6$functions         <- caretFuncs
outerctrl.6$functions$summary <- fiveStats

#### 3.2.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM
system.time(rfe.6  <- rfe(matrix.train.6, labels.train.6, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.6,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.6  # 19 optVars found 
write.table(rfe.6$results, file = "TrainingSet06_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.6, type = c("g", "o"))
plot(rfe.6, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.6 <- cbind(rfe.6$optVariables, Annotation[rfe.6$optVariables,])
write.table(optFeatures.rfe.6, file = "TrainingSet06_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 06 ######################################################
##############################################################################################################

rfeResamples.6 <- resamples(list("SVM_full.6" = svmFull.6,"SVM_RFE.6" = rfe.6))
sink("TrainingSet06_resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.6)
sink()

modelDifferences.6 <- diff(rfeResamples.6)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet06_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.6)
sink()


#############################################################################################################
##### 5  PCA on optimal variables ###########################################################################
#############################################################################################################
matrix.train.rfe.6     <- t(matrix.train.6[,rfe.6$optVariables])
matrix.test.rfe.6      <- matrix.test.6[rfe.6$optVariables,]
matrix.opt.6           <- cbind(matrix.train.rfe.6, matrix.test.rfe.6)
pData.opt.6            <- rbind(pData.train.6,pData.test.6)
pData.opt.6$Set        <- ifelse(row.names(pData.opt.6)%in% row.names(pData.train.6),"train","test")
pData.opt.6$Design_Color <- ifelse(pData.opt.6$Set == "train", pData.opt.6$Design_Color,"#000000")

pca.6         <- prcomp(t(matrix.opt.6))           
plot(pca.6$x, pch=16, col=pData.opt.6$Design_Color, cex=1.5, asp=1)
legend(-5,3.1, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.6$Design_Color), pch=16, bty="n", cex=0.8)

##############################################################################################################
##### 6 Train SVM on optVars determined by rfe and GA and predict independent TestSet 06 #####################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.train.rfe.6 <- matrix.train.6[,rfe.6$optVariables]   # subset fot the 19 optVars
matrix.test.rfe.6  <- t(matrix.test.6[rfe.6$optVariables,])

# train SVM and predict TestSet06 on 33 optVars found by rfe
set.seed(721)
svmOpt.rfe.6  <- train(matrix.train.rfe.6,labels.train.6,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.6)
svmOpt.rfe.6  

Prediction_rfe.6 <- predict(svmOpt.rfe.6,matrix.test.rfe.6, type = "prob")
Prediction_rfe.6$Prediction_rfe.6 <- ifelse(Prediction_rfe.6$transforming>0.50,"transforming","untransforming")
Prediction_rfe.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_rfe.6)
write.table(Prediction_rfe.6, file = paste("TestSet06_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet06 with all predictors
matrix.test.6.full <- matrix.test.6[row.names(t(matrix.train.6)),]
Prediction_SVM_full.6 <- predict(svmFull.6,t(matrix.test.6.full), type = "prob")
Prediction_SVM_full.6$Prediction_SVM_full <- ifelse(Prediction_SVM_full.6$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_SVM_full.6)
write.table(Prediction_SVM_full.6, file = paste("TestSet06_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 06 #####################################################
#############################################################################################################

#### 7.1 Confusion matrices  ################################################################################  
sink("TestSet06_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.6$Prediction_rfe.6), as.factor(Prediction_rfe.6$TrueLabel))
sink()

sink("TestSet06_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.6$Prediction_SVM_full), as.factor(Prediction_SVM_full.6$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 06 ########################################################

Prediction_rfe.6$Class <- as.factor(ifelse(Prediction_rfe.6$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.6$Class,                    
               Prediction_rfe.6$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_SVM_full.6$Class <- as.factor(ifelse(Prediction_SVM_full.6$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.6$Class,                    # response vector (factor or character)
               Prediction_SVM_full.6$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### Split 07 / TestSet 07 ##################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(543)
split.7 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.7 <- eset.batch[, split.7] 
matrix.test.7  <- eset.batch[,-split.7] 
pData.train.7  <- pData[split.7,]
pData.test.7   <- pData[-split.7,]

table(pData.train.7$Design)
table(pData.test.7$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering TrainingSet 07 ######################################################
##############################################################################################################
fselect.7  <- genefilter(matrix.train.7, filterfun(f1))
summary(fselect.7)
matrix.train.7 <-matrix.train.7[fselect.7,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 07 #######################################################################
##############################################################################################################

matrix.train.7 <- (t(matrix.train.7))
labels.train.7 <- as.factor(pData.train.7$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(847)
index.7 <- createMultiFolds(labels.train.7, k=10, times = 20)  

fullCtrl.7 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.7,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.7 <- train(matrix.train.7,labels.train.7,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.7)

svmFull.7  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 07 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.7      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               index = index.7,
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.7$functions         <- caretFuncs
outerctrl.7$functions$summary <- fiveStats


#### 3.2 SVM-RFE #############################################################################################
##############################################################################################################

system.time(rfe.7  <- rfe(matrix.train.7, labels.train.7, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.7,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.7    # 13 optVars found by rfe
write.table(rfe.7$results, file = "TrainingSet07_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.7, type = c("g", "o"))
plot(rfe.7, type = c("g", "o"), xlim = c(0,61))

optFeatures.7 <- cbind(rfe.7$optVariables, Annotation[rfe.7$optVariables,])
write.table(optFeatures.7, file = "TrainingSet07_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 07 ######################################################
##############################################################################################################

rfeResamples.7 <- resamples(list("SVM_full.7" = svmFull.7,"SVM_RFE.7" = rfe.7))
sink("TrainingSet07_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.7)
sink()

modelDifferences.7 <- diff(rfeResamples.7)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet07_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.7)
sink()

##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.7     <- t(matrix.train.7[,rfe.7$optVariables])
matrix.test.rfe.7      <- matrix.test.7[rfe.7$optVariables,]
matrix.opt.7           <- cbind(matrix.train.rfe.7, matrix.test.rfe.7)
pData.opt.7            <- rbind(pData.train.7,pData.test.7)
pData.opt.7$Set        <- ifelse(row.names(pData.opt.7)%in% row.names(pData.train.7),"train","test")
pData.opt.7$Design_Color <- ifelse(pData.opt.7$Set == "train", pData.opt.7$Design_Color,"#000000")

pca.7         <- prcomp(t(matrix.opt.7))           
plot(pca.7$x, pch=16, col=pData.opt.7$Design_Color, cex=1.5, asp=1)
legend(3,3.7, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.7$Design_Color), pch=16, bty="n", cex=0.8)


##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and predict independent TestSet 07 #############################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.train.rfe.7 <- matrix.train.7[,rfe.7$optVariables]   # subset fot the 33 optVars
matrix.test.rfe.7  <- t(matrix.test.7[rfe.7$optVariables,])

# train SVM and predict TestSet02 on 26 optVars found by rfe
set.seed(721)
svmOpt.rfe.7  <- train(matrix.train.rfe.7,labels.train.7,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.7)
svmOpt.rfe.7  

Prediction_rfe.7 <- predict(svmOpt.rfe.7,matrix.test.rfe.7, type = "prob")
Prediction_rfe.7$Prediction_rfe.7 <- ifelse(Prediction_rfe.7$transforming>0.50,"transforming","untransforming")
Prediction_rfe.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_rfe.7)
write.table(Prediction_rfe.7, file = paste("TestSet07_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet07 with all predictors
matrix.test.7.full <- matrix.test.7[row.names(t(matrix.train.7)),]
Prediction_SVM_full.7 <- predict(svmFull.7,t(matrix.test.7.full), type = "prob")
Prediction_SVM_full.7$Prediction_SVM_full <- ifelse(Prediction_SVM_full.7$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_SVM_full.7)
write.table(Prediction_SVM_full.7, file = paste("TestSet07_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 7 Performance of optVars Classifier on  TestSet 07 #####################################################
#############################################################################################################

#### 7.1 Confusion matrices  ################################################################################  
sink("TestSet07_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.7$Prediction_rfe.7), as.factor(Prediction_rfe.7$TrueLabel))
sink()

sink("TestSet07_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.7$Prediction_SVM_full), as.factor(Prediction_SVM_full.7$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 07 ########################################################

Prediction_rfe.7$Class <- as.factor(ifelse(Prediction_rfe.7$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.7$Class,                    
               Prediction_rfe.7$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_SVM_full.7$Class <- as.factor(ifelse(Prediction_SVM_full.7$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.7$Class,                    # response vector (factor or character)
               Prediction_SVM_full.7$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### Split 08/ TestSet 08 ###################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(1528)
split.8 <- createDataPartition(as.factor(pData$Class), p = .8, list = FALSE) # group-stratified sampling: 108 Training, 45 TestSamples

matrix.train.8 <- eset.batch[, split.8] 
matrix.test.8  <- eset.batch[,-split.8] 
pData.train.8  <- pData[split.8,]
pData.test.8   <- pData[-split.8,]

table(pData.train.8$Design)
table(pData.test.8$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering TrainingSet 08 ######################################################
##############################################################################################################
fselect.8  <- genefilter(matrix.train.8, filterfun(f1))
summary(fselect.8)
matrix.train.8 <-matrix.train.8[fselect.8,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 08 #######################################################################
##############################################################################################################
matrix.train.8 <- (t(matrix.train.8))
labels.train.8 <- as.factor(pData.train.8$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.8 <- createMultiFolds(labels.train.8, k=10, times = 20)  

fullCtrl.8 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.8,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.8 <- train(matrix.train.8,labels.train.8,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.8)

svmFull.8  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 08 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.8      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               index = index.8,
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.8$functions         <- caretFuncs
outerctrl.8$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 08 ##############################################################################
##############################################################################################################

system.time(rfe.8  <- rfe(matrix.train.8, labels.train.8, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.8,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.8   # 16 variables found by SVM-rfe
write.table(rfe.8$results, file = "TrainingSet08_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.8, type = c("g", "o"))
plot(rfe.8, type = c("g", "o"), xlim = c(0,61))

optFeatures.8 <- cbind(rfe.8$optVariables, Annotation[rfe.8$optVariables,])
write.table(optFeatures.8, file = "TrainingSet08_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 08 ######################################################
##############################################################################################################

rfeResamples.8 <- resamples(list("SVM_full.8" = svmFull.8,"SVM_RFE.8" = rfe.8))
sink("TrainingSet08_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.8)
sink()

modelDifferences.8 <- diff(rfeResamples.8)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet08_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.8)
sink()


##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.8     <- t(matrix.train.8[,rfe.8$optVariables])
matrix.test.rfe.8      <- matrix.test.8[rfe.8$optVariables,]
matrix.opt.8           <- cbind(matrix.train.rfe.8, matrix.test.rfe.8)
pData.opt.8            <- rbind(pData.train.8,pData.test.8)
pData.opt.8$Set        <- ifelse(row.names(pData.opt.8)%in% row.names(pData.train.8),"train","test")
pData.opt.8$Design_Color <- ifelse(pData.opt.8$Set == "train", pData.opt.8$Design_Color,"#000000")

pca.8         <- prcomp(t(matrix.opt.8))           
plot(pca.8$x, pch=16, col=pData.opt.8$Design_Color, cex=1.5, asp=1)
legend(2,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.8$Design_Color), pch=16, bty="n", cex=0.8)


##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet 08 ######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.train.rfe.8 <- matrix.train.8[,rfe.8$optVariables]   # subset fot the 20 optVars from rfe
matrix.test.rfe.8  <- t(matrix.test.8[rfe.8$optVariables,])

# train SVM and predict TestSet08 on 20 optVars found by rfe
set.seed(721)
svmOpt.rfe.8  <- train(matrix.train.rfe.8,labels.train.8,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.8)
svmOpt.rfe.8  

Prediction_rfe.8 <- predict(svmOpt.rfe.8,matrix.test.rfe.8, type = "prob")
Prediction_rfe.8$Prediction_rfe.8 <- ifelse(Prediction_rfe.8$transforming>0.50,"transforming","untransforming")
Prediction_rfe.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_rfe.8)
write.table(Prediction_rfe.8, file = paste("TestSet08_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


##### Predict TestSet08 with all predictors
matrix.test.8.full <- matrix.test.8[row.names(t(matrix.train.8)),]
Prediction_SVM_full.8 <- predict(svmFull.8,t(matrix.test.8.full), type = "prob")
Prediction_SVM_full.8$Prediction_SVM_full <- ifelse(Prediction_SVM_full.8$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_SVM_full.8)
write.table(Prediction_SVM_full.8, file = paste("TestSet08_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 7 Performance of optVars Classifier on  TestSet 08 #####################################################
#############################################################################################################

#### 7.1 Confusion matrices  ################################################################################  
sink("TestSet08_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.8$Prediction_rfe.8), as.factor(Prediction_rfe.8$TrueLabel))
sink()

sink("TestSet08_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.8$Prediction_SVM_full), as.factor(Prediction_SVM_full.8$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 08 ########################################################

Prediction_rfe.8$Class <- as.factor(ifelse(Prediction_rfe.8$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.8$Class,                    
               Prediction_rfe.8$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)

Prediction_SVM_full.8$Class <- as.factor(ifelse(Prediction_SVM_full.8$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.8$Class,                    # response vector (factor or character)
               Prediction_SVM_full.8$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### Split 09 / TestSet 09 ##################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(1218)
split.9 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.9 <- eset.batch[, split.9] 
matrix.test.9  <- eset.batch[,-split.9] 
pData.train.9  <- pData[split.9,]
pData.test.9   <- pData[-split.9,]

table(pData.train.9$Design)
table(pData.test.9$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering TrainingSet 09  #####################################################
##############################################################################################################
fselect.9  <- genefilter(matrix.train.9, filterfun(f1))
summary(fselect.9)
matrix.train.9 <-matrix.train.9[fselect.9,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 09 #######################################################################
##############################################################################################################
matrix.train.9 <- (t(matrix.train.9))
labels.train.9 <- as.factor(pData.train.9$Class)

set.seed(4567)
index.9 <- createMultiFolds(labels.train.9, k=10, times = 20)  

fullCtrl.9 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.9,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(721)
svmFull.9 <- train(matrix.train.9,labels.train.9,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.9)
svmFull.9  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 09 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.9      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               index = index.9,
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.9$functions         <- caretFuncs
outerctrl.9$functions$summary <- fiveStats

#### 3.2 SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.9  <- rfe(matrix.train.9, labels.train.9, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.9,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.9 # 6 optVars found by rfe
write.table(rfe.9$results, file = "TrainingSet09_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.9, type = c("g", "o"))
plot(rfe.9, type = c("g", "o"), xlim = c(0,61))

optFeatures.9 <- cbind(rfe.9$optVariables, Annotation[rfe.9$optVariables,])
write.table(optFeatures.9, file = "TrainingSet09_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 09 ######################################################
##############################################################################################################
rfeResamples.9 <- resamples(list("SVM_full.9" = svmFull.9,"SVM_RFE.9" = rfe.9))
sink("TrainingSet09_Resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.9)
sink()

modelDifferences.9 <- diff(rfeResamples.9)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet09_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.9)
sink()

##############################################################################################################
##### 4  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.9     <- t(matrix.train.9[,row.names(optFeatures.9)])
matrix.test.rfe.9      <- matrix.test.9[row.names(optFeatures.9),]
matrix.opt.9           <- cbind(matrix.train.rfe.9, matrix.test.rfe.9)
pData.opt.9            <- rbind(pData.train.9,pData.test.9)
pData.opt.9$Set        <- ifelse(row.names(pData.opt.9)%in% row.names(pData.train.9),"train","test")
pData.opt.9$Design_Color <- ifelse(pData.opt.9$Set == "train", pData.opt.9$Design_Color,"#000000")

pca.9         <- prcomp(t(matrix.opt.9))           
plot(pca.9$x, pch=16, col=pData.opt.9$Design_Color, cex=1.5, asp=1)
legend(3,3, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.9$Design_Color), pch=16, bty="n", cex=0.8)

##############################################################################################################
##### 5 Train SVM on optVars from rfe and predict independent TestSet 09 #####################################
##############################################################################################################

set.seed(721)
svmOpt.9  <- train(t(matrix.train.rfe.9),labels.train.9,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.9)
svmOpt.9  

Prediction_rfe.9 <- predict(svmOpt.9,t(matrix.test.rfe.9), type = "prob")
Prediction_rfe.9$Prediction_rfe.9 <- ifelse(Prediction_rfe.9$transforming>0.50,"transforming","untransforming")
Prediction_rfe.9 <- cbind(pData.test.9[,c(1:3)],TrueLabel=pData.test.9$Class,Prediction_rfe.9)
write.table(Prediction_rfe.9, file = paste("TestSet09_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet09 with all predictors
matrix.test.9.full <- matrix.test.9[row.names(t(matrix.train.9)),]
Prediction_SVM_full.9 <- predict(svmFull.9,t(matrix.test.9.full), type = "prob")
Prediction_SVM_full.9$Prediction_SVM_full <- ifelse(Prediction_SVM_full.9$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.9 <- cbind(pData.test.9[,c(1:3)],TrueLabel=pData.test.9$Class,Prediction_SVM_full.9)
write.table(Prediction_SVM_full.9, file = paste("TestSet09_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 6 Performance of opVars Classifier on test samples TestSet 09 ##########################################
#############################################################################################################

#### 6.1 Confusion matrix  ##################################################################################  
sink("TestSet09_ConfusionMatrix_optVars_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.9$Prediction_rfe.9), as.factor(Prediction_rfe.9$TrueLabel))
sink()

sink("TestSet09_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.9$Prediction_SVM_full), as.factor(Prediction_SVM_full.9$TrueLabel))
sink()

#### 6.2 ROC on probability "transforming" TestSet 09 ########################################################
Prediction_rfe.9$Class <- as.factor(ifelse(Prediction_rfe.9$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.9$Class,                    
               Prediction_rfe.9$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_SVM_full.9$Class <- as.factor(ifelse(Prediction_SVM_full.9$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.9$Class,                    # response vector (factor or character)
               Prediction_SVM_full.9$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


#############################################################################################################################################
#############################################################################################################################################
#### Split10 / TestSet 10 ###################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(1914)
split.10 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.10 <- eset.batch[, split.10] 
matrix.test.10  <- eset.batch[,-split.10] 
pData.train.10  <- pData[split.10,]
pData.test.10   <- pData[-split.10,]

table(pData.train.10$Design)
table(pData.test.10$Design)

##############################################################################################################
#### 1. nonspecific feature prefiltering #####################################################################
##############################################################################################################
fselect.10  <- genefilter(matrix.train.10, filterfun(f1))
summary(fselect.10)
matrix.train.10 <-matrix.train.10[fselect.10,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 10 #######################################################################
##############################################################################################################
matrix.train.10 <- (t(matrix.train.10))
labels.train.10 <- as.factor(pData.train.10$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1920)
index.10 <- createMultiFolds(labels.train.10, k=10, times = 20)  

fullCtrl.10 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.10,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.10 <- train(matrix.train.10,labels.train.10,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = fullCtrl.10)

svmFull.10  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 10 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################
outerctrl.10      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.10,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.10$functions         <- caretFuncs
outerctrl.10$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 10 ##################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.10  <- rfe(matrix.train.10, labels.train.10, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.10,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.10   # 17 optVars selected
write.table(rfe.10$results, file = "TrainingSet10_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.10, type = c("g", "o"))
plot(rfe.10, type = c("g", "o"), xlim = c(0,61))

optFeatures.10 <- cbind(rfe.10$optVariables, Annotation[rfe.10$optVariables,])
write.table(optFeatures.10, file = "TrainingSet10_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 10  #####################################################
##############################################################################################################
rfeResamples.10 <- resamples(list("SVM_full.10" = svmFull.10,"SVM_RFE.10" = rfe.10))
sink("TrainingSet10_Resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.10)
sink()

modelDifferences.10 <- diff(rfeResamples.10)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet10_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.10)
sink()


##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.10     <- t(matrix.train.10[,rfe.10$optVariables])
matrix.test.rfe.10      <- matrix.test.10[rfe.10$optVariables,]
matrix.opt.10           <- cbind(matrix.train.rfe.10, matrix.test.rfe.10)
pData.opt.10            <- rbind(pData.train.10,pData.test.10)
pData.opt.10$Set        <- ifelse(row.names(pData.opt.10)%in% row.names(pData.train.10),"train","test")
pData.opt.10$Design_Color <- ifelse(pData.opt.10$Set == "train", pData.opt.10$Design_Color,"#000000")

pca.10         <- prcomp(t(matrix.opt.10))           
plot(pca.10$x, pch=16, col=pData.opt.10$Design_Color, cex=1.5, asp=1)
legend(3.8,3, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.10$Design_Color), pch=16, bty="n", cex=0.8)


##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet 10 ######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.train.rfe.10 <- matrix.train.10[,rfe.10$optVariables]
matrix.test.rfe.10  <- t(matrix.test.10[rfe.10$optVariables,])

# train SVM and predict TestSet10 on 17 optVars found by rfe
set.seed(721)
svmOpt.rfe.10  <- train(matrix.train.rfe.10,labels.train.10,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.10)
svmOpt.rfe.10  

Prediction_rfe.10 <- predict(svmOpt.rfe.10,matrix.test.rfe.10, type = "prob")
Prediction_rfe.10$Prediction_rfe.10 <- ifelse(Prediction_rfe.10$transforming>0.50,"transforming","untransforming")
Prediction_rfe.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_rfe.10)
write.table(Prediction_rfe.10, file = paste("TestSet10_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet10 with all predictors
matrix.test.10.full <- matrix.test.10[row.names(t(matrix.train.10)),]
Prediction_SVM_full.10 <- predict(svmFull.10,t(matrix.test.10.full), type = "prob")
Prediction_SVM_full.10$Prediction_SVM_full <- ifelse(Prediction_SVM_full.10$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_SVM_full.10)
write.table(Prediction_SVM_full.10, file = paste("TestSet10_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 10 #####################################################
#############################################################################################################

#### 7.1 Confusion matrices  ################################################################################  
sink("TestSet10_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.10$Prediction_rfe.10), as.factor(Prediction_rfe.10$TrueLabel))
sink()

sink("TestSet10_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.10$Prediction_SVM_full), as.factor(Prediction_SVM_full.10$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 10 ########################################################
Prediction_rfe.10$Class <- as.factor(ifelse(Prediction_rfe.10$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe <- roc(Prediction_rfe.10$Class,                    
               Prediction_rfe.10$transforming,             
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


Prediction_SVM_full.10$Class <- as.factor(ifelse(Prediction_SVM_full.10$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.10$Class,                    # response vector (factor or character)
               Prediction_SVM_full.10$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)

#############################################################################################################################################
#############################################################################################################################################
#### FINAL MODEL WITHOUT INDEPENDENT TEST SET: Performance estimated only from resampling ###################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. nonspecific feature prefiltering FINAL  ##############################################################
##############################################################################################################

fselect.FINAL  <- genefilter(eset.batch, filterfun(f1))
summary(fselect.FINAL)
matrix.train.FINAL <- eset.batch[fselect.FINAL,]

##############################################################################################################
#### 2. SVM: FULL MODEL FINAL ################################################################################
##############################################################################################################
matrix.train.FINAL <- (t(matrix.train.FINAL))
labels.train.FINAL <- as.factor(pData$Class)

## create 200 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for comparison
set.seed(123)
index.FINAL <- createMultiFolds(labels.train.FINAL, k=10, times = 20)  

fullCtrl.FINAL <- trainControl(method = "repeatedcv",
                               repeats = 20,
                               summaryFunction = fiveStats,
                               classProbs = TRUE,
                               index = index.FINAL,
                               allowParallel = TRUE)

set.seed(721)
svmFull.FINAL <- train(matrix.train.FINAL,labels.train.FINAL,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.FINAL)

svmFull.FINAL  

## export Variable importances from the final set
varImp_FINAL <- varImp(svmFull.FINAL)$importance
varImp_FINAL <- varImp_FINAL[order(varImp_FINAL$transforming,decreasing = T),]
varImp_FINAL <- cbind(varImp_FINAL,Annotation[row.names(varImp_FINAL),])
write.table(varImp_FINAL, file = "FINAL_VariableImportance_ROC.txt", sep="\t",col.names=NA)

##############################################################################################################
#### 3. SVM-RFE FINAL ########################################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################
outerctrl.FINAL      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                   index = index.FINAL,
                                   saveDetails = TRUE,
                                   returnResamp="final", 
                                   verbose = TRUE, 
                                   seeds = seeds.rfe,
                                   allowParallel = TRUE)

outerctrl.FINAL$functions         <- caretFuncs
outerctrl.FINAL$functions$summary <- fiveStats

#### 3.2 SVM-RFE FINAL #######################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.FINAL  <- rfe(matrix.train.FINAL, labels.train.FINAL, 
                              sizes=FeatureNumbers,
                              rfeControl=outerctrl.FINAL,
                              metric = "Accuracy",
                              ## Options to train()
                              method="svmRadial",
                              tuneLength = 20,
                              trControl = innerctrl))

rfe.FINAL   # 32 optVars found by rfe
write.table(rfe.FINAL$results, file = "FinalSet152_Results_IQR1.2.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.FINAL, type = c("g", "o"))
plot(rfe.FINAL, type = c("g", "o"), xlim = c(0,61))

# Plot ROC over FeatureNumber
plot(rfe.FINAL$results$Variables,rfe.FINAL$results$ROC, type = "o",panel.first = grid())

optFeatures.FINAL <- cbind(rfe.FINAL$optVariables, Annotation[rfe.FINAL$optVariables,])
write.table(optFeatures.FINAL, file = "FinalSet152_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances FINAL ###############################################################
##############################################################################################################

rfeResamples.FINAL <- resamples(list("SVM_full.FINAL" = svmFull.FINAL,"SVM_RFE.FINAL" = rfe.FINAL))
sink("FinalSet152_Resamples_rfe.txt", append = TRUE)
summary(rfeResamples.FINAL)
sink()

modelDifferences.FINAL <- diff(rfeResamples.FINAL)  # paired t-test for H0: difference = 0 between the different models. 
sink("FinalSet152_ModelDifferences_rfe.txt", append = TRUE)
summary(modelDifferences.FINAL)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION FINALSET ###############################
################################################################################################
matrix.train.rfe.FINAL <- matrix.train.FINAL[,rfe.FINAL$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(1820)
index.GA.FINAL <- createMultiFolds(labels.train.FINAL, k=10, times = 5)  

set.seed(1821)
seeds.GA.FINAL <- vector(mode = "integer", length = length(index.GA.FINAL)+1)    # B+1 elements where B is the number of resamples = 51
for(i in 1:length(index.GA.FINAL)+1) seeds.GA.FINAL[[i]] <- sample.int(10000, 1)

outerctrl.GA.FINAL <- gafsControl(functions = svmGA,
                                  method = "repeatedcv", repeats = 5,
                                  index = index.GA.FINAL,                                       
                                  seeds = seeds.GA.FINAL,                                      
                                  returnResamp="all", 
                                  verbose = TRUE,
                                  maximize = c(internal = TRUE,
                                               external = TRUE),
                                  allowParallel = TRUE)                                  

#### 4.2 run GA  ###############################################################################
################################################################################################

system.time(GA.FINAL<- gafs(matrix.train.rfe.FINAL, labels.train.FINAL, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.FINAL,
                         ## Now we pass options to `train` via "svmGA":               
                         metric = "Accuracy",
                         method = "svmRadial",
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.3 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.FINAL  # 15 Features selected by GA
optVars.GA.FINAL <- Annotation[GA.FINAL$optVariables,]  # export optimal variables
write.table(optVars.GA.FINAL, file = "FinalSet152_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.FINAL <- GA.FINAL$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.FINAL <- arrange(performance.external.FINAL, Iter)

performance.FINAL <- GA.FINAL$ga$internal               # export internal accuracy (within the resample)
performance.FINAL$AccuracyExternal <- aggregate(performance.external.FINAL$Accuracy, by=list(Iter=performance.external.FINAL$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.FINAL <- data.frame(Iter = c(performance.FINAL$Iter,performance.FINAL$Iter), Accuracy = c(performance.FINAL$AccuracyExternal,performance.FINAL$Accuracy), Group=c(rep("external",40),rep("internal",40)))

## plot internal and external accuracy over the iterations 
ggplot(performance.long.FINAL, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.5,se = T) +
  theme_bw()

##############################################################################################################
##### 5  PCA on optimal variables FINALSET152 ################################################################
##############################################################################################################
matrix.train.rfe.FINAL <- t(matrix.train.FINAL[,GA.FINAL$optVariables])
pca.FINAL              <- prcomp(t(matrix.train.rfe.FINAL))           
plot(pca.FINAL$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(5,4.6, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

#############################################################################################################################################
#############################################################################################################################################
#### Analyze aggregated results over the different splits ###################################################################################
#############################################################################################################################################
#############################################################################################################################################

### 1. proportion of samples that were predicted at least once:  ############################################################################
#############################################################################################################################################
a <- unique(c(row.names(pData.test.1),row.names(pData.test.2),row.names(pData.test.3),row.names(pData.test.4),row.names(pData.test.5),
              row.names(pData.test.6),row.names(pData.test.7),row.names(pData.test.8),row.names(pData.test.9),row.names(pData.test.10)))
length(a)  # 141/152 = 93 % of the samples have been predicted at least once


### 2. Supplementary Table X: Resampling and TestSet performances for IQR = 0.5 #############################################################
#############################################################################################################################################

### create data.frame 11x9 filled with NAs
Results <- data.frame(Testset = c(seq(1,10,by=1),"FinalSet152"), FeatureNumber_full = rep(NA,11), ResamplingAccuracy_full=rep(NA,11), 
                      FeatureNumber_rfe = rep(NA,11), ResamplingAccuracy_rfe = rep(NA,11), FeatureNumber_GA = rep(NA,11),
                      Accuracy_TestSet_full = rep(NA,11),Accuracy_TestSet_rfe = rep(NA,11),Accuracy_TestSet_GA = rep(NA,11))

### fill in the number of features selected by IQR filter for each training/test split
for (i in 1:10) {Results[i,2] <- dim(eval(parse(text=paste("matrix.train.",i, sep=""))))[2] }
Results[11,2] <- dim(matrix.train.FINAL)[2]

### fill in the resampling accuracy for the full model
for (i in 1:10) { Results[i,3] <- eval(parse(text=paste("svmFull.",i, sep="")))$results[eval(parse(text=paste("svmFull.",i, sep="")))$results$C == eval(parse(text=paste("svmFull.",i, sep="")))$bestTune$C,"Accuracy"] }
Results[11,3] <- svmFull.FINAL$results[svmFull.FINAL$results$C == svmFull.FINAL$bestTune$C,"Accuracy"]


### fill in the number of features found by SVM-rfe for each training/test split
for (i in 1:10) {Results[i,4] <- eval(parse(text=paste("rfe.",i, sep="")))$optsize }
Results[11,4] <- rfe.FINAL$optsize

### fill in the resampling accuracy for the svm-rfe model
for (i in 1:10) { Results[i,5] <- eval(parse(text=paste("rfe.",i, sep="")))$results[eval(parse(text=paste("rfe.",i, sep="")))$results$Variables == eval(parse(text=paste("rfe.",i, sep="")))$optsize,5] }
Results[11,5] <- rfe.FINAL$results[rfe.FINAL$results$Variables == rfe.FINAL$optsize,5]

### fill in the number of features found by SVM-GA 
for (i in c(3)) { Results[i,6] <- length(eval(parse(text=paste("GA.",i, sep="")))$optVariables)  }
Results[11,6] <- 15

### fill in the TestSet accuracy for the full model
for (i in 1:10) {Results[i,7] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$Prediction_SVM_full), as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$TrueLabel))$overall[1] }

### fill in the TestSet accuracy for the SVM-rfe model
for (i in 1:10) { Results[i,8] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))$TrueLabel))$overall[1] }

### fill in the TestSet accuracy for the SVM-GA model
for (i in c(3)) { Results[i,9] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))$TrueLabel))$overall[1] }

write.table(Results, file = "AllResults_rfe_GA40_10TestSets_FinalSet152_IQR1.2.txt", sep="\t",col.names=NA)


