#############################################################################################################################################
#############################################################################################################################################
#### 19 Test Batches with TrainingSet processed with bapred upfront #########################################################################
#############################################################################################################################################
#############################################################################################################################################
library(genefilter)
library(Rtsne)
library(RColorBrewer)
library(caret)
library(kernlab)
library(ggplot2)
library(dplyr)
library(pROC)
library(bapred)
library(PRROC)
library(doMC)
registerDoMC(cores = 10)    

#############################################################################################################################################
#############################################################################################################################################
#### A) create RAW dataset of 154 arrays / 19 batches #######################################################################################
#############################################################################################################################################
#############################################################################################################################################

# all 169 arrays were read in and combined into a EListRaw object without further modification
# 15 samples with unknown ground truth were removed from the dataset, resulting in 154 assays
# IVIM ID 180523 was splitted due to severe class imbalance, the two mock controls were duplicated and used for 180523A and 180523B 
# One separate IVIM assay is treated as one batch except for IVIM #160706 and IVIM #160525 which comprise a single batch (7) since they have been processed together
# quantile normalization, averaging, log2 and batch correction are performed on the training and test sets separately to prevent data leakage

################################################################################################
#### 1. Read in data into EListRaw objects #####################################################
################################################################################################
pData.full       <- read.delim("SAGA_Targets_FINAL_169.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation       <- read.delim("Annotation_SAGA_FINAL_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)       # all 39428 probes 
Annotation.known <- read.delim("Annotation_SAGA_FINAL_KNOWN_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1) # 36226 probes with Annotation available

#### 1.1 Read in SAGA Assays on Agilent-048306  ################################################
################################################################################################
files.048306 <- pData.full$Filename[1:64]
system.time(RAW.048306 <- read.maimages(files=files.048306, path=".", source="agilent.median", green.only=T,columns=list(G="gMedianSignal"), annotation=c("ProbeName")))
colnames(RAW.048306)   <- row.names(pData.full)[1:64]  

## subset for the 39,428 probes from Agilent ID026655  
index1             <- RAW.048306$genes$ProbeName %in% row.names(Annotation)
RAW.048306_Agilent <- RAW.048306[index1,]                                            # 157,712 probes / 4 = 39,428 quadruplicate probes
RAW.048306_Agilent <- RAW.048306_Agilent[order(RAW.048306_Agilent$genes$ProbeName),] # order prior cbind

#### 1.2 Read in SAGA Assays on Agilent-066423  ################################################
################################################################################################
files.066423 <- pData.full$Filename[65:67]
RAW.066423   <- read.maimages(files=files.066423, path=".", source="agilent.median", green.only=T,columns=list(G="gMedianSignal"), annotation=c("ProbeName"))
colnames(RAW.066423) <- row.names(pData.full)[65:67]  

## subset for the 39428 probes from Agilent ID026655  
index2             <- RAW.066423$genes$ProbeName %in% row.names(Annotation)
RAW.066423_Agilent <- RAW.066423[index2,]  
RAW.066423_Agilent <- RAW.066423_Agilent[order(RAW.066423_Agilent$genes$ProbeName),] # order prior cbind

#### 1.3 Read in SAGA Assays on Agilent-084107  ################################################
################################################################################################
files.084107 <- pData.full$Filename[68:91]
RAW.084107   <- read.maimages(files=files.084107, path=".", source="agilent.median", green.only=T,columns=list(G="gMedianSignal"), annotation=c("ProbeName"))
colnames(RAW.084107) <- row.names(pData.full)[68:91]  

## subset for the 39,428 probes from Agilent ID026655  
index3             <- RAW.084107$genes$ProbeName %in% row.names(Annotation)
RAW.084107_Agilent <- RAW.084107[index3,]   
RAW.084107_Agilent <- RAW.084107_Agilent[order(RAW.084107_Agilent$genes$ProbeName),] # order prior cbind

#### 1.4 Read in SAGA Assays on Agilent-084956  ################################################
################################################################################################
files.084956 <- pData.full$Filename[92:169]
RAW.084956   <- read.maimages(files=files.084956, path=".", source="agilent.median", green.only=T, columns=list(G="gMedianSignal"), annotation=c("ProbeName"))
colnames(RAW.084956) <- row.names(pData.full)[92:169]  

## subset for the 39,428 probes from Agilent ID026655  
index4             <- RAW.084956$genes$ProbeName %in% row.names(Annotation)
RAW.084956_Agilent <- RAW.084956[index4,]  
RAW.084956_Agilent <- RAW.084956_Agilent[order(RAW.084956_Agilent$genes$ProbeName),] # order prior cbind

################################################################################################
#### 2. combine datasets and remove vectors with unknown ground truth (to few IVIMs) ###########
################################################################################################
table(RAW.048306_Agilent$genes$ProbeName == RAW.066423_Agilent$genes$ProbeName) # check whether all probeNames match
table(RAW.048306_Agilent$genes$ProbeName == RAW.084107_Agilent$genes$ProbeName)
table(RAW.048306_Agilent$genes$ProbeName == RAW.084956_Agilent$genes$ProbeName)

#### 2.1 RAW dataset complete ################################################################## 
################################################################################################
SAGA_RAW.full <- cbind.EList(RAW.048306_Agilent,RAW.066423_Agilent,RAW.084107_Agilent,RAW.084956_Agilent)
boxplot(log2(SAGA_RAW.full$E), col=pData.full$IVIM_Color, names=row.names(pData.full),boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)  

#### 2.2 RAW dataset with vectors of unknown class label removed ############################### 
################################################################################################
pData    <- subset(pData.full, !is.na(pData.full$Class)) # remove 15 samples
SAGA_RAW <- SAGA_RAW.full[,row.names(pData)]   # EListRaw object with 154 samples in 19 batches /  157,712 probes (39,428 unique probes in quadruplicate)       

all(row.names(pData) == colnames(SAGA_RAW))   
write.table(pData, file = "SAGA_154_Targets_batchwise.txt", sep="\t",row.names = TRUE, col.names = NA)

#### 2.3 clean up to free memory  ############################################################# 
###############################################################################################
rm(pData.full)
rm(SAGA_RAW.full)
rm(RAW.048306,RAW.048306_Agilent,RAW.066423,RAW.066423_Agilent,RAW.084107,RAW.084107_Agilent,RAW.084956,RAW.084956_Agilent)
rm(files.048306,files.066423,files.084107,files.084956)
rm(index1,index2,index3,index4)



#############################################################################################################################################
#############################################################################################################################################
#### I. TestSet 01 = Batch 01 = IVIM #120411 ################################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 01 and TestSet 01, ###########################################################
##############################################################################################################
pData.test.1  <- pData[pData$Batch==1,]               # select batch 1 / IVIM #120411 
pData.train.1 <- pData[pData$Batch!=1,]               # select all remaining assays 
RAW_train.1   <- SAGA_RAW[,row.names(pData.train.1)]  # use all remaining assays as training set
RAW_test.1    <- SAGA_RAW[,row.names(pData.test.1)]   # set aside batch 1 / IVIM #120411 as test set

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.1$E), col=pData.train.1$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.1 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.1    <- qunormtrain(t(RAW_train.1$E))
matrix.train.1.qn <- log2(t(qunorm.train.1$xnorm))
matrix.train.1.qn <- avereps(matrix.train.1.qn, ID= RAW_train.1$genes$ProbeName)  
colnames(matrix.train.1.qn) <- row.names(pData.train.1)
boxplot(matrix.train.1.qn, col=pData.train.1$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.1.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.1         <- as.factor(pData.train.1$Batch-1)
combat.train.1        <- combatba(t(matrix.train.1.qn), batch = batch.train.1)
matrix.train.1.batch  <- t(combat.train.1$xadj)
colnames(matrix.train.1.batch) <- row.names(pData.train.1)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(3)
plot(Rtsne(t(matrix.train.1.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$Design_Color, pch=16, cex=1.3) 


##############################################################################################################
#### 3. nonspecific feature prefiltering  ####################################################################
##############################################################################################################
f1       <- function(x) (IQR(x) > 0.8)    
fselect.1  <- genefilter(matrix.train.1.batch, filterfun(f1))
summary(fselect.1)                                  # 1309 probes selected with interquartile range of log2-int >0.8
matrix.train.1 <-matrix.train.1.batch[fselect.1,]   # subset 

##############################################################################################################
#### 4. SVM: FULL MODEL (1309 predictors) ####################################################################
##############################################################################################################
matrix.train.1 <- (t(matrix.train.1))
labels.train.1 <- as.factor(pData.train.1$Class)

# calculate performance measures (accuracy, sensitivity, specificity, ROC) of external resamples
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))   

## create 200 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for direct model comparison
set.seed(123)
index.1 <- createMultiFolds(labels.train.1, k=10, times = 20)  

## define parameters of the train function  
fullCtrl.1 <- trainControl(method = "repeatedcv",repeats = 20, # 10foldCVn20  
                           index = index.1,                    # fixed index to compare with the rfe model
                           summaryFunction = fiveStats,        # define summary function
                           classProbs = TRUE,                  # calculate class probabilities
                           allowParallel = TRUE)               # enable parallel computation 

## build full model on complete training data with all predictors    
set.seed(721)
system.time(svmFull.1 <- train(matrix.train.1,labels.train.1,  # define training set  
                               method = "svmRadial",                       # support vector machine with radial kernel
                               metric = "Accuracy",                        # use accuracy to select the best model
                               tuneLength = 20,                            # number of cost values to test (Caret creates a range of values and uses a single value of sigma that is calculated internally with kernlab “sigest” function) 
                               trControl = fullCtrl.1))
svmFull.1  
svmFull.1$results  

##############################################################################################################
#### 5. SVM-RFE on TrainingSet 01 ############################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

##  set predictors subset sizes to test: 1,2,…,40,45,50  = 43 subsets in total
FeatureNumbers <- c(seq(1,40,by=1),45,50)                 

## set all seeds for reproducibility: 52 seeds for each of the 200 resamples + 1 for the complete set 
set.seed(123)
seeds.rfe <- vector(mode = "list", length = length(index.1)+1)                            
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

#### 5.2 Parameters for inner (nested) cross-validation loop for hyperparameter tuning within each resample and for each subset size #######
##############################################################################################################
innerctrl <- trainControl(method = "repeatedcv",repeats = 3,           # 10CVn3 = 30 resamples
                          verboseIter = FALSE,
                          classProbs = TRUE,
                          allowParallel = FALSE)                      

#### 5.3.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of external holdouts x 30 internal resamples for tuning x 20 values for cost parameter = 6,240,000 steps
# this will take around 4-5 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.1  <- rfe(matrix.train.1, labels.train.1, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.1,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.1   # yields 22 features
write.table(rfe.1$results, file = "TrainingSet01_Results_rfe.txt", sep="\t",col.names=NA)
optFeatures.rfe.1 <- cbind(rfe.1$optVariables, Annotation[rfe.1$optVariables,])
write.table(optFeatures.rfe.1, file = "TrainingSet01_OptVars_rfe.txt", sep="\t",col.names=NA)

# Plot Accuracy over FeatureNumber
trellis.par.set(caretTheme())
plot(rfe.1, type = c("g", "o"))
plot(rfe.1, type = c("g", "o"), xlim = c(0,51))

#### 5.4.compare resampling performances between full model and rfe ##########################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.1 <- resamples(list("SVM_full.1" = svmFull.1,"SVM_RFE.1" = rfe.1))
sink("TrainingSet01_Resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.1)
sink()

modelDifferences.1 <- diff(rfeResamples.1)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet01_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.1)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 01 #########################
################################################################################################
matrix.train.rfe.1 <- matrix.train.1[,rfe.1$optVariables]   # subset fot the 20 optVars from rfe

#### 4.0 Set global SVM-GA Parameters:##########################################################
################################################################################################
svmGA <- caretGA  # predefined helper functions for the genetic algorithm

# define custom function to create an initial population with individuals that consist of 40% of the features on average to ensure efficient reduction in feature number
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

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.1 <- createMultiFolds(labels.train.1, k=10, times = 5)  

outerctrl.GA.1 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.1,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  


system.time(GA.1<- gafs(matrix.train.rfe.1, labels.train.1, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.1,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.1   # yields 10 features
optVars.GA.1 <- Annotation[GA.1$optVariables,]  # export optimal variables
write.table(optVars.GA.1, file = "TrainingSet01_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.1 <- GA.1$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.1 <- arrange(performance.external.1, Iter)

performance.1 <- GA.1$ga$internal               # export internal accuracy (within the resample)
performance.1$AccuracyExternal <- aggregate(performance.external.1$Accuracy, by=list(Iter=performance.external.1$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.1 <- data.frame(Iter = c(performance.1$Iter,performance.1$Iter), Accuracy = c(performance.1$AccuracyExternal,performance.1$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.1 <- GA.1$averages[GA.1$optIter,2]
accuracy.external.1

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.1 <- subset(performance.external.1,performance.external.1$Iter == GA.1$optIter)
accuracy.external.opt.1 <- accuracy.external.opt.1$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.1, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

### plot external accuracy over the iterations ################################################################
performance.aggr.1 <- subset(performance.long.1,performance.long.1$Group =="external")

ggplot(performance.aggr.1, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())



##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 01 ########################
##############################################################################################################
matrix.train.1.bapred.GA <- t(matrix.train.1.batch[GA.1$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.1  <- train(matrix.train.1.bapred.GA,labels.train.1,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.1)
svmOpt.GA.1  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.1.qn  <- log2(t(qunormaddon(qunorm.train.1, t(RAW_test.1$E))))   # use qunormtrain object to normalize test set
matrix.test.1.qn  <- avereps(matrix.test.1.qn, ID= RAW_test.1$genes$ProbeName)
boxplot(cbind(matrix.train.1.qn,matrix.test.1.qn),col=c(pData.train.1$IVIM_Color,pData.test.1$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.1.bapred  <- t(combatbaaddon(combat.train.1, t(matrix.test.1.qn), batch = as.factor(pData.test.1$Batch)))

plot(Rtsne(t(cbind(matrix.train.1.bapred,matrix.test.1.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
           col=c(pData.train.1$Design_Color,pData.test.1$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.1.bapred.GA   <- t(matrix.test.1.bapred[GA.1$optVariables,])
matrix.test.1.bapred.full <- t(matrix.test.1.bapred[colnames(matrix.train.1),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.1.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.1$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.1$Design), col=unique(pData.train.1$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.1.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.1$Design_Color,pData.test.1$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.1), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.1 <- predict(svmOpt.GA.1,matrix.test.1.bapred.GA, type = "prob")
Prediction_GA.1$Prediction_GA.1 <- ifelse(Prediction_GA.1$transforming>0.50,"transforming","untransforming")
Prediction_GA.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_GA.1)
write.table(Prediction_GA.1, file = paste("TestSet01_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.1 <- predict(svmFull.1, matrix.test.1.bapred.full, type = "prob")
Prediction_SVM_full.1$Prediction_SVM_full <- ifelse(Prediction_SVM_full.1$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_SVM_full.1)
write.table(Prediction_SVM_full.1, file = paste("TestSet01_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 01 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.1  ###########################  
sink("TestSet01_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.1$Prediction_GA.1), as.factor(Prediction_GA.1$TrueLabel))
sink()

sink("TestSet01_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.1$Prediction_SVM_full), as.factor(Prediction_SVM_full.1$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 01 ########################################################

Prediction_GA.1$Class <- as.factor(ifelse(Prediction_GA.1$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.1 <- roc(Prediction_GA.1$Class,                    # response vector (factor or character)
              Prediction_GA.1$transforming,             # predictor vector (numeric)
              percent=TRUE, levels=c("nontransforming","transforming"),
              plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
              print.auc=T)

Prediction_SVM_full.1$Class <- as.factor(ifelse(Prediction_SVM_full.1$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.1<- roc(Prediction_SVM_full.1$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.1$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 01 ####################################
Prediction_GA.1$Class_Code <- ifelse(Prediction_GA.1$TrueLabel == "transforming",1,0)
pr.GA.1 <- pr.curve(scores.class0 = Prediction_GA.1$transforming , weights.class0 = Prediction_GA.1$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.1, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.1$Class_Code <- ifelse(Prediction_SVM_full.1$TrueLabel == "transforming",1,0)
pr.full.1 <- pr.curve(scores.class0 = Prediction_SVM_full.1$transforming , weights.class0 = Prediction_SVM_full.1$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.1, rand.plot = TRUE, legend = F, color = 1,main = "")



#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 02 = Batch 02 = IVIM #150128 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 02 and TestSet 02, ###########################################################
##############################################################################################################
pData.test.2  <- pData[pData$Batch==2,]    # set aside batch 2 / IVIM #220422
pData.train.2 <- pData[pData$Batch!=2,]    # use all remaining assays as training set
RAW_train.2   <- SAGA_RAW[,row.names(pData.train.2)]
RAW_test.2    <- SAGA_RAW[,row.names(pData.test.2)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.2$E), col=pData.train.2$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.2    <- qunormtrain(t(RAW_train.2$E))
matrix.train.2.qn <- log2(t(qunorm.train.2$xnorm))
matrix.train.2.qn <- avereps(matrix.train.2.qn, ID= RAW_train.2$genes$ProbeName)  
colnames(matrix.train.2.qn) <- row.names(pData.train.2)
boxplot(matrix.train.2.qn, col=pData.train.2$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.2.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.2$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.2         <- as.factor(ifelse(pData.train.2$Batch>2,pData.train.2$Batch-1,pData.train.2$Batch))                       
combat.train.2        <- combatba(t(matrix.train.2.qn), batch = batch.train.2)
matrix.train.2.batch  <- t(combat.train.2$xadj)
colnames(matrix.train.2.batch) <- row.names(pData.train.2)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.2.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.2$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 02 ######################################################
##############################################################################################################

fselect.2  <- genefilter(matrix.train.2.batch, filterfun(f1))
summary(fselect.2)
matrix.train.2 <-matrix.train.2.batch[fselect.2,]

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

system.time(rfe.2 <- rfe(matrix.train.2, labels.train.2, 
                         sizes=FeatureNumbers,
                         rfeControl=outerctrl.2,
                         metric = "Accuracy",
                         method="svmRadial",
                         tuneLength = 20,
                         trControl = innerctrl))

rfe.2
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
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 02 ########################
##############################################################################################################
matrix.train.2.bapred.rfe  <- t(matrix.train.2.bapred[rfe.2$optVariables,])# subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.rfe.2  <- train(matrix.train.2.bapred.rfe,labels.train.2,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.2)
svmOpt.rfe.2  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.2.qn  <- log2(t(qunormaddon(qunorm.train.2, t(RAW_test.2$E))))   # use qunormtrain object to normalize test set
matrix.test.2.qn  <- avereps(matrix.test.2.qn, ID= RAW_test.2$genes$ProbeName)
boxplot(cbind(matrix.train.2.qn,matrix.test.2.qn),col=c(pData.train.2$IVIM_Color,pData.test.2$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.2.bapred  <- t(combatbaaddon(combat.train.2, t(matrix.test.2.qn), batch = as.factor(pData.test.2$Batch)))

plot(Rtsne(t(cbind(matrix.train.2.bapred,matrix.test.2.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.2$Design_Color,pData.test.2$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.2.bapred.rfe   <- t(matrix.test.2.bapred[rfe.2$optVariables,])
matrix.test.2.bapred.full <- t(matrix.test.2.bapred[colnames(matrix.train.2),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.2.bapred.rfe, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.2$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.2$Design), col=unique(pData.train.2$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.2.bapred.rfe)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.2$Design_Color,pData.test.2$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.2), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_rfe.2 <- predict(svmOpt.rfe.2,matrix.test.2.bapred.rfe, type = "prob")
Prediction_rfe.2$Prediction_rfe.2 <- ifelse(Prediction_rfe.2$transforming>0.50,"transforming","untransforming")
Prediction_rfe.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_rfe.2)
write.table(Prediction_GA.2, file = paste("TestSet02_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.2 <- predict(svmFull.2, matrix.test.2.bapred.full, type = "prob")
Prediction_SVM_full.2$Prediction_SVM_full <- ifelse(Prediction_SVM_full.2$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_SVM_full.2)
write.table(Prediction_SVM_full.2, file = paste("TestSet02_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 02 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet02_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.2$Prediction_rfe.2), as.factor(Prediction_rfe.2$TrueLabel))
sink()

sink("TestSet02_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.2$Prediction_SVM_full), as.factor(Prediction_SVM_full.2$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 02 ########################################################

Prediction_rfe.2$Class <- as.factor(ifelse(Prediction_rfe.2$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe.2 <- roc(Prediction_rfe.2$Class,                    # response vector (factor or character)
                 Prediction_rfe.2$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.2$Class <- as.factor(ifelse(Prediction_SVM_full.2$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.2<- roc(Prediction_SVM_full.2$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.2$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 02 ####################################
Prediction_rfe.2$Class_Code <- ifelse(Prediction_rfe.2$TrueLabel == "transforming",1,0)
pr.rfe.2 <- pr.curve(scores.class0 = Prediction_rfe.2$transforming , weights.class0 = Prediction_rfe.2$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.rfe.2, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.2$Class_Code <- ifelse(Prediction_SVM_full.2$TrueLabel == "transforming",1,0)
pr.full.2 <- pr.curve(scores.class0 = Prediction_SVM_full.2$transforming , weights.class0 = Prediction_SVM_full.2$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.2, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### III. TestSet 03 = batch 03 / IVIM ID 150304 ############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 03 and TestSet 03 ############################################################
##############################################################################################################
pData.test.3  <- pData[pData$Batch==3,]    # set aside batch 3 / IVIM #150304
pData.train.3 <- pData[pData$Batch!=3,]    # use all remaining assays as training set
RAW_train.3   <- SAGA_RAW[,row.names(pData.train.3)]
RAW_test.3    <- SAGA_RAW[,row.names(pData.test.3)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.3$E), col=pData.train.3$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.3 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.3    <- qunormtrain(t(RAW_train.3$E))
matrix.train.3.qn <- log2(t(qunorm.train.3$xnorm))
matrix.train.3.qn <- avereps(matrix.train.3.qn, ID= RAW_train.3$genes$ProbeName)  
colnames(matrix.train.3.qn) <- row.names(pData.train.3)
boxplot(matrix.train.3.qn, col=pData.train.3$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.3.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.3$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.3         <- as.factor(ifelse(pData.train.3$Batch>3,pData.train.3$Batch-1,pData.train.3$Batch))                       
combat.train.3        <- combatba(t(matrix.train.3.qn), batch = batch.train.3)
matrix.train.3.batch  <- t(combat.train.3$xadj)
colnames(matrix.train.3.batch) <- row.names(pData.train.3)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.3.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.3$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 1. nonspecific feature prefiltering in TrainingSet 03 ###################################################
##############################################################################################################
fselect.3  <- genefilter(matrix.train.3.batch, filterfun(f1))
summary(fselect.3)
matrix.train.3 <-matrix.train.3.batch[fselect.3,]

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

system.time(rfe.3  <- rfe(matrix.train.3, labels.train.3, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.3,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.3    # 15 predictors found 
write.table(rfe.3$results, file = "TrainingSet03_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.3, type = c("g", "o"))
plot(rfe.3, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.3 <- cbind(rfe.3$optVariables, Annotation[rfe.3$optVariables,])
write.table(optFeatures.rfe.3, file = "TrainingSet03_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3.compare resampling performances #####################################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

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
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
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

# extract average resampling accuracy at optimal iteration 
accuracy.external.3 <- GA.3$averages[GA.3$optIter,2]

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.3 <- subset(performance.external.3,performance.external.3$Iter == GA.3$optIter)
accuracy.external.opt.3 <- accuracy.external.opt.3$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.3, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.5,se = T) +
  theme_bw()

### plot external accuracy over the iterations ################################################################
performance.aggr.3 <- subset(performance.long.3,performance.long.3$Group =="external")

ggplot(performance.aggr.3, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())



##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 03 ########################
##############################################################################################################
matrix.train.3.bapred.GA <- t(matrix.train.3.batch[GA.3$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.3  <- train(matrix.train.3.bapred.GA,labels.train.3,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.3)
svmOpt.GA.3  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.3.qn  <- log2(t(qunormaddon(qunorm.train.3, t(RAW_test.3$E))))   # use qunormtrain object to normalize test set
matrix.test.3.qn  <- avereps(matrix.test.3.qn, ID= RAW_test.3$genes$ProbeName)
boxplot(cbind(matrix.train.3.qn,matrix.test.3.qn),col=c(pData.train.3$IVIM_Color,pData.test.3$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.3.bapred  <- t(combatbaaddon(combat.train.3, t(matrix.test.3.qn), batch = as.factor(pData.test.3$Batch)))

plot(Rtsne(t(cbind(matrix.train.3.bapred,matrix.test.3.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.3$Design_Color,pData.test.3$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.3.bapred.GA   <- t(matrix.test.3.bapred[GA.3$optVariables,])
matrix.test.3.bapred.full <- t(matrix.test.3.bapred[colnames(matrix.train.3),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.3.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.3$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.3$Design), col=unique(pData.train.3$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.3.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.3$Design_Color,pData.test.3$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.3), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.3 <- predict(svmOpt.GA.3,matrix.test.3.bapred.GA, type = "prob")
Prediction_GA.3$Prediction_GA.3 <- ifelse(Prediction_GA.3$transforming>0.50,"transforming","untransforming")
Prediction_GA.3 <- cbind(pData.test.3[,c(1:3)],TrueLabel=pData.test.3$Class,Prediction_GA.3)
write.table(Prediction_GA.3, file = paste("TestSet03_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.3 <- predict(svmFull.3, matrix.test.3.bapred.full, type = "prob")
Prediction_SVM_full.3$Prediction_SVM_full <- ifelse(Prediction_SVM_full.3$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.3 <- cbind(pData.test.3[,c(1:3)],TrueLabel=pData.test.3$Class,Prediction_SVM_full.3)
write.table(Prediction_SVM_full.3, file = paste("TestSet03_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 03 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet03_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.3$Prediction_GA.3), as.factor(Prediction_GA.3$TrueLabel))
sink()

sink("TestSet03_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.3$Prediction_SVM_full), as.factor(Prediction_SVM_full.3$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 03 ########################################################

Prediction_GA.3$Class <- as.factor(ifelse(Prediction_GA.3$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.3 <- roc(Prediction_GA.3$Class,                    # response vector (factor or character)
                Prediction_GA.3$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.3$Class <- as.factor(ifelse(Prediction_SVM_full.3$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.3<- roc(Prediction_SVM_full.3$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.3$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 03 ####################################
Prediction_GA.3$Class_Code <- ifelse(Prediction_GA.3$TrueLabel == "transforming",1,0)
pr.GA.3 <- pr.curve(scores.class0 = Prediction_GA.3$transforming , weights.class0 = Prediction_GA.3$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.3, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.3$Class_Code <- ifelse(Prediction_SVM_full.3$TrueLabel == "transforming",1,0)
pr.full.3 <- pr.curve(scores.class0 = Prediction_SVM_full.3$transforming , weights.class0 = Prediction_SVM_full.3$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.3, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### IV. TestSet 04 = Batch 04 / IVIM #151014 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 04 and TestSet 04, ###########################################################
##############################################################################################################
pData.test.4  <- pData[pData$Batch==4,]    # set aside batch 4 / IVIM #151014
pData.train.4 <- pData[pData$Batch!=4,]    # use all remaining assays as training set
RAW_train.4   <- SAGA_RAW[,row.names(pData.train.4)]
RAW_test.4    <- SAGA_RAW[,row.names(pData.test.4)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.4$E), col=pData.train.4$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.4 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.4    <- qunormtrain(t(RAW_train.4$E))
matrix.train.4.qn <- log2(t(qunorm.train.4$xnorm))
matrix.train.4.qn <- avereps(matrix.train.4.qn, ID= RAW_train.4$genes$ProbeName)  
colnames(matrix.train.4.qn) <- row.names(pData.train.4)
boxplot(matrix.train.4.qn, col=pData.train.4$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.4.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.4$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.4         <- as.factor(ifelse(pData.train.4$Batch>4,pData.train.4$Batch-1,pData.train.4$Batch))                       
combat.train.4        <- combatba(t(matrix.train.4.qn), batch = batch.train.4)
matrix.train.4.batch  <- t(combat.train.4$xadj)
colnames(matrix.train.4.batch) <- row.names(pData.train.4)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.4.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.4$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 04 ######################################################
##############################################################################################################
fselect.4  <- genefilter(matrix.train.4.batch, filterfun(f1))
summary(fselect.4)
matrix.train.4 <-matrix.train.4.batch[fselect.4,]

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

system.time(rfe.4  <- rfe(matrix.train.4, labels.train.4, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.4,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.4    # 22 predictors chosen
write.table(rfe.4$results, file = "TrainingSet04_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.4, type = c("g", "o"))
plot(rfe.4, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.4 <- cbind(rfe.4$optVariables, Annotation[rfe.4$optVariables,])
write.table(optFeatures.rfe.4, file = "TrainingSet04_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances #####################################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.4 <- resamples(list("SVM_full.4" = svmFull.4,"SVM_RFE.4" = rfe.4))
sink("TrainingSet04_Resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.4)
sink()

modelDifferences.4 <- diff(rfeResamples.4)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet04_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.4)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 04 #########################
################################################################################################

matrix.train.rfe.4 <- matrix.train.4[,rfe.4$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.4 <- createMultiFolds(labels.train.4, k=10, times = 5)  

outerctrl.GA.4 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.4,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  


system.time(GA.4<- gafs(matrix.train.rfe.4, labels.train.4, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.4,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.4   # yields 10 features
optVars.GA.4 <- Annotation[GA.4$optVariables,]  # export optimal variables
write.table(optVars.GA.4, file = "TrainingSet04_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.4 <- GA.4$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.4 <- arrange(performance.external.4, Iter)

performance.4 <- GA.4$ga$internal               # export internal accuracy (within the resample)
performance.4$AccuracyExternal <- aggregate(performance.external.4$Accuracy, by=list(Iter=performance.external.4$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.4 <- data.frame(Iter = c(performance.4$Iter,performance.4$Iter), Accuracy = c(performance.4$AccuracyExternal,performance.4$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.4 <- GA.4$averages[GA.4$optIter,2]
accuracy.external.4

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.4 <- subset(performance.external.4,performance.external.4$Iter == GA.4$optIter)
accuracy.external.opt.4 <- accuracy.external.opt.4$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.4, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

### plot external accuracy over the iterations ################################################################
performance.aggr.4 <- subset(performance.long.4,performance.long.4$Group =="external")

ggplot(performance.aggr.4, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())



##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 04 ########################
##############################################################################################################
matrix.train.4.bapred.GA <- t(matrix.train.4.batch[GA.4$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.4  <- train(matrix.train.4.bapred.GA,labels.train.4,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.4)
svmOpt.GA.4  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.4.qn  <- log2(t(qunormaddon(qunorm.train.4, t(RAW_test.4$E))))   # use qunormtrain object to normalize test set
matrix.test.4.qn  <- avereps(matrix.test.4.qn, ID= RAW_test.4$genes$ProbeName)
boxplot(cbind(matrix.train.4.qn,matrix.test.4.qn),col=c(pData.train.4$IVIM_Color,pData.test.4$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.4.bapred  <- t(combatbaaddon(combat.train.4, t(matrix.test.4.qn), batch = as.factor(pData.test.4$Batch)))

plot(Rtsne(t(cbind(matrix.train.4.bapred,matrix.test.4.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.4$Design_Color,pData.test.4$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.4.bapred.GA   <- t(matrix.test.4.bapred[GA.4$optVariables,])
matrix.test.4.bapred.full <- t(matrix.test.4.bapred[colnames(matrix.train.4),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.4.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.4$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.4$Design), col=unique(pData.train.4$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.4.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.4$Design_Color,pData.test.4$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.4), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.4 <- predict(svmOpt.GA.4,matrix.test.4.bapred.GA, type = "prob")
Prediction_GA.4$Prediction_GA.4 <- ifelse(Prediction_GA.4$transforming>0.50,"transforming","untransforming")
Prediction_GA.4 <- cbind(pData.test.4[,c(1:3)],TrueLabel=pData.test.4$Class,Prediction_GA.4)
write.table(Prediction_GA.4, file = paste("TestSet04_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.4 <- predict(svmFull.4, matrix.test.4.bapred.full, type = "prob")
Prediction_SVM_full.4$Prediction_SVM_full <- ifelse(Prediction_SVM_full.4$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.4 <- cbind(pData.test.4[,c(1:3)],TrueLabel=pData.test.4$Class,Prediction_SVM_full.4)
write.table(Prediction_SVM_full.4, file = paste("TestSet04_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 04 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet04_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.4$Prediction_GA.4), as.factor(Prediction_GA.4$TrueLabel))
sink()

sink("TestSet04_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.4$Prediction_SVM_full), as.factor(Prediction_SVM_full.4$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 04 ########################################################

Prediction_GA.4$Class <- as.factor(ifelse(Prediction_GA.4$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.4 <- roc(Prediction_GA.4$Class,                    # response vector (factor or character)
                Prediction_GA.4$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.4$Class <- as.factor(ifelse(Prediction_SVM_full.4$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.4<- roc(Prediction_SVM_full.4$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.4$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 04 ####################################
Prediction_GA.4$Class_Code <- ifelse(Prediction_GA.4$TrueLabel == "transforming",1,0)
pr.GA.4 <- pr.curve(scores.class0 = Prediction_GA.4$transforming , weights.class0 = Prediction_GA.4$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.4, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.4$Class_Code <- ifelse(Prediction_SVM_full.4$TrueLabel == "transforming",1,0)
pr.full.4 <- pr.curve(scores.class0 = Prediction_SVM_full.4$transforming , weights.class0 = Prediction_SVM_full.4$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.4, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 05 = Batch 05 = IVIM #160210 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 05 and TestSet 05 ############################################################
##############################################################################################################
pData.test.5  <- pData[pData$Batch==5,]    # set aside batch 5 / IVIM #160210
pData.train.5 <- pData[pData$Batch!=5,]    # use all remaining assays as training set
RAW_train.5   <- SAGA_RAW[,row.names(pData.train.5)]
RAW_test.5    <- SAGA_RAW[,row.names(pData.test.5)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.5$E), col=pData.train.5$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.5 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.5    <- qunormtrain(t(RAW_train.5$E))
matrix.train.5.qn <- log2(t(qunorm.train.5$xnorm))
matrix.train.5.qn <- avereps(matrix.train.5.qn, ID= RAW_train.5$genes$ProbeName)  
colnames(matrix.train.5.qn) <- row.names(pData.train.5)
boxplot(matrix.train.5.qn, col=pData.train.5$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.5.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.5$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.5         <- as.factor(ifelse(pData.train.5$Batch>5,pData.train.5$Batch-1,pData.train.5$Batch))                       
combat.train.5        <- combatba(t(matrix.train.5.qn), batch = batch.train.5)
matrix.train.5.batch  <- t(combat.train.5$xadj)
colnames(matrix.train.5.batch) <- row.names(pData.train.5)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.5.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.5$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 02 ######################################################
##############################################################################################################
fselect.5  <- genefilter(matrix.train.5.batch, filterfun(f1))
summary(fselect.5)
matrix.train.5 <-matrix.train.5.batch[fselect.5,]
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

system.time(rfe.5  <- rfe(matrix.train.5, labels.train.5, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.5,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.5    # 50 optimal predictors found 
write.table(rfe.5$results, file = "TrainingSet05_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.5, type = c("g", "o"))
plot(rfe.5, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.5 <- cbind(rfe.5$optVariables, Annotation[rfe.5$optVariables,])
write.table(optFeatures.rfe.5, file = "TrainingSet05_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances on TrainingSet 05 ###################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.5 <- resamples(list("SVM_full.5" = svmFull.5,"SVM_RFE.5" = rfe.5))
sink("TrainingSet05_resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.5)
sink()

modelDifferences.5 <- diff(rfeResamples.5)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet05_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.5)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 05 #########################
################################################################################################

matrix.train.rfe.5 <- matrix.train.5[,rfe.5$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.5 <- createMultiFolds(labels.train.5, k=10, times = 5)  

outerctrl.GA.5 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.5,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  


system.time(GA.5<- gafs(matrix.train.rfe.5, labels.train.5, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.5,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.5   # yields 10 features
optVars.GA.5 <- Annotation[GA.5$optVariables,]  # export optimal variables
write.table(optVars.GA.5, file = "TrainingSet05_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.5 <- GA.5$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.5 <- arrange(performance.external.5, Iter)

performance.5 <- GA.5$ga$internal               # export internal accuracy (within the resample)
performance.5$AccuracyExternal <- aggregate(performance.external.5$Accuracy, by=list(Iter=performance.external.5$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.5 <- data.frame(Iter = c(performance.5$Iter,performance.5$Iter), Accuracy = c(performance.5$AccuracyExternal,performance.5$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.5 <- GA.5$averages[GA.5$optIter,2]
accuracy.external.5

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.5 <- subset(performance.external.5,performance.external.5$Iter == GA.5$optIter)
accuracy.external.opt.5 <- accuracy.external.opt.5$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.5, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()



##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 05 ########################
##############################################################################################################
matrix.train.5.bapred.GA <- t(matrix.train.5.batch[GA.5$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.5  <- train(matrix.train.5.bapred.GA,labels.train.5,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.5)
svmOpt.GA.5  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.5.qn  <- log2(t(qunormaddon(qunorm.train.5, t(RAW_test.5$E))))   # use qunormtrain object to normalize test set
matrix.test.5.qn  <- avereps(matrix.test.5.qn, ID= RAW_test.5$genes$ProbeName)
boxplot(cbind(matrix.train.5.qn,matrix.test.5.qn),col=c(pData.train.5$IVIM_Color,pData.test.5$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.5.bapred  <- t(combatbaaddon(combat.train.5, t(matrix.test.5.qn), batch = as.factor(pData.test.5$Batch)))

plot(Rtsne(t(cbind(matrix.train.5.bapred,matrix.test.5.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.5$Design_Color,pData.test.5$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.5.bapred.GA   <- t(matrix.test.5.bapred[GA.5$optVariables,])
matrix.test.5.bapred.full <- t(matrix.test.5.bapred[colnames(matrix.train.5),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.5.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.5$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.5$Design), col=unique(pData.train.5$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.5.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.5$Design_Color,pData.test.5$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.5), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.5 <- predict(svmOpt.GA.5,matrix.test.5.bapred.GA, type = "prob")
Prediction_GA.5$Prediction_GA.5 <- ifelse(Prediction_GA.5$transforming>0.50,"transforming","untransforming")
Prediction_GA.5 <- cbind(pData.test.5[,c(1:3)],TrueLabel=pData.test.5$Class,Prediction_GA.5)
write.table(Prediction_GA.5, file = paste("TestSet05_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.5 <- predict(svmFull.5, matrix.test.5.bapred.full, type = "prob")
Prediction_SVM_full.5$Prediction_SVM_full <- ifelse(Prediction_SVM_full.5$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.5 <- cbind(pData.test.5[,c(1:3)],TrueLabel=pData.test.5$Class,Prediction_SVM_full.5)
write.table(Prediction_SVM_full.5, file = paste("TestSet05_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 05 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet05_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.5$Prediction_GA.5), as.factor(Prediction_GA.5$TrueLabel))
sink()

sink("TestSet05_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.5$Prediction_SVM_full), as.factor(Prediction_SVM_full.5$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 05 ########################################################

Prediction_GA.5$Class <- as.factor(ifelse(Prediction_GA.5$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.5 <- roc(Prediction_GA.5$Class,                    # response vector (factor or character)
                Prediction_GA.5$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.5$Class <- as.factor(ifelse(Prediction_SVM_full.5$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.5<- roc(Prediction_SVM_full.5$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.5$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 05 ####################################
Prediction_GA.5$Class_Code <- ifelse(Prediction_GA.5$TrueLabel == "transforming",1,0)
pr.GA.5 <- pr.curve(scores.class0 = Prediction_GA.5$transforming , weights.class0 = Prediction_GA.5$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.5, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.5$Class_Code <- ifelse(Prediction_SVM_full.5$TrueLabel == "transforming",1,0)
pr.full.5 <- pr.curve(scores.class0 = Prediction_SVM_full.5$transforming , weights.class0 = Prediction_SVM_full.5$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.5, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 06 = Batch 06 = IVIM #160413 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 06 and TestSet 06 ############################################################
##############################################################################################################
pData.test.6  <- pData[pData$Batch==6,]    # set aside batch 6 / IVIM #160413
pData.train.6 <- pData[pData$Batch!=6,]    # use all remaining assays as training set
RAW_train.6   <- SAGA_RAW[,row.names(pData.train.6)]
RAW_test.6    <- SAGA_RAW[,row.names(pData.test.6)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.6$E), col=pData.train.6$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.6 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.6    <- qunormtrain(t(RAW_train.6$E))
matrix.train.6.qn <- log2(t(qunorm.train.6$xnorm))
matrix.train.6.qn <- avereps(matrix.train.6.qn, ID= RAW_train.6$genes$ProbeName)  
colnames(matrix.train.6.qn) <- row.names(pData.train.6)
boxplot(matrix.train.6.qn, col=pData.train.6$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.6.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.6$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.6         <- as.factor(ifelse(pData.train.6$Batch>6,pData.train.6$Batch-1,pData.train.6$Batch))                       
combat.train.6        <- combatba(t(matrix.train.6.qn), batch = batch.train.6)
matrix.train.6.batch  <- t(combat.train.6$xadj)
colnames(matrix.train.6.batch) <- row.names(pData.train.6)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.6.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.6$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 06 ######################################################
##############################################################################################################
fselect.6  <- genefilter(matrix.train.6.batch, filterfun(f1))
summary(fselect.6)
matrix.train.6 <-matrix.train.6.batch[fselect.6,]

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

system.time(rfe.6  <- rfe(matrix.train.6, labels.train.6, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.6,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.6  # 34 optVars found 
write.table(rfe.6$results, file = "TrainingSet06_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.6, type = c("g", "o"))
plot(rfe.6, type = c("g", "o"), xlim = c(0,61))

optFeatures.rfe.6 <- cbind(rfe.6$optVariables, Annotation[rfe.6$optVariables,])
write.table(optFeatures.rfe.6, file = "TrainingSet06_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 06 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets
rfeResamples.6 <- resamples(list("SVM_full.6" = svmFull.6,"SVM_RFE.6" = rfe.6))
sink("TrainingSet06_resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.6)
sink()

modelDifferences.6 <- diff(rfeResamples.6)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet06_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.6)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION ########################################
################################################################################################
matrix.train.rfe.6 <- matrix.train.6[,rfe.6$optVariables]   # subset fot the 33 optVars

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2210)
index.GA.6 <- createMultiFolds(labels.train.6, k=10, times = 5)  

outerctrl.GA.6 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.6,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  

#### 4.2 run GA  ###############################################################################
################################################################################################

system.time(GA.6<- gafs(matrix.train.rfe.6, labels.train.6, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.6,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.3 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.6
optVars.GA.6 <- Annotation[GA.6$optVariables,]  # export optimal variables
write.table(optVars.GA.6, file = "TrainingSet06_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.6 <- GA.6$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.6 <- arrange(performance.external.6, Iter)

performance.6 <- GA.6$ga$internal               # export internal accuracy (within the resample)
performance.6$AccuracyExternal <- aggregate(performance.external.6$Accuracy, by=list(Iter=performance.external.6$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.6 <- data.frame(Iter = c(performance.6$Iter,performance.6$Iter), Accuracy = c(performance.6$AccuracyExternal,performance.6$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.6 <- GA.6$averages[GA.6$optIter,2]
accuracy.external.6

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.6 <- subset(performance.external.6,performance.external.6$Iter == GA.6$optIter)
accuracy.external.opt.6 <- accuracy.external.opt.6$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.6, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 06 ########################
##############################################################################################################
matrix.train.6.bapred.GA <- t(matrix.train.6.batch[GA.6$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.6  <- train(matrix.train.6.bapred.GA,labels.train.6,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.6)
svmOpt.GA.6  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.6.qn  <- log2(t(qunormaddon(qunorm.train.6, t(RAW_test.6$E))))   # use qunormtrain object to normalize test set
matrix.test.6.qn  <- avereps(matrix.test.6.qn, ID= RAW_test.6$genes$ProbeName)
boxplot(cbind(matrix.train.6.qn,matrix.test.6.qn),col=c(pData.train.6$IVIM_Color,pData.test.6$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.6.bapred  <- t(combatbaaddon(combat.train.6, t(matrix.test.6.qn), batch = as.factor(pData.test.6$Batch)))

plot(Rtsne(t(cbind(matrix.train.6.bapred,matrix.test.6.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.6$Design_Color,pData.test.6$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.6.bapred.GA   <- t(matrix.test.6.bapred[GA.6$optVariables,])
matrix.test.6.bapred.full <- t(matrix.test.6.bapred[colnames(matrix.train.6),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.6.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.6$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.6$Design), col=unique(pData.train.6$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.6.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.6$Design_Color,pData.test.6$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.6), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.6 <- predict(svmOpt.GA.6,matrix.test.6.bapred.GA, type = "prob")
Prediction_GA.6$Prediction_GA.6 <- ifelse(Prediction_GA.6$transforming>0.50,"transforming","untransforming")
Prediction_GA.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_GA.6)
write.table(Prediction_GA.6, file = paste("TestSet06_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.6 <- predict(svmFull.6, matrix.test.6.bapred.full, type = "prob")
Prediction_SVM_full.6$Prediction_SVM_full <- ifelse(Prediction_SVM_full.6$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_SVM_full.6)
write.table(Prediction_SVM_full.6, file = paste("TestSet06_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 06 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet06_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.6$Prediction_GA.6), as.factor(Prediction_GA.6$TrueLabel))
sink()

sink("TestSet06_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.6$Prediction_SVM_full), as.factor(Prediction_SVM_full.6$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 06 ########################################################

Prediction_GA.6$Class <- as.factor(ifelse(Prediction_GA.6$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.6 <- roc(Prediction_GA.6$Class,                    # response vector (factor or character)
                Prediction_GA.6$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.6$Class <- as.factor(ifelse(Prediction_SVM_full.6$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.6<- roc(Prediction_SVM_full.6$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.6$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 06 ####################################
Prediction_GA.6$Class_Code <- ifelse(Prediction_GA.6$TrueLabel == "transforming",1,0)
pr.GA.6 <- pr.curve(scores.class0 = Prediction_GA.6$transforming , weights.class0 = Prediction_GA.6$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.6, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.6$Class_Code <- ifelse(Prediction_SVM_full.6$TrueLabel == "transforming",1,0)
pr.full.6 <- pr.curve(scores.class0 = Prediction_SVM_full.6$transforming , weights.class0 = Prediction_SVM_full.6$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.6, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 07 = Batch 07 = IVIM #160525 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 07 and TestSet 07 ############################################################
##############################################################################################################
pData.test.7  <- pData[pData$Batch==7,]    # set aside batch 7 / IVIM #160525
pData.train.7 <- pData[pData$Batch!=7,]    # use all remaining assays as training set
RAW_train.7   <- SAGA_RAW[,row.names(pData.train.7)]
RAW_test.7    <- SAGA_RAW[,row.names(pData.test.7)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.7$E), col=pData.train.7$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.7 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.7    <- qunormtrain(t(RAW_train.7$E))
matrix.train.7.qn <- log2(t(qunorm.train.7$xnorm))
matrix.train.7.qn <- avereps(matrix.train.7.qn, ID= RAW_train.7$genes$ProbeName)  
colnames(matrix.train.7.qn) <- row.names(pData.train.7)
boxplot(matrix.train.7.qn, col=pData.train.7$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.7.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.7$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.7         <- as.factor(ifelse(pData.train.7$Batch>7,pData.train.7$Batch-1,pData.train.7$Batch))                       
combat.train.7        <- combatba(t(matrix.train.7.qn), batch = batch.train.7)
matrix.train.7.batch  <- t(combat.train.7$xadj)
colnames(matrix.train.7.batch) <- row.names(pData.train.7)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.7.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.7$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 07 ######################################################
##############################################################################################################
fselect.7  <- genefilter(matrix.train.7.batch, filterfun(f1))
summary(fselect.7)
matrix.train.7 <-matrix.train.7.batch[fselect.7,]

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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.7    # 34 optVars found by rfe
write.table(rfe.7$results, file = "TrainingSet07_Results_rfe.txt", sep="\t",col.names=NA)


### Figure 3b ################################################################################################
trellis.par.set(caretTheme())
plot(rfe.7, type = c("g", "o"))
plot(rfe.7, type = c("g", "o"), xlim = c(0,61))    # max at 26 optVars found by rfe

optFeatures.7 <- cbind(rfe.7$optVariables, Annotation[rfe.7$optVariables,])
write.table(optFeatures.7, file = "TrainingSet07_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 07 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.7 <- resamples(list("SVM_full.7" = svmFull.7,"SVM_RFE.7" = rfe.7))
sink("TrainingSet07_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.7)
sink()

modelDifferences.7 <- diff(rfeResamples.7)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet07_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.7)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 07 #########################
################################################################################################
matrix.train.rfe.7 <- matrix.train.7[,rfe.7$optVariables]   # subset fot the 33 optVars

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2155)
index.GA.7 <- createMultiFolds(labels.train.7, k=10, times = 5)  

outerctrl.GA.7 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.7,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  

#### 4.2 run GA ################################################################################
################################################################################################
system.time(GA.7<- gafs(matrix.train.rfe.7, labels.train.7, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.7,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.3 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.7    # 15 optVars found by GA
optVars.GA.7 <- Annotation[GA.7$optVariables,]  # export optimal variables
write.table(optVars.GA.7, file = "TrainingSet07_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.7 <- GA.7$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.7 <- arrange(performance.external.7, Iter)

performance.7 <- GA.7$ga$internal               # export internal accuracy (within the resample)
performance.7$AccuracyExternal <- aggregate(performance.external.7$Accuracy, by=list(Iter=performance.external.7$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.7 <- data.frame(Iter = c(performance.7$Iter,performance.7$Iter), Accuracy = c(performance.7$AccuracyExternal,performance.7$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.7 <- GA.7$averages[GA.7$optIter,2]
accuracy.external.7

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.7 <- subset(performance.external.7,performance.external.7$Iter == GA.7$optIter)
accuracy.external.opt.7 <- accuracy.external.opt.7$Accuracy  

# plot internal and external accuracy over the iterations 
ggplot(performance.long.7, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())


##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 07 ########################
##############################################################################################################
matrix.train.7.bapred.GA <- t(matrix.train.7.batch[GA.7$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.7  <- train(matrix.train.7.bapred.GA,labels.train.7,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.7)
svmOpt.GA.7  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.7.qn  <- log2(t(qunormaddon(qunorm.train.7, t(RAW_test.7$E))))   # use qunormtrain object to normalize test set
matrix.test.7.qn  <- avereps(matrix.test.7.qn, ID= RAW_test.7$genes$ProbeName)
boxplot(cbind(matrix.train.7.qn,matrix.test.7.qn),col=c(pData.train.7$IVIM_Color,pData.test.7$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.7.bapred  <- t(combatbaaddon(combat.train.7, t(matrix.test.7.qn), batch = as.factor(pData.test.7$Batch)))

plot(Rtsne(t(cbind(matrix.train.7.bapred,matrix.test.7.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.7$Design_Color,pData.test.7$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.7.bapred.GA   <- t(matrix.test.7.bapred[GA.7$optVariables,])
matrix.test.7.bapred.full <- t(matrix.test.7.bapred[colnames(matrix.train.7),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.7.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.7$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.7$Design), col=unique(pData.train.7$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.7.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.7$Design_Color,pData.test.7$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.7), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.7 <- predict(svmOpt.GA.7,matrix.test.7.bapred.GA, type = "prob")
Prediction_GA.7$Prediction_GA.7 <- ifelse(Prediction_GA.7$transforming>0.50,"transforming","untransforming")
Prediction_GA.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_GA.7)
write.table(Prediction_GA.7, file = paste("TestSet07_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.7 <- predict(svmFull.7, matrix.test.7.bapred.full, type = "prob")
Prediction_SVM_full.7$Prediction_SVM_full <- ifelse(Prediction_SVM_full.7$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_SVM_full.7)
write.table(Prediction_SVM_full.7, file = paste("TestSet07_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 07 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet07_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.7$Prediction_GA.7), as.factor(Prediction_GA.7$TrueLabel))
sink()

sink("TestSet07_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.7$Prediction_SVM_full), as.factor(Prediction_SVM_full.7$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 07 ########################################################

Prediction_GA.7$Class <- as.factor(ifelse(Prediction_GA.7$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.7 <- roc(Prediction_GA.7$Class,                    # response vector (factor or character)
                Prediction_GA.7$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.7$Class <- as.factor(ifelse(Prediction_SVM_full.7$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.7<- roc(Prediction_SVM_full.7$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.7$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 07 ####################################
Prediction_GA.7$Class_Code <- ifelse(Prediction_GA.7$TrueLabel == "transforming",1,0)
pr.GA.7 <- pr.curve(scores.class0 = Prediction_GA.7$transforming , weights.class0 = Prediction_GA.7$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.7, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.7$Class_Code <- ifelse(Prediction_SVM_full.7$TrueLabel == "transforming",1,0)
pr.full.7 <- pr.curve(scores.class0 = Prediction_SVM_full.7$transforming , weights.class0 = Prediction_SVM_full.7$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.7, rand.plot = TRUE, legend = F, color = 1,main = "")


#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 08 = Batch 08 = IVIM #150318 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 08 and TestSet 08 ############################################################
##############################################################################################################
pData.test.8  <- pData[pData$Batch==8,]    # set aside batch 8 / IVIM #150318
pData.train.8 <- pData[pData$Batch!=8,]    # use all remaining assays as training set
RAW_train.8   <- SAGA_RAW[,row.names(pData.train.8)]
RAW_test.8    <- SAGA_RAW[,row.names(pData.test.8)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.8$E), col=pData.train.8$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.1 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.8    <- qunormtrain(t(RAW_train.8$E))
matrix.train.8.qn <- log2(t(qunorm.train.8$xnorm))
matrix.train.8.qn <- avereps(matrix.train.8.qn, ID= RAW_train.8$genes$ProbeName)  
colnames(matrix.train.8.qn) <- row.names(pData.train.8)
boxplot(matrix.train.8.qn, col=pData.train.8$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.8.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.8$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.8         <- as.factor(ifelse(pData.train.8$Batch>8,pData.train.8$Batch-1,pData.train.8$Batch))                       
combat.train.8        <- combatba(t(matrix.train.8.qn), batch = batch.train.8)
matrix.train.8.batch  <- t(combat.train.8$xadj)
colnames(matrix.train.8.batch) <- row.names(pData.train.8)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.8.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.8$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 08 ######################################################
##############################################################################################################
fselect.8  <- genefilter(matrix.train.8.batch, filterfun(f1))
summary(fselect.8)
matrix.train.8 <-matrix.train.8.batch[fselect.8,]

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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.8   # 23 variables found by SVM-rfe
write.table(rfe.8$results, file = "TrainingSet08_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.8, type = c("g", "o"))
plot(rfe.8, type = c("g", "o"), xlim = c(0,61))

optFeatures.8 <- cbind(rfe.8$optVariables, Annotation[rfe.8$optVariables,])
write.table(optFeatures.8, file = "TrainingSet08_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 08 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.8 <- resamples(list("SVM_full.8" = svmFull.8,"SVM_RFE.8" = rfe.8))
sink("TrainingSet08_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.8)
sink()

modelDifferences.8 <- diff(rfeResamples.8)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet08_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.8)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 08 #########################
################################################################################################

matrix.train.rfe.8 <- matrix.train.8[,rfe.8$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.8 <- createMultiFolds(labels.train.8, k=10, times = 5)  

outerctrl.GA.8 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.8,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  


system.time(GA.8<- gafs(matrix.train.rfe.8, labels.train.8, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.8,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.8   # yields 10 features
optVars.GA.8 <- Annotation[GA.8$optVariables,]  # export optimal variables
write.table(optVars.GA.8, file = "TrainingSet08_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.8 <- GA.8$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.8 <- arrange(performance.external.8, Iter)

performance.8 <- GA.8$ga$internal               # export internal accuracy (within the resample)
performance.8$AccuracyExternal <- aggregate(performance.external.8$Accuracy, by=list(Iter=performance.external.8$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.8 <- data.frame(Iter = c(performance.8$Iter,performance.8$Iter), Accuracy = c(performance.8$AccuracyExternal,performance.8$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.8 <- GA.8$averages[GA.8$optIter,2]
accuracy.external.8

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.8 <- subset(performance.external.8,performance.external.8$Iter == GA.8$optIter)
accuracy.external.opt.8 <- accuracy.external.opt.8$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.8, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 08 ########################
##############################################################################################################
matrix.train.8.bapred.GA <- t(matrix.train.8.batch[GA.8$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.8  <- train(matrix.train.8.bapred.GA,labels.train.8,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.8)
svmOpt.GA.8  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.8.qn  <- log2(t(qunormaddon(qunorm.train.8, t(RAW_test.8$E))))   # use qunormtrain object to normalize test set
matrix.test.8.qn  <- avereps(matrix.test.8.qn, ID= RAW_test.8$genes$ProbeName)
boxplot(cbind(matrix.train.8.qn,matrix.test.8.qn),col=c(pData.train.8$IVIM_Color,pData.test.8$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.8.bapred  <- t(combatbaaddon(combat.train.8, t(matrix.test.8.qn), batch = as.factor(pData.test.8$Batch)))

plot(Rtsne(t(cbind(matrix.train.8.bapred,matrix.test.8.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.8$Design_Color,pData.test.8$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.8.bapred.GA   <- t(matrix.test.8.bapred[GA.8$optVariables,])
matrix.test.8.bapred.full <- t(matrix.test.8.bapred[colnames(matrix.train.8),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.8.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.8$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.8$Design), col=unique(pData.train.8$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.8.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.8$Design_Color,pData.test.8$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.8), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.8 <- predict(svmOpt.GA.8,matrix.test.8.bapred.GA, type = "prob")
Prediction_GA.8$Prediction_GA.8 <- ifelse(Prediction_GA.8$transforming>0.50,"transforming","untransforming")
Prediction_GA.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_GA.8)
write.table(Prediction_GA.8, file = paste("TestSet08_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.8 <- predict(svmFull.8, matrix.test.8.bapred.full, type = "prob")
Prediction_SVM_full.8$Prediction_SVM_full <- ifelse(Prediction_SVM_full.8$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_SVM_full.8)
write.table(Prediction_SVM_full.8, file = paste("TestSet08_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 08 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet08_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.8$Prediction_GA.8), as.factor(Prediction_GA.8$TrueLabel))
sink()

sink("TestSet08_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.8$Prediction_SVM_full), as.factor(Prediction_SVM_full.8$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 08 ########################################################

Prediction_GA.8$Class <- as.factor(ifelse(Prediction_GA.8$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.8 <- roc(Prediction_GA.8$Class,                    # response vector (factor or character)
                Prediction_GA.8$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.8$Class <- as.factor(ifelse(Prediction_SVM_full.8$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.8<- roc(Prediction_SVM_full.8$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.8$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 08 ####################################
Prediction_GA.8$Class_Code <- ifelse(Prediction_GA.8$TrueLabel == "transforming",1,0)
pr.GA.8 <- pr.curve(scores.class0 = Prediction_GA.8$transforming , weights.class0 = Prediction_GA.8$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.8, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.8$Class_Code <- ifelse(Prediction_SVM_full.8$TrueLabel == "transforming",1,0)
pr.full.8 <- pr.curve(scores.class0 = Prediction_SVM_full.8$transforming , weights.class0 = Prediction_SVM_full.8$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.8, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 09 = Batch 09 = IVIM #161102 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 09 and TestSet 09 ############################################################
##############################################################################################################
pData.test.9  <- pData[pData$Batch==9,]    # set aside batch 9 / IVIM #161102
pData.train.9 <- pData[pData$Batch!=9,]    # use all remaining assays as training set
RAW_train.9   <- SAGA_RAW[,row.names(pData.train.9)]
RAW_test.9    <- SAGA_RAW[,row.names(pData.test.9)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.9$E), col=pData.train.9$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.9 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.9    <- qunormtrain(t(RAW_train.9$E))
matrix.train.9.qn <- log2(t(qunorm.train.9$xnorm))
matrix.train.9.qn <- avereps(matrix.train.9.qn, ID= RAW_train.9$genes$ProbeName)  
colnames(matrix.train.9.qn) <- row.names(pData.train.9)
boxplot(matrix.train.9.qn, col=pData.train.9$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.9.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.9$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.9         <- as.factor(ifelse(pData.train.9$Batch>9,pData.train.9$Batch-1,pData.train.9$Batch))                       
combat.train.9        <- combatba(t(matrix.train.9.qn), batch = batch.train.9)
matrix.train.9.batch  <- t(combat.train.9$xadj)
colnames(matrix.train.9.batch) <- row.names(pData.train.9)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.9.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.9$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 09 ######################################################
##############################################################################################################
fselect.9  <- genefilter(matrix.train.9.batch, filterfun(f1))
summary(fselect.9)
matrix.train.9 <-matrix.train.9.batch[fselect.9,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 09 #######################################################################
##############################################################################################################
matrix.train.9 <- (t(matrix.train.9))
labels.train.9 <- as.factor(pData.train.9$Class)

set.seed(2351)
index.9 <- createMultiFolds(labels.train.9, k=10, times = 20)  

fullCtrl.9 <- trainControl(method = "repeatedcv",repeats = 20,
                           index = index.9,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           allowParallel = TRUE)

set.seed(123)
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

system.time(rfe.9  <- rfe(matrix.train.9, labels.train.9, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.9,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.9 # 50 optVars found by rfe
write.table(rfe.9$results, file = "TrainingSet09_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.9, type = c("g", "o"))
plot(rfe.9, type = c("g", "o"), xlim = c(0,61))

optFeatures.9 <- cbind(rfe.9$optVariables, Annotation[rfe.9$optVariables,])
write.table(optFeatures.9, file = "TrainingSet09_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 09 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.9 <- resamples(list("SVM_full.9" = svmFull.9,"SVM_RFE.9" = rfe.9))
sink("TrainingSet09_Resamples_rfe_vs_full.txt", append = TRUE)
summary(rfeResamples.9)
sink()

modelDifferences.9 <- diff(rfeResamples.9)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet09_ModelDifferences_rfe_vs_full.txt", append = TRUE)
summary(modelDifferences.9)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 09 #########################
################################################################################################

matrix.train.rfe.9 <- matrix.train.9[,rfe.9$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.9 <- createMultiFolds(labels.train.9, k=10, times = 5)  

outerctrl.GA.9 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.9,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  


system.time(GA.9<- gafs(matrix.train.rfe.9, labels.train.9, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.9,
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.9   # yields 10 features
optVars.GA.9 <- Annotation[GA.9$optVariables,]  # export optimal variables
write.table(optVars.GA.9, file = "TrainingSet09_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.9 <- GA.9$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.9 <- arrange(performance.external.9, Iter)

performance.9 <- GA.9$ga$internal               # export internal accuracy (within the resample)
performance.9$AccuracyExternal <- aggregate(performance.external.9$Accuracy, by=list(Iter=performance.external.9$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.9 <- data.frame(Iter = c(performance.9$Iter,performance.9$Iter), Accuracy = c(performance.9$AccuracyExternal,performance.9$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.9 <- GA.9$averages[GA.9$optIter,2]
accuracy.external.9

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.9 <- subset(performance.external.9,performance.external.9$Iter == GA.9$optIter)
accuracy.external.opt.9 <- accuracy.external.opt.9$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.9, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 09 ########################
##############################################################################################################
matrix.train.9.bapred.GA <- t(matrix.train.9.batch[GA.9$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.9  <- train(matrix.train.9.bapred.GA,labels.train.9,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.9)
svmOpt.GA.9  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.9.qn  <- log2(t(qunormaddon(qunorm.train.9, t(RAW_test.9$E))))   # use qunormtrain object to normalize test set
matrix.test.9.qn  <- avereps(matrix.test.9.qn, ID= RAW_test.9$genes$ProbeName)
boxplot(cbind(matrix.train.9.qn,matrix.test.9.qn),col=c(pData.train.9$IVIM_Color,pData.test.9$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.9.bapred  <- t(combatbaaddon(combat.train.9, t(matrix.test.9.qn), batch = as.factor(pData.test.9$Batch)))

plot(Rtsne(t(cbind(matrix.train.9.bapred,matrix.test.9.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.9$Design_Color,pData.test.9$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.9.bapred.GA   <- t(matrix.test.9.bapred[GA.9$optVariables,])
matrix.test.9.bapred.full <- t(matrix.test.9.bapred[colnames(matrix.train.9),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.9.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.9$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.9$Design), col=unique(pData.train.9$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.9.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.9$Design_Color,pData.test.9$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.9), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.9 <- predict(svmOpt.GA.9,matrix.test.9.bapred.GA, type = "prob")
Prediction_GA.9$Prediction_GA.9 <- ifelse(Prediction_GA.9$transforming>0.50,"transforming","untransforming")
Prediction_GA.9 <- cbind(pData.test.9[,c(1:3)],TrueLabel=pData.test.9$Class,Prediction_GA.9)
write.table(Prediction_GA.9, file = paste("TestSet09_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.9 <- predict(svmFull.9, matrix.test.9.bapred.full, type = "prob")
Prediction_SVM_full.9$Prediction_SVM_full <- ifelse(Prediction_SVM_full.9$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.9 <- cbind(pData.test.9[,c(1:3)],TrueLabel=pData.test.9$Class,Prediction_SVM_full.9)
write.table(Prediction_SVM_full.9, file = paste("TestSet09_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 09 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet09_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.9$Prediction_GA.9), as.factor(Prediction_GA.9$TrueLabel))
sink()

sink("TestSet09_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.9$Prediction_SVM_full), as.factor(Prediction_SVM_full.9$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 09 ########################################################

Prediction_GA.9$Class <- as.factor(ifelse(Prediction_GA.9$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.9 <- roc(Prediction_GA.9$Class,                    # response vector (factor or character)
                Prediction_GA.9$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)

Prediction_SVM_full.9$Class <- as.factor(ifelse(Prediction_SVM_full.9$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.9<- roc(Prediction_SVM_full.9$Class,                    # response vector (factor or character)
                 Prediction_SVM_full.9$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 09 ####################################
Prediction_GA.9$Class_Code <- ifelse(Prediction_GA.9$TrueLabel == "transforming",1,0)
pr.GA.9 <- pr.curve(scores.class0 = Prediction_GA.9$transforming , weights.class0 = Prediction_GA.9$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.9, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.9$Class_Code <- ifelse(Prediction_SVM_full.9$TrueLabel == "transforming",1,0)
pr.full.9 <- pr.curve(scores.class0 = Prediction_SVM_full.9$transforming , weights.class0 = Prediction_SVM_full.9$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.9, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 010 = Batch 010 = IVIM #170125 #############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 010 and TestSet 010 ##########################################################
##############################################################################################################
pData.test.10  <- pData[pData$Batch==10,]    # set aside batch 10 / IVIM #160210
pData.train.10 <- pData[pData$Batch!=10,]    # use all remaining assays as training set
RAW_train.10   <- SAGA_RAW[,row.names(pData.train.10)]
RAW_test.10    <- SAGA_RAW[,row.names(pData.test.10)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.10$E), col=pData.train.10$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.1 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.10    <- qunormtrain(t(RAW_train.10$E))
matrix.train.10.qn <- log2(t(qunorm.train.10$xnorm))
matrix.train.10.qn <- avereps(matrix.train.10.qn, ID= RAW_train.10$genes$ProbeName)  
colnames(matrix.train.10.qn) <- row.names(pData.train.10)
boxplot(matrix.train.10.qn, col=pData.train.10$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.10.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.10$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.10         <- as.factor(ifelse(pData.train.10$Batch>10,pData.train.10$Batch-1,pData.train.10$Batch))                       
combat.train.10        <- combatba(t(matrix.train.10.qn), batch = batch.train.10)
matrix.train.10.batch  <- t(combat.train.10$xadj)
colnames(matrix.train.10.batch) <- row.names(pData.train.10)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.10.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.10$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 010 ######################################################
##############################################################################################################
fselect.10  <- genefilter(matrix.train.10.batch, filterfun(f1))
summary(fselect.10)
matrix.train.10 <-matrix.train.10.batch[fselect.10,]

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

system.time(rfe.10  <- rfe(matrix.train.10, labels.train.10, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.10,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.10   # 18 optVars selected
write.table(rfe.10$results, file = "TrainingSet10_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.10, type = c("g", "o"))
plot(rfe.10, type = c("g", "o"), xlim = c(0,61))

optFeatures.10 <- cbind(rfe.10$optVariables, Annotation[rfe.10$optVariables,])
write.table(optFeatures.10, file = "TrainingSet10_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 10  #####################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (2008) since performances were measured on identically resampled sets

rfeResamples.10 <- resamples(list("SVM_full.10" = svmFull.10,"SVM_RFE.10" = rfe.10))
sink("TrainingSet10_Resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.10)
sink()

modelDifferences.10 <- diff(rfeResamples.10)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet10_ModelDifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.10)
sink()

################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 10 #########################
################################################################################################
matrix.train.rfe.10 <- matrix.train.10[,rfe.10$optVariables]   # subset fot the 45 optVars

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2312)
index.GA.10 <- createMultiFolds(labels.train.10, k=10, times = 5)  

outerctrl.GA.10 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.10,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  

#### 4.2 run GA  ###############################################################################
################################################################################################
system.time(GA.10<- gafs(matrix.train.rfe.10, labels.train.10, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.10,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.3 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.10   # 19 features selected
optVars.GA.10 <- Annotation[GA.10$optVariables,]  # export optimal variables
write.table(optVars.GA.10, file = "TrainingSet10_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.10 <- GA.10$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.10 <- arrange(performance.external.10, Iter)

performance.10 <- GA.10$ga$internal               # export internal accuracy (within the resample)
performance.10$AccuracyExternal <- aggregate(performance.external.10$Accuracy, by=list(Iter=performance.external.10$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.10 <- data.frame(Iter = c(performance.10$Iter,performance.10$Iter), Accuracy = c(performance.10$AccuracyExternal,performance.10$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.10 <- GA.10$averages[GA.10$optIter,2]
accuracy.external.10

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.10 <- subset(performance.external.10,performance.external.10$Iter == GA.10$optIter)
accuracy.external.opt.10 <- accuracy.external.opt.10$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.10, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.5,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 010 #######################
##############################################################################################################
matrix.train.10.bapred.GA <- t(matrix.train.10.batch[GA.10$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.10  <- train(matrix.train.10.bapred.GA,labels.train.10,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.10)
svmOpt.GA.10  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.10.qn  <- log2(t(qunormaddon(qunorm.train.10, t(RAW_test.10$E))))   # use qunormtrain object to normalize test set
matrix.test.10.qn  <- avereps(matrix.test.10.qn, ID= RAW_test.10$genes$ProbeName)
boxplot(cbind(matrix.train.10.qn,matrix.test.10.qn),col=c(pData.train.10$IVIM_Color,pData.test.10$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.10.bapred  <- t(combatbaaddon(combat.train.10, t(matrix.test.10.qn), batch = as.factor(pData.test.10$Batch)))

plot(Rtsne(t(cbind(matrix.train.10.bapred,matrix.test.10.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.10$Design_Color,pData.test.10$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.10.bapred.GA   <- t(matrix.test.10.bapred[GA.10$optVariables,])
matrix.test.10.bapred.full <- t(matrix.test.10.bapred[colnames(matrix.train.10),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.10.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.10$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.10$Design), col=unique(pData.train.10$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.10.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.10$Design_Color,pData.test.10$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.10), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.10 <- predict(svmOpt.GA.10,matrix.test.10.bapred.GA, type = "prob")
Prediction_GA.10$Prediction_GA.10 <- ifelse(Prediction_GA.10$transforming>0.50,"transforming","untransforming")
Prediction_GA.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_GA.10)
write.table(Prediction_GA.10, file = paste("TestSet010_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.10 <- predict(svmFull.10, matrix.test.10.bapred.full, type = "prob")
Prediction_SVM_full.10$Prediction_SVM_full <- ifelse(Prediction_SVM_full.10$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_SVM_full.10)
write.table(Prediction_SVM_full.10, file = paste("TestSet010_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 010 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet010_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.10$Prediction_GA.10), as.factor(Prediction_GA.10$TrueLabel))
sink()

sink("TestSet010_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.10$Prediction_SVM_full), as.factor(Prediction_SVM_full.10$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 010 ########################################################

Prediction_GA.10$Class <- as.factor(ifelse(Prediction_GA.10$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.10 <- roc(Prediction_GA.10$Class,                    # response vector (factor or character)
                 Prediction_GA.10$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.10$Class <- as.factor(ifelse(Prediction_SVM_full.10$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.10<- roc(Prediction_SVM_full.10$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.10$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 010 ####################################
Prediction_GA.10$Class_Code <- ifelse(Prediction_GA.10$TrueLabel == "transforming",1,0)
pr.GA.10 <- pr.curve(scores.class0 = Prediction_GA.10$transforming , weights.class0 = Prediction_GA.10$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.10, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.10$Class_Code <- ifelse(Prediction_SVM_full.10$TrueLabel == "transforming",1,0)
pr.full.10 <- pr.curve(scores.class0 = Prediction_SVM_full.10$transforming , weights.class0 = Prediction_SVM_full.10$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.10, rand.plot = TRUE, legend = F, color = 1,main = "")


#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 011 = Batch 011 = IVIM #171102 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 011 and TestSet 011 ############################################################
##############################################################################################################
pData.test.11  <- pData[pData$Batch==11,]    # set aside batch 11 / IVIM #171102
pData.train.11 <- pData[pData$Batch!=11,]    # use all remaining assays as training set
RAW_train.11   <- SAGA_RAW[,row.names(pData.train.11)]
RAW_test.11    <- SAGA_RAW[,row.names(pData.test.11)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.11$E), col=pData.train.11$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.11 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.11    <- qunormtrain(t(RAW_train.11$E))
matrix.train.11.qn <- log2(t(qunorm.train.11$xnorm))
matrix.train.11.qn <- avereps(matrix.train.11.qn, ID= RAW_train.11$genes$ProbeName)  
colnames(matrix.train.11.qn) <- row.names(pData.train.11)
boxplot(matrix.train.11.qn, col=pData.train.11$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.11.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.11$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.11         <- as.factor(ifelse(pData.train.11$Batch>11,pData.train.11$Batch-1,pData.train.11$Batch))                       
combat.train.11        <- combatba(t(matrix.train.11.qn), batch = batch.train.11)
matrix.train.11.batch  <- t(combat.train.11$xadj)
colnames(matrix.train.11.batch) <- row.names(pData.train.11)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.11.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.11$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 011 ######################################################
##############################################################################################################
fselect.11  <- genefilter(matrix.train.11.batch, filterfun(f1))
summary(fselect.11)
matrix.train.11 <-matrix.train.11.batch[fselect.11,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 011 #######################################################################
##############################################################################################################
matrix.train.11 <- (t(matrix.train.11))
labels.train.11 <- as.factor(pData.train.11$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.11 <- createMultiFolds(labels.train.11, k=10, times = 20)  

fullCtrl.11 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.11,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.11 <- train(matrix.train.11,labels.train.11,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.11)

svmFull.11  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 011 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.11      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.11,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.11$functions         <- caretFuncs
outerctrl.11$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 011 ##############################################################################
##############################################################################################################

system.time(rfe.11  <- rfe(matrix.train.11, labels.train.11, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.11,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.11   # 20 variables found by SVM-rfe
write.table(rfe.11$results, file = "TrainingSet011_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.11, type = c("g", "o"))
plot(rfe.11, type = c("g", "o"), xlim = c(0,61))

optFeatures.11 <- cbind(rfe.11$optVariables, Annotation[rfe.11$optVariables,])
write.table(optFeatures.11, file = "TrainingSet011_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 011 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20011) since performances were measured on identically resampled sets

rfeResamples.11 <- resamples(list("SVM_full.11" = svmFull.11,"SVM_RFE.11" = rfe.11))
sink("TrainingSet011_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.11)
sink()

modelDifferences.11 <- diff(rfeResamples.11)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet011_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.11)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 011 #########################
################################################################################################

matrix.train.rfe.11 <- matrix.train.11[,rfe.11$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.11 <- createMultiFolds(labels.train.11, k=10, times = 5)  

outerctrl.GA.11 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.11,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.11<- gafs(matrix.train.rfe.11, labels.train.11, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.11,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.11  # yields 9 features
optVars.GA.11 <- Annotation[GA.11$optVariables,]  # export optimal variables
write.table(optVars.GA.11, file = "TrainingSet011_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.11 <- GA.11$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.11 <- arrange(performance.external.11, Iter)

performance.11 <- GA.11$ga$internal               # export internal accuracy (within the resample)
performance.11$AccuracyExternal <- aggregate(performance.external.11$Accuracy, by=list(Iter=performance.external.11$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.11 <- data.frame(Iter = c(performance.11$Iter,performance.11$Iter), Accuracy = c(performance.11$AccuracyExternal,performance.11$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.11 <- GA.11$averages[GA.11$optIter,2]
accuracy.external.11

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.11 <- subset(performance.external.11,performance.external.11$Iter == GA.11$optIter)
accuracy.external.opt.11 <- accuracy.external.opt.11$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.11, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 011 ########################
##############################################################################################################
matrix.train.11.bapred.GA <- t(matrix.train.11.batch[GA.11$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.11  <- train(matrix.train.11.bapred.GA,labels.train.11,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.11)
svmOpt.GA.11  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.11.qn  <- log2(t(qunormaddon(qunorm.train.11, t(RAW_test.11$E))))   # use qunormtrain object to normalize test set
matrix.test.11.qn  <- avereps(matrix.test.11.qn, ID= RAW_test.11$genes$ProbeName)
boxplot(cbind(matrix.train.11.qn,matrix.test.11.qn),col=c(pData.train.11$IVIM_Color,pData.test.11$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.11.bapred  <- t(combatbaaddon(combat.train.11, t(matrix.test.11.qn), batch = as.factor(pData.test.11$Batch)))

plot(Rtsne(t(cbind(matrix.train.11.bapred,matrix.test.11.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.11$Design_Color,pData.test.11$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.11.bapred.GA   <- t(matrix.test.11.bapred[GA.11$optVariables,])
matrix.test.11.bapred.full <- t(matrix.test.11.bapred[colnames(matrix.train.11),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.11.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.11$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.11$Design), col=unique(pData.train.11$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.11.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.11$Design_Color,pData.test.11$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.11), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.11 <- predict(svmOpt.GA.11,matrix.test.11.bapred.GA, type = "prob")
Prediction_GA.11$Prediction_GA.11 <- ifelse(Prediction_GA.11$transforming>0.50,"transforming","untransforming")
Prediction_GA.11 <- cbind(pData.test.11[,c(1:3)],TrueLabel=pData.test.11$Class,Prediction_GA.11)
write.table(Prediction_GA.11, file = paste("TestSet011_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.11 <- predict(svmFull.11, matrix.test.11.bapred.full, type = "prob")
Prediction_SVM_full.11$Prediction_SVM_full <- ifelse(Prediction_SVM_full.11$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.11 <- cbind(pData.test.11[,c(1:3)],TrueLabel=pData.test.11$Class,Prediction_SVM_full.11)
write.table(Prediction_SVM_full.11, file = paste("TestSet011_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 011 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet011_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.11$Prediction_GA.11), as.factor(Prediction_GA.11$TrueLabel))
sink()

sink("TestSet011_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.11$Prediction_SVM_full), as.factor(Prediction_SVM_full.11$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 011 ########################################################

Prediction_GA.11$Class <- as.factor(ifelse(Prediction_GA.11$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.11 <- roc(Prediction_GA.11$Class,                    # response vector (factor or character)
                 Prediction_GA.11$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.11$Class <- as.factor(ifelse(Prediction_SVM_full.11$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.11<- roc(Prediction_SVM_full.11$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.11$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 011 ####################################
Prediction_GA.11$Class_Code <- ifelse(Prediction_GA.11$TrueLabel == "transforming",1,0)
pr.GA.11 <- pr.curve(scores.class0 = Prediction_GA.11$transforming , weights.class0 = Prediction_GA.11$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.11, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.11$Class_Code <- ifelse(Prediction_SVM_full.11$TrueLabel == "transforming",1,0)
pr.full.11 <- pr.curve(scores.class0 = Prediction_SVM_full.11$transforming , weights.class0 = Prediction_SVM_full.11$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.11, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 012 = Batch 012 = IVIM #171115 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 012 and TestSet 012 ############################################################
##############################################################################################################
pData.test.12  <- pData[pData$Batch==12,]    # set aside batch 12 / IVIM #171115
pData.train.12 <- pData[pData$Batch!=12,]    # use all remaining assays as training set
RAW_train.12   <- SAGA_RAW[,row.names(pData.train.12)]
RAW_test.12    <- SAGA_RAW[,row.names(pData.test.12)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.12$E), col=pData.train.12$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.12 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.12    <- qunormtrain(t(RAW_train.12$E))
matrix.train.12.qn <- log2(t(qunorm.train.12$xnorm))
matrix.train.12.qn <- avereps(matrix.train.12.qn, ID= RAW_train.12$genes$ProbeName)  
colnames(matrix.train.12.qn) <- row.names(pData.train.12)
boxplot(matrix.train.12.qn, col=pData.train.12$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.12.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.12$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.12         <- as.factor(ifelse(pData.train.12$Batch>12,pData.train.12$Batch-1,pData.train.12$Batch))                       
combat.train.12        <- combatba(t(matrix.train.12.qn), batch = batch.train.12)
matrix.train.12.batch  <- t(combat.train.12$xadj)
colnames(matrix.train.12.batch) <- row.names(pData.train.12)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.12.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.12$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 012 ######################################################
##############################################################################################################
fselect.12  <- genefilter(matrix.train.12.batch, filterfun(f1))
summary(fselect.12)
matrix.train.12 <-matrix.train.12.batch[fselect.12,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 012 #######################################################################
##############################################################################################################
matrix.train.12 <- (t(matrix.train.12))
labels.train.12 <- as.factor(pData.train.12$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1657)
index.12 <- createMultiFolds(labels.train.12, k=10, times = 20)  

fullCtrl.12 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.12,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.12 <- train(matrix.train.12,labels.train.12,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.12)

svmFull.12  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 012 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.12      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.12,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.12$functions         <- caretFuncs
outerctrl.12$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 012 ##############################################################################
##############################################################################################################

system.time(rfe.12  <- rfe(matrix.train.12, labels.train.12, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.12,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.12   # 21 variables found by SVM-rfe
write.table(rfe.12$results, file = "TrainingSet012_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.12, type = c("g", "o"))
plot(rfe.12, type = c("g", "o"), xlim = c(0,61))

optFeatures.12 <- cbind(rfe.12$optVariables, Annotation[rfe.12$optVariables,])
write.table(optFeatures.12, file = "TrainingSet012_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 012 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20012) since performances were measured on identically resampled sets

rfeResamples.12 <- resamples(list("SVM_full.12" = svmFull.12,"SVM_RFE.12" = rfe.12))
sink("TrainingSet012_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.12)
sink()

modelDifferences.12 <- diff(rfeResamples.12)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet012_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.12)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 012 #########################
################################################################################################

matrix.train.rfe.12 <- matrix.train.12[,rfe.12$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.12 <- createMultiFolds(labels.train.12, k=10, times = 5)  

outerctrl.GA.12 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.12,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.12<- gafs(matrix.train.rfe.12, labels.train.12, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.12,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.12   # yields 10 features
optVars.GA.12 <- Annotation[GA.12$optVariables,]  # export optimal variables
write.table(optVars.GA.12, file = "TrainingSet012_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.12 <- GA.12$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.12 <- arrange(performance.external.12, Iter)

performance.12 <- GA.12$ga$internal               # export internal accuracy (within the resample)
performance.12$AccuracyExternal <- aggregate(performance.external.12$Accuracy, by=list(Iter=performance.external.12$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.12 <- data.frame(Iter = c(performance.12$Iter,performance.12$Iter), Accuracy = c(performance.12$AccuracyExternal,performance.12$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.12 <- GA.12$averages[GA.12$optIter,2]
accuracy.external.12

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.12 <- subset(performance.external.12,performance.external.12$Iter == GA.12$optIter)
accuracy.external.opt.12 <- accuracy.external.opt.12$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.12, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 012 #######################
##############################################################################################################
matrix.train.12.bapred.GA <- t(matrix.train.12.batch[GA.12$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.12  <- train(matrix.train.12.bapred.GA,labels.train.12,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.12)
svmOpt.GA.12  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.12.qn  <- log2(t(qunormaddon(qunorm.train.12, t(RAW_test.12$E))))   # use qunormtrain object to normalize test set
matrix.test.12.qn  <- avereps(matrix.test.12.qn, ID= RAW_test.12$genes$ProbeName)
boxplot(cbind(matrix.train.12.qn,matrix.test.12.qn),col=c(pData.train.12$IVIM_Color,pData.test.12$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.12.bapred  <- t(combatbaaddon(combat.train.12, t(matrix.test.12.qn), batch = as.factor(pData.test.12$Batch)))

plot(Rtsne(t(cbind(matrix.train.12.bapred,matrix.test.12.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.12$Design_Color,pData.test.12$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.12.bapred.GA   <- t(matrix.test.12.bapred[GA.12$optVariables,])
matrix.test.12.bapred.full <- t(matrix.test.12.bapred[colnames(matrix.train.12),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.12.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.12$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.12$Design), col=unique(pData.train.12$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.12.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.12$Design_Color,pData.test.12$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.12), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.12 <- predict(svmOpt.GA.12,matrix.test.12.bapred.GA, type = "prob")
Prediction_GA.12$Prediction_GA.12 <- ifelse(Prediction_GA.12$transforming>0.50,"transforming","untransforming")
Prediction_GA.12 <- cbind(pData.test.12[,c(1:3)],TrueLabel=pData.test.12$Class,Prediction_GA.12)
write.table(Prediction_GA.12, file = paste("TestSet012_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.12 <- predict(svmFull.12, matrix.test.12.bapred.full, type = "prob")
Prediction_SVM_full.12$Prediction_SVM_full <- ifelse(Prediction_SVM_full.12$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.12 <- cbind(pData.test.12[,c(1:3)],TrueLabel=pData.test.12$Class,Prediction_SVM_full.12)
write.table(Prediction_SVM_full.12, file = paste("TestSet012_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 012 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet012_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.12$Prediction_GA.12), as.factor(Prediction_GA.12$TrueLabel))
sink()

sink("TestSet012_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.12$Prediction_SVM_full), as.factor(Prediction_SVM_full.12$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 012 ########################################################

Prediction_GA.12$Class <- as.factor(ifelse(Prediction_GA.12$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.12 <- roc(Prediction_GA.12$Class,                    # response vector (factor or character)
                 Prediction_GA.12$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.12$Class <- as.factor(ifelse(Prediction_SVM_full.12$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.12<- roc(Prediction_SVM_full.12$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.12$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 012 ####################################
Prediction_GA.12$Class_Code <- ifelse(Prediction_GA.12$TrueLabel == "transforming",1,0)
pr.GA.12 <- pr.curve(scores.class0 = Prediction_GA.12$transforming , weights.class0 = Prediction_GA.12$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.12, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.12$Class_Code <- ifelse(Prediction_SVM_full.12$TrueLabel == "transforming",1,0)
pr.full.12 <- pr.curve(scores.class0 = Prediction_SVM_full.12$transforming , weights.class0 = Prediction_SVM_full.12$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.12, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 013 = Batch 013 = IVIM #180110 #############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 013 and TestSet 013 ############################################################
##############################################################################################################
pData.test.13  <- pData[pData$Batch==13,]    # set aside batch 13 / IVIM #180110
pData.train.13 <- pData[pData$Batch!=13,]    # use all remaining assays as training set
RAW_train.13   <- SAGA_RAW[,row.names(pData.train.13)]
RAW_test.13    <- SAGA_RAW[,row.names(pData.test.13)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.13$E), col=pData.train.13$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.13 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.13    <- qunormtrain(t(RAW_train.13$E))
matrix.train.13.qn <- log2(t(qunorm.train.13$xnorm))
matrix.train.13.qn <- avereps(matrix.train.13.qn, ID= RAW_train.13$genes$ProbeName)  
colnames(matrix.train.13.qn) <- row.names(pData.train.13)
boxplot(matrix.train.13.qn, col=pData.train.13$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.13.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.13$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.13         <- as.factor(ifelse(pData.train.13$Batch>13,pData.train.13$Batch-1,pData.train.13$Batch))                       
combat.train.13        <- combatba(t(matrix.train.13.qn), batch = batch.train.13)
matrix.train.13.batch  <- t(combat.train.13$xadj)
colnames(matrix.train.13.batch) <- row.names(pData.train.13)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.13.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.13$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 013 ######################################################
##############################################################################################################
fselect.13  <- genefilter(matrix.train.13.batch, filterfun(f1))
summary(fselect.13)
matrix.train.13 <-matrix.train.13.batch[fselect.13,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 013 #######################################################################
##############################################################################################################
matrix.train.13 <- (t(matrix.train.13))
labels.train.13 <- as.factor(pData.train.13$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.13 <- createMultiFolds(labels.train.13, k=10, times = 20)  

fullCtrl.13 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.13,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.13 <- train(matrix.train.13,labels.train.13,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.13)

svmFull.13  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 013 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.13      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.13,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.13$functions         <- caretFuncs
outerctrl.13$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 013 ##############################################################################
##############################################################################################################

system.time(rfe.13  <- rfe(matrix.train.13, labels.train.13, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.13,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.13   # 28 variables found by SVM-rfe
write.table(rfe.13$results, file = "TrainingSet013_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.13, type = c("g", "o"))
plot(rfe.13, type = c("g", "o"), xlim = c(0,61))

optFeatures.13 <- cbind(rfe.13$optVariables, Annotation[rfe.13$optVariables,])
write.table(optFeatures.13, file = "TrainingSet013_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 013 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20013) since performances were measured on identically resampled sets

rfeResamples.13 <- resamples(list("SVM_full.13" = svmFull.13,"SVM_RFE.13" = rfe.13))
sink("TrainingSet013_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.13)
sink()

modelDifferences.13 <- diff(rfeResamples.13)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet013_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.13)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 013 #########################
################################################################################################

matrix.train.rfe.13 <- matrix.train.13[,rfe.13$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2320)
index.GA.13 <- createMultiFolds(labels.train.13, k=10, times = 5)  

outerctrl.GA.13 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.13,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.13<- gafs(matrix.train.rfe.13, labels.train.13, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.13,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.13   # yields 12 features
optVars.GA.13 <- Annotation[GA.13$optVariables,]  # export optimal variables
write.table(optVars.GA.13, file = "TrainingSet013_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.13 <- GA.13$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.13 <- arrange(performance.external.13, Iter)

performance.13 <- GA.13$ga$internal               # export internal accuracy (within the resample)
performance.13$AccuracyExternal <- aggregate(performance.external.13$Accuracy, by=list(Iter=performance.external.13$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.13 <- data.frame(Iter = c(performance.13$Iter,performance.13$Iter), Accuracy = c(performance.13$AccuracyExternal,performance.13$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.13 <- GA.13$averages[GA.13$optIter,2]
accuracy.external.13

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.13 <- subset(performance.external.13,performance.external.13$Iter == GA.13$optIter)
accuracy.external.opt.13 <- accuracy.external.opt.13$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.13, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 013 ########################
##############################################################################################################
matrix.train.13.bapred.GA <- t(matrix.train.13.batch[GA.13$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.13  <- train(matrix.train.13.bapred.GA,labels.train.13,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.13)
svmOpt.GA.13  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.13.qn  <- log2(t(qunormaddon(qunorm.train.13, t(RAW_test.13$E))))   # use qunormtrain object to normalize test set
matrix.test.13.qn  <- avereps(matrix.test.13.qn, ID= RAW_test.13$genes$ProbeName)
boxplot(cbind(matrix.train.13.qn,matrix.test.13.qn),col=c(pData.train.13$IVIM_Color,pData.test.13$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.13.bapred  <- t(combatbaaddon(combat.train.13, t(matrix.test.13.qn), batch = as.factor(pData.test.13$Batch)))

plot(Rtsne(t(cbind(matrix.train.13.bapred,matrix.test.13.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.13$Design_Color,pData.test.13$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.13.bapred.GA   <- t(matrix.test.13.bapred[GA.13$optVariables,])
matrix.test.13.bapred.full <- t(matrix.test.13.bapred[colnames(matrix.train.13),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.13.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.13$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.13$Design), col=unique(pData.train.13$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.13.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.13$Design_Color,pData.test.13$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.13), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.13 <- predict(svmOpt.GA.13,matrix.test.13.bapred.GA, type = "prob")
Prediction_GA.13$Prediction_GA.13 <- ifelse(Prediction_GA.13$transforming>0.50,"transforming","untransforming")
Prediction_GA.13 <- cbind(pData.test.13[,c(1:3)],TrueLabel=pData.test.13$Class,Prediction_GA.13)
write.table(Prediction_GA.13, file = paste("TestSet013_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.13 <- predict(svmFull.13, matrix.test.13.bapred.full, type = "prob")
Prediction_SVM_full.13$Prediction_SVM_full <- ifelse(Prediction_SVM_full.13$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.13 <- cbind(pData.test.13[,c(1:3)],TrueLabel=pData.test.13$Class,Prediction_SVM_full.13)
write.table(Prediction_SVM_full.13, file = paste("TestSet013_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 013 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet013_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.13$Prediction_GA.13), as.factor(Prediction_GA.13$TrueLabel))
sink()

sink("TestSet013_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.13$Prediction_SVM_full), as.factor(Prediction_SVM_full.13$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 013 ########################################################

Prediction_GA.13$Class <- as.factor(ifelse(Prediction_GA.13$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.13 <- roc(Prediction_GA.13$Class,                    # response vector (factor or character)
                 Prediction_GA.13$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.13$Class <- as.factor(ifelse(Prediction_SVM_full.13$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.13<- roc(Prediction_SVM_full.13$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.13$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 013 ####################################
Prediction_GA.13$Class_Code <- ifelse(Prediction_GA.13$TrueLabel == "transforming",1,0)
pr.GA.13 <- pr.curve(scores.class0 = Prediction_GA.13$transforming , weights.class0 = Prediction_GA.13$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.13, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.13$Class_Code <- ifelse(Prediction_SVM_full.13$TrueLabel == "transforming",1,0)
pr.full.13 <- pr.curve(scores.class0 = Prediction_SVM_full.13$transforming , weights.class0 = Prediction_SVM_full.13$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.13, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 014 = Batch 014 = IVIM #180131 #############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 014 and TestSet 014 ##########################################################
##############################################################################################################
pData.test.14  <- pData[pData$Batch==14,]    # set aside batch 14 / IVIM #180131
pData.train.14 <- pData[pData$Batch!=14,]    # use all remaining assays as training set
RAW_train.14   <- SAGA_RAW[,row.names(pData.train.14)]
RAW_test.14    <- SAGA_RAW[,row.names(pData.test.14)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.14$E), col=pData.train.14$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.14 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.14    <- qunormtrain(t(RAW_train.14$E))
matrix.train.14.qn <- log2(t(qunorm.train.14$xnorm))
matrix.train.14.qn <- avereps(matrix.train.14.qn, ID= RAW_train.14$genes$ProbeName)  
colnames(matrix.train.14.qn) <- row.names(pData.train.14)
boxplot(matrix.train.14.qn, col=pData.train.14$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.14.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.14$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.14         <- as.factor(ifelse(pData.train.14$Batch>14,pData.train.14$Batch-1,pData.train.14$Batch))                       
combat.train.14        <- combatba(t(matrix.train.14.qn), batch = batch.train.14)
matrix.train.14.batch  <- t(combat.train.14$xadj)
colnames(matrix.train.14.batch) <- row.names(pData.train.14)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.14.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.14$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 014 ######################################################
##############################################################################################################
fselect.14  <- genefilter(matrix.train.14.batch, filterfun(f1))
summary(fselect.14)
matrix.train.14 <-matrix.train.14.batch[fselect.14,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 014 #######################################################################
##############################################################################################################
matrix.train.14 <- (t(matrix.train.14))
labels.train.14 <- as.factor(pData.train.14$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.14 <- createMultiFolds(labels.train.14, k=10, times = 20)  

fullCtrl.14 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.14,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.14 <- train(matrix.train.14,labels.train.14,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.14)

svmFull.14  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 014 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.14      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.14,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.14$functions         <- caretFuncs
outerctrl.14$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 014 ##############################################################################
##############################################################################################################

system.time(rfe.14  <- rfe(matrix.train.14, labels.train.14, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.14,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.14   # 25 variables found by SVM-rfe
write.table(rfe.14$results, file = "TrainingSet014_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.14, type = c("g", "o"))
plot(rfe.14, type = c("g", "o"), xlim = c(0,61))

optFeatures.14 <- cbind(rfe.14$optVariables, Annotation[rfe.14$optVariables,])
write.table(optFeatures.14, file = "TrainingSet014_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 014 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20014) since performances were measured on identically resampled sets

rfeResamples.14 <- resamples(list("SVM_full.14" = svmFull.14,"SVM_RFE.14" = rfe.14))
sink("TrainingSet014_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.14)
sink()

modelDifferences.14 <- diff(rfeResamples.14)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet014_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.14)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 014 #########################
################################################################################################

matrix.train.rfe.14 <- matrix.train.14[,rfe.14$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(213)
index.GA.14 <- createMultiFolds(labels.train.14, k=10, times = 5)  

outerctrl.GA.14 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.14,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.14<- gafs(matrix.train.rfe.14, labels.train.14, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.14,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.14   # yields 10 features
optVars.GA.14 <- Annotation[GA.14$optVariables,]  # export optimal variables
write.table(optVars.GA.14, file = "TrainingSet014_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.14 <- GA.14$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.14 <- arrange(performance.external.14, Iter)

performance.14 <- GA.14$ga$internal               # export internal accuracy (within the resample)
performance.14$AccuracyExternal <- aggregate(performance.external.14$Accuracy, by=list(Iter=performance.external.14$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.14 <- data.frame(Iter = c(performance.14$Iter,performance.14$Iter), Accuracy = c(performance.14$AccuracyExternal,performance.14$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.14 <- GA.14$averages[GA.14$optIter,2]
accuracy.external.14

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.14 <- subset(performance.external.14,performance.external.14$Iter == GA.14$optIter)
accuracy.external.opt.14 <- accuracy.external.opt.14$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.14, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 014 ########################
##############################################################################################################
matrix.train.14.bapred.GA <- t(matrix.train.14.batch[GA.14$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.14  <- train(matrix.train.14.bapred.GA,labels.train.14,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.14)
svmOpt.GA.14  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.14.qn  <- log2(t(qunormaddon(qunorm.train.14, t(RAW_test.14$E))))   # use qunormtrain object to normalize test set
matrix.test.14.qn  <- avereps(matrix.test.14.qn, ID= RAW_test.14$genes$ProbeName)
boxplot(cbind(matrix.train.14.qn,matrix.test.14.qn),col=c(pData.train.14$IVIM_Color,pData.test.14$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.14.bapred  <- t(combatbaaddon(combat.train.14, t(matrix.test.14.qn), batch = as.factor(pData.test.14$Batch)))

plot(Rtsne(t(cbind(matrix.train.14.bapred,matrix.test.14.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.14$Design_Color,pData.test.14$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.14.bapred.GA   <- t(matrix.test.14.bapred[GA.14$optVariables,])
matrix.test.14.bapred.full <- t(matrix.test.14.bapred[colnames(matrix.train.14),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.14.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.14$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.14$Design), col=unique(pData.train.14$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.14.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.14$Design_Color,pData.test.14$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.14), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.14 <- predict(svmOpt.GA.14,matrix.test.14.bapred.GA, type = "prob")
Prediction_GA.14$Prediction_GA.14 <- ifelse(Prediction_GA.14$transforming>0.50,"transforming","untransforming")
Prediction_GA.14 <- cbind(pData.test.14[,c(1:3)],TrueLabel=pData.test.14$Class,Prediction_GA.14)
write.table(Prediction_GA.14, file = paste("TestSet014_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.14 <- predict(svmFull.14, matrix.test.14.bapred.full, type = "prob")
Prediction_SVM_full.14$Prediction_SVM_full <- ifelse(Prediction_SVM_full.14$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.14 <- cbind(pData.test.14[,c(1:3)],TrueLabel=pData.test.14$Class,Prediction_SVM_full.14)
write.table(Prediction_SVM_full.14, file = paste("TestSet014_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 014 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet014_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.14$Prediction_GA.14), as.factor(Prediction_GA.14$TrueLabel))
sink()

sink("TestSet014_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.14$Prediction_SVM_full), as.factor(Prediction_SVM_full.14$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 014 ########################################################

Prediction_GA.14$Class <- as.factor(ifelse(Prediction_GA.14$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.14 <- roc(Prediction_GA.14$Class,                    # response vector (factor or character)
                 Prediction_GA.14$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.14$Class <- as.factor(ifelse(Prediction_SVM_full.14$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.14<- roc(Prediction_SVM_full.14$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.14$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 014 ####################################
Prediction_GA.14$Class_Code <- ifelse(Prediction_GA.14$TrueLabel == "transforming",1,0)
pr.GA.14 <- pr.curve(scores.class0 = Prediction_GA.14$transforming , weights.class0 = Prediction_GA.14$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.14, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.14$Class_Code <- ifelse(Prediction_SVM_full.14$TrueLabel == "transforming",1,0)
pr.full.14 <- pr.curve(scores.class0 = Prediction_SVM_full.14$transforming , weights.class0 = Prediction_SVM_full.14$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.14, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 015 = Batch 015 = IVIM #180620 ############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 015 and TestSet 015 ############################################################
##############################################################################################################
pData.test.15  <- pData[pData$Batch==15,]    # set aside batch 15 / IVIM #180620
pData.train.15 <- pData[pData$Batch!=15,]    # use all remaining assays as training set
RAW_train.15   <- SAGA_RAW[,row.names(pData.train.15)]
RAW_test.15    <- SAGA_RAW[,row.names(pData.test.15)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.15$E), col=pData.train.15$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.15 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.15    <- qunormtrain(t(RAW_train.15$E))
matrix.train.15.qn <- log2(t(qunorm.train.15$xnorm))
matrix.train.15.qn <- avereps(matrix.train.15.qn, ID= RAW_train.15$genes$ProbeName)  
colnames(matrix.train.15.qn) <- row.names(pData.train.15)
boxplot(matrix.train.15.qn, col=pData.train.15$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.15.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.15$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.15         <- as.factor(ifelse(pData.train.15$Batch>15,pData.train.15$Batch-1,pData.train.15$Batch))                       
combat.train.15        <- combatba(t(matrix.train.15.qn), batch = batch.train.15)
matrix.train.15.batch  <- t(combat.train.15$xadj)
colnames(matrix.train.15.batch) <- row.names(pData.train.15)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.15.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.15$Design_Color, pch=16, cex=1.3) 
##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 015 ######################################################
##############################################################################################################
fselect.15  <- genefilter(matrix.train.15.batch, filterfun(f1))
summary(fselect.15)
matrix.train.15 <-matrix.train.15.batch[fselect.15,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 015 #######################################################################
##############################################################################################################
matrix.train.15 <- (t(matrix.train.15))
labels.train.15 <- as.factor(pData.train.15$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.15 <- createMultiFolds(labels.train.15, k=10, times = 20)  

fullCtrl.15 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.15,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.15 <- train(matrix.train.15,labels.train.15,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.15)

svmFull.15  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 015 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.15      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.15,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.15$functions         <- caretFuncs
outerctrl.15$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 015 ##############################################################################
##############################################################################################################

system.time(rfe.15  <- rfe(matrix.train.15, labels.train.15, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.15,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.15   # 21 variables found by SVM-rfe
write.table(rfe.15$results, file = "TrainingSet015_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.15, type = c("g", "o"))
plot(rfe.15, type = c("g", "o"), xlim = c(0,61))

optFeatures.15 <- cbind(rfe.15$optVariables, Annotation[rfe.15$optVariables,])
write.table(optFeatures.15, file = "TrainingSet015_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 015 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20015) since performances were measured on identically resampled sets

rfeResamples.15 <- resamples(list("SVM_full.15" = svmFull.15,"SVM_RFE.15" = rfe.15))
sink("TrainingSet015_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.15)
sink()

modelDifferences.15 <- diff(rfeResamples.15)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet015_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.15)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 015 #########################
################################################################################################

matrix.train.rfe.15 <- matrix.train.15[,rfe.15$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(123)
index.GA.15 <- createMultiFolds(labels.train.15, k=10, times = 5)  

outerctrl.GA.15 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.15,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.15<- gafs(matrix.train.rfe.15, labels.train.15, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.15,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.15   # yields 10 features
optVars.GA.15 <- Annotation[GA.15$optVariables,]  # export optimal variables
write.table(optVars.GA.15, file = "TrainingSet015_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.15 <- GA.15$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.15 <- arrange(performance.external.15, Iter)

performance.15 <- GA.15$ga$internal               # export internal accuracy (within the resample)
performance.15$AccuracyExternal <- aggregate(performance.external.15$Accuracy, by=list(Iter=performance.external.15$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.15 <- data.frame(Iter = c(performance.15$Iter,performance.15$Iter), Accuracy = c(performance.15$AccuracyExternal,performance.15$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.15 <- GA.15$averages[GA.15$optIter,2]
accuracy.external.15

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.15 <- subset(performance.external.15,performance.external.15$Iter == GA.15$optIter)
accuracy.external.opt.15 <- accuracy.external.opt.15$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.15, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 015 ########################
##############################################################################################################
matrix.train.15.bapred.GA <- t(matrix.train.15.batch[GA.15$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.15  <- train(matrix.train.15.bapred.GA,labels.train.15,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.15)
svmOpt.GA.15  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.15.qn  <- log2(t(qunormaddon(qunorm.train.15, t(RAW_test.15$E))))   # use qunormtrain object to normalize test set
matrix.test.15.qn  <- avereps(matrix.test.15.qn, ID= RAW_test.15$genes$ProbeName)
boxplot(cbind(matrix.train.15.qn,matrix.test.15.qn),col=c(pData.train.15$IVIM_Color,pData.test.15$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.15.bapred  <- t(combatbaaddon(combat.train.15, t(matrix.test.15.qn), batch = as.factor(pData.test.15$Batch)))

plot(Rtsne(t(cbind(matrix.train.15.bapred,matrix.test.15.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.15$Design_Color,pData.test.15$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.15.bapred.GA   <- t(matrix.test.15.bapred[GA.15$optVariables,])
matrix.test.15.bapred.full <- t(matrix.test.15.bapred[colnames(matrix.train.15),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.15.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.15$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.15$Design), col=unique(pData.train.15$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.15.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.15$Design_Color,pData.test.15$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.15), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.15 <- predict(svmOpt.GA.15,matrix.test.15.bapred.GA, type = "prob")
Prediction_GA.15$Prediction_GA.15 <- ifelse(Prediction_GA.15$transforming>0.50,"transforming","untransforming")
Prediction_GA.15 <- cbind(pData.test.15[,c(1:3)],TrueLabel=pData.test.15$Class,Prediction_GA.15)
write.table(Prediction_GA.15, file = paste("TestSet015_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.15 <- predict(svmFull.15, matrix.test.15.bapred.full, type = "prob")
Prediction_SVM_full.15$Prediction_SVM_full <- ifelse(Prediction_SVM_full.15$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.15 <- cbind(pData.test.15[,c(1:3)],TrueLabel=pData.test.15$Class,Prediction_SVM_full.15)
write.table(Prediction_SVM_full.15, file = paste("TestSet015_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 015 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet015_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.15$Prediction_GA.15), as.factor(Prediction_GA.15$TrueLabel))
sink()

sink("TestSet015_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.15$Prediction_SVM_full), as.factor(Prediction_SVM_full.15$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 015 ########################################################

Prediction_GA.15$Class <- as.factor(ifelse(Prediction_GA.15$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.15 <- roc(Prediction_GA.15$Class,                    # response vector (factor or character)
                 Prediction_GA.15$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.15$Class <- as.factor(ifelse(Prediction_SVM_full.15$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.15<- roc(Prediction_SVM_full.15$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.15$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 015 ####################################
Prediction_GA.15$Class_Code <- ifelse(Prediction_GA.15$TrueLabel == "transforming",1,0)
pr.GA.15 <- pr.curve(scores.class0 = Prediction_GA.15$transforming , weights.class0 = Prediction_GA.15$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.15, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.15$Class_Code <- ifelse(Prediction_SVM_full.15$TrueLabel == "transforming",1,0)
pr.full.15 <- pr.curve(scores.class0 = Prediction_SVM_full.15$transforming , weights.class0 = Prediction_SVM_full.15$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.15, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 016 = Batch 016 = IVIM #180523A ############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 016 and TestSet 016 ##########################################################
##############################################################################################################
pData.test.16  <- pData[pData$Batch==16,]    # set aside batch 16 / IVIM #180523A
pData.train.16 <- pData[pData$Batch!=16,]    # use all remaining assays as training set
RAW_train.16   <- SAGA_RAW[,row.names(pData.train.16)]
RAW_test.16    <- SAGA_RAW[,row.names(pData.test.16)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.16$E), col=pData.train.16$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.16 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.16    <- qunormtrain(t(RAW_train.16$E))
matrix.train.16.qn <- log2(t(qunorm.train.16$xnorm))
matrix.train.16.qn <- avereps(matrix.train.16.qn, ID= RAW_train.16$genes$ProbeName)  
colnames(matrix.train.16.qn) <- row.names(pData.train.16)
boxplot(matrix.train.16.qn, col=pData.train.16$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.16.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.16$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.16         <- as.factor(ifelse(pData.train.16$Batch>16,pData.train.16$Batch-1,pData.train.16$Batch))                       
combat.train.16        <- combatba(t(matrix.train.16.qn), batch = batch.train.16)
matrix.train.16.batch  <- t(combat.train.16$xadj)
colnames(matrix.train.16.batch) <- row.names(pData.train.16)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.16.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.16$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 016 ######################################################
##############################################################################################################
fselect.16  <- genefilter(matrix.train.16.batch, filterfun(f1))
summary(fselect.16)
matrix.train.16 <-matrix.train.16.batch[fselect.16,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 016 #######################################################################
##############################################################################################################
matrix.train.16 <- (t(matrix.train.16))
labels.train.16 <- as.factor(pData.train.16$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.16 <- createMultiFolds(labels.train.16, k=10, times = 20)  

fullCtrl.16 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.16,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.16 <- train(matrix.train.16,labels.train.16,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.16)

svmFull.16  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 016 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.16      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.16,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.16$functions         <- caretFuncs
outerctrl.16$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 016 ##############################################################################
##############################################################################################################

system.time(rfe.16  <- rfe(matrix.train.16, labels.train.16, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.16,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.16   # 20 variables found by SVM-rfe
write.table(rfe.16$results, file = "TrainingSet016_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.16, type = c("g", "o"))
plot(rfe.16, type = c("g", "o"), xlim = c(0,61))

optFeatures.16 <- cbind(rfe.16$optVariables, Annotation[rfe.16$optVariables,])
write.table(optFeatures.16, file = "TrainingSet016_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 016 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20016) since performances were measured on identically resampled sets

rfeResamples.16 <- resamples(list("SVM_full.16" = svmFull.16,"SVM_RFE.16" = rfe.16))
sink("TrainingSet016_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.16)
sink()

modelDifferences.16 <- diff(rfeResamples.16)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet016_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.16)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 016 #########################
################################################################################################

matrix.train.rfe.16 <- matrix.train.16[,rfe.16$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.16 <- createMultiFolds(labels.train.16, k=10, times = 5)  

outerctrl.GA.16 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.16,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.16<- gafs(matrix.train.rfe.16, labels.train.16, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.16,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.16   # yields 10 features
optVars.GA.16 <- Annotation[GA.16$optVariables,]  # export optimal variables
write.table(optVars.GA.16, file = "TrainingSet016_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.16 <- GA.16$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.16 <- arrange(performance.external.16, Iter)

performance.16 <- GA.16$ga$internal               # export internal accuracy (within the resample)
performance.16$AccuracyExternal <- aggregate(performance.external.16$Accuracy, by=list(Iter=performance.external.16$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.16 <- data.frame(Iter = c(performance.16$Iter,performance.16$Iter), Accuracy = c(performance.16$AccuracyExternal,performance.16$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.16 <- GA.16$averages[GA.16$optIter,2]
accuracy.external.16

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.16 <- subset(performance.external.16,performance.external.16$Iter == GA.16$optIter)
accuracy.external.opt.16 <- accuracy.external.opt.16$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.16, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 016 ########################
##############################################################################################################
matrix.train.16.bapred.GA <- t(matrix.train.16.batch[GA.16$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.16  <- train(matrix.train.16.bapred.GA,labels.train.16,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.16)
svmOpt.GA.16  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.16.qn  <- log2(t(qunormaddon(qunorm.train.16, t(RAW_test.16$E))))   # use qunormtrain object to normalize test set
matrix.test.16.qn  <- avereps(matrix.test.16.qn, ID= RAW_test.16$genes$ProbeName)
boxplot(cbind(matrix.train.16.qn,matrix.test.16.qn),col=c(pData.train.16$IVIM_Color,pData.test.16$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.16.bapred  <- t(combatbaaddon(combat.train.16, t(matrix.test.16.qn), batch = as.factor(pData.test.16$Batch)))

plot(Rtsne(t(cbind(matrix.train.16.bapred,matrix.test.16.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.16$Design_Color,pData.test.16$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.16.bapred.GA   <- t(matrix.test.16.bapred[GA.16$optVariables,])
matrix.test.16.bapred.full <- t(matrix.test.16.bapred[colnames(matrix.train.16),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.16.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.16$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.16$Design), col=unique(pData.train.16$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.16.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.16$Design_Color,pData.test.16$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.16), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.16 <- predict(svmOpt.GA.16,matrix.test.16.bapred.GA, type = "prob")
Prediction_GA.16$Prediction_GA.16 <- ifelse(Prediction_GA.16$transforming>0.50,"transforming","untransforming")
Prediction_GA.16 <- cbind(pData.test.16[,c(1:3)],TrueLabel=pData.test.16$Class,Prediction_GA.16)
write.table(Prediction_GA.16, file = paste("TestSet016_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.16 <- predict(svmFull.16, matrix.test.16.bapred.full, type = "prob")
Prediction_SVM_full.16$Prediction_SVM_full <- ifelse(Prediction_SVM_full.16$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.16 <- cbind(pData.test.16[,c(1:3)],TrueLabel=pData.test.16$Class,Prediction_SVM_full.16)
write.table(Prediction_SVM_full.16, file = paste("TestSet016_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 016 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet016_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.16$Prediction_GA.16), as.factor(Prediction_GA.16$TrueLabel))
sink()

sink("TestSet016_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.16$Prediction_SVM_full), as.factor(Prediction_SVM_full.16$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 016 ########################################################

Prediction_GA.16$Class <- as.factor(ifelse(Prediction_GA.16$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.16 <- roc(Prediction_GA.16$Class,                    # response vector (factor or character)
                 Prediction_GA.16$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.16$Class <- as.factor(ifelse(Prediction_SVM_full.16$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.16<- roc(Prediction_SVM_full.16$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.16$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 016 ####################################
Prediction_GA.16$Class_Code <- ifelse(Prediction_GA.16$TrueLabel == "transforming",1,0)
pr.GA.16 <- pr.curve(scores.class0 = Prediction_GA.16$transforming , weights.class0 = Prediction_GA.16$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.16, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.16$Class_Code <- ifelse(Prediction_SVM_full.16$TrueLabel == "transforming",1,0)
pr.full.16 <- pr.curve(scores.class0 = Prediction_SVM_full.16$transforming , weights.class0 = Prediction_SVM_full.16$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.16, rand.plot = TRUE, legend = F, color = 1,main = "")


#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 017 = Batch 017 = IVIM #180523B ############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 017 and TestSet 017 ############################################################
##############################################################################################################
pData.test.17  <- pData[pData$Batch==17,]    # set aside batch 17 / IVIM #180523B
pData.train.17 <- pData[pData$Batch!=17,]    # use all remaining assays as training set
RAW_train.17   <- SAGA_RAW[,row.names(pData.train.17)]
RAW_test.17    <- SAGA_RAW[,row.names(pData.test.17)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.17$E), col=pData.train.17$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.17 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.17    <- qunormtrain(t(RAW_train.17$E))
matrix.train.17.qn <- log2(t(qunorm.train.17$xnorm))
matrix.train.17.qn <- avereps(matrix.train.17.qn, ID= RAW_train.17$genes$ProbeName)  
colnames(matrix.train.17.qn) <- row.names(pData.train.17)
boxplot(matrix.train.17.qn, col=pData.train.17$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.17.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.17$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.17         <- as.factor(ifelse(pData.train.17$Batch>17,pData.train.17$Batch-1,pData.train.17$Batch))                       
combat.train.17        <- combatba(t(matrix.train.17.qn), batch = batch.train.17)
matrix.train.17.batch  <- t(combat.train.17$xadj)
colnames(matrix.train.17.batch) <- row.names(pData.train.17)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.17.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.17$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 017 ######################################################
##############################################################################################################
fselect.17  <- genefilter(matrix.train.17.batch, filterfun(f1))
summary(fselect.17)
matrix.train.17 <-matrix.train.17.batch[fselect.17,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 017 #######################################################################
##############################################################################################################
matrix.train.17 <- (t(matrix.train.17))
labels.train.17 <- as.factor(pData.train.17$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(2101)
index.17 <- createMultiFolds(labels.train.17, k=10, times = 20)  

fullCtrl.17 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.17,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(123)
svmFull.17 <- train(matrix.train.17,labels.train.17,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.17)

svmFull.17  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 017 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.17      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.17,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.17$functions         <- caretFuncs
outerctrl.17$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 017 ##############################################################################
##############################################################################################################

system.time(rfe.17  <- rfe(matrix.train.17, labels.train.17, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.17,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.17   # 1336 predictors found by SVM-rfe ==> rfe failed to reduce predictors 
write.table(rfe.17$results, file = "TrainingSet017_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.17, type = c("g", "o"))
plot(rfe.17, type = c("g", "o"), xlim = c(0,61))  # second best solution = 22 predictors 

# manually select the second best solution (number of predictors=22)
# 1336 predictors: CV-accuracy = 0.9078
# 22 predictors :  CV-accuracy = 0.9008

varImp_17 <- varImp(svmFull.17)$importance                              # extract the variable importance (AUCROC) of all 1336 predictors 
varImp_17 <- varImp_17[order(varImp_17$transforming,decreasing = T),]   # rank by AUCROC
varImp_17 <- cbind(varImp_17,Annotation[row.names(varImp_17),])        
optFeatures.17 <- varImp_17[c(1:22),]                                   # select the 22 most important predictors
write.table(optFeatures.17, file = "TrainingSet017_optVars_rfe.txt", sep="\t",col.names=NA)



################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 017 #########################
################################################################################################

matrix.train.rfe.17 <- matrix.train.17[,row.names(optFeatures.17)]   

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.17 <- createMultiFolds(labels.train.17, k=10, times = 5)  

outerctrl.GA.17 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.17,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.17<- gafs(matrix.train.rfe.17, labels.train.17, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.17,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.17   # yields 12 features
optVars.GA.17 <- Annotation[GA.17$optVariables,]  # export optimal variables
write.table(optVars.GA.17, file = "TrainingSet017_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.17 <- GA.17$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.17 <- arrange(performance.external.17, Iter)

performance.17 <- GA.17$ga$internal               # export internal accuracy (within the resample)
performance.17$AccuracyExternal <- aggregate(performance.external.17$Accuracy, by=list(Iter=performance.external.17$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.17 <- data.frame(Iter = c(performance.17$Iter,performance.17$Iter), Accuracy = c(performance.17$AccuracyExternal,performance.17$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.17 <- GA.17$averages[GA.17$optIter,2]
accuracy.external.17

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.17 <- subset(performance.external.17,performance.external.17$Iter == GA.17$optIter)
accuracy.external.opt.17 <- accuracy.external.opt.17$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.17, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()


##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 017 ########################
##############################################################################################################
matrix.train.17.bapred.GA <- t(matrix.train.17.batch[GA.17$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.17  <- train(matrix.train.17.bapred.GA,labels.train.17,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.17)
svmOpt.GA.17  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.17.qn  <- log2(t(qunormaddon(qunorm.train.17, t(RAW_test.17$E))))   # use qunormtrain object to normalize test set
matrix.test.17.qn  <- avereps(matrix.test.17.qn, ID= RAW_test.17$genes$ProbeName)
boxplot(cbind(matrix.train.17.qn,matrix.test.17.qn),col=c(pData.train.17$IVIM_Color,pData.test.17$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.17.bapred  <- t(combatbaaddon(combat.train.17, t(matrix.test.17.qn), batch = as.factor(pData.test.17$Batch)))

plot(Rtsne(t(cbind(matrix.train.17.bapred,matrix.test.17.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.17$Design_Color,pData.test.17$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.17.bapred.GA   <- t(matrix.test.17.bapred[GA.17$optVariables,])
matrix.test.17.bapred.full <- t(matrix.test.17.bapred[colnames(matrix.train.17),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.17.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.17$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.17$Design), col=unique(pData.train.17$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.17.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.17$Design_Color,pData.test.17$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.17), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.17 <- predict(svmOpt.GA.17,matrix.test.17.bapred.GA, type = "prob")
Prediction_GA.17$Prediction_GA.17 <- ifelse(Prediction_GA.17$transforming>0.50,"transforming","untransforming")
Prediction_GA.17 <- cbind(pData.test.17[,c(1:3)],TrueLabel=pData.test.17$Class,Prediction_GA.17)
write.table(Prediction_GA.17, file = paste("TestSet017_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.17 <- predict(svmFull.17, matrix.test.17.bapred.full, type = "prob")
Prediction_SVM_full.17$Prediction_SVM_full <- ifelse(Prediction_SVM_full.17$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.17 <- cbind(pData.test.17[,c(1:3)],TrueLabel=pData.test.17$Class,Prediction_SVM_full.17)
write.table(Prediction_SVM_full.17, file = paste("TestSet017_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 017 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet017_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.17$Prediction_GA.17), as.factor(Prediction_GA.17$TrueLabel))
sink()

sink("TestSet017_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.17$Prediction_SVM_full), as.factor(Prediction_SVM_full.17$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 017 ########################################################

Prediction_GA.17$Class <- as.factor(ifelse(Prediction_GA.17$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.17 <- roc(Prediction_GA.17$Class,                    # response vector (factor or character)
                 Prediction_GA.17$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.17$Class <- as.factor(ifelse(Prediction_SVM_full.17$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.17<- roc(Prediction_SVM_full.17$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.17$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 017 ####################################
Prediction_GA.17$Class_Code <- ifelse(Prediction_GA.17$TrueLabel == "transforming",1,0)
pr.GA.17 <- pr.curve(scores.class0 = Prediction_GA.17$transforming , weights.class0 = Prediction_GA.17$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.17, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.17$Class_Code <- ifelse(Prediction_SVM_full.17$TrueLabel == "transforming",1,0)
pr.full.17 <- pr.curve(scores.class0 = Prediction_SVM_full.17$transforming , weights.class0 = Prediction_SVM_full.17$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.17, rand.plot = TRUE, legend = F, color = 1,main = "")


#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 018 = Batch 018 = IVIM #180801 ###############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 018 and TestSet 018 ############################################################
##############################################################################################################
pData.test.18  <- pData[pData$Batch==18,]    # set aside batch 18 / IVIM #180801
pData.train.18 <- pData[pData$Batch!=18,]    # use all remaining assays as training set
RAW_train.18   <- SAGA_RAW[,row.names(pData.train.18)]
RAW_test.18    <- SAGA_RAW[,row.names(pData.test.18)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.18$E), col=pData.train.18$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.18 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.18    <- qunormtrain(t(RAW_train.18$E))
matrix.train.18.qn <- log2(t(qunorm.train.18$xnorm))
matrix.train.18.qn <- avereps(matrix.train.18.qn, ID= RAW_train.18$genes$ProbeName)  
colnames(matrix.train.18.qn) <- row.names(pData.train.18)
boxplot(matrix.train.18.qn, col=pData.train.18$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.18.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.18$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.18         <- as.factor(ifelse(pData.train.18$Batch>18,pData.train.18$Batch-1,pData.train.18$Batch))                       
combat.train.18        <- combatba(t(matrix.train.18.qn), batch = batch.train.18)
matrix.train.18.batch  <- t(combat.train.18$xadj)
colnames(matrix.train.18.batch) <- row.names(pData.train.18)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.18.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.18$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 018 ######################################################
##############################################################################################################
fselect.18  <- genefilter(matrix.train.18.batch, filterfun(f1))
summary(fselect.18)
matrix.train.18 <-matrix.train.18.batch[fselect.18,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 018 #######################################################################
##############################################################################################################
matrix.train.18 <- (t(matrix.train.18))
labels.train.18 <- as.factor(pData.train.18$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.18 <- createMultiFolds(labels.train.18, k=10, times = 20)  

fullCtrl.18 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.18,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.18 <- train(matrix.train.18,labels.train.18,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.18)

svmFull.18  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 018 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.18      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.18,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.18$functions         <- caretFuncs
outerctrl.18$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 018 ##############################################################################
##############################################################################################################

system.time(rfe.18  <- rfe(matrix.train.18, labels.train.18, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.18,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.18   # 19 variables found by SVM-rfe
write.table(rfe.18$results, file = "TrainingSet018_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.18, type = c("g", "o"))
plot(rfe.18, type = c("g", "o"), xlim = c(0,61))

optFeatures.18 <- cbind(rfe.18$optVariables, Annotation[rfe.18$optVariables,])
write.table(optFeatures.18, file = "TrainingSet018_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 018 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20018) since performances were measured on identically resampled sets

rfeResamples.18 <- resamples(list("SVM_full.18" = svmFull.18,"SVM_RFE.18" = rfe.18))
sink("TrainingSet018_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.18)
sink()

modelDifferences.18 <- diff(rfeResamples.18)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet018_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.18)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 018 #########################
################################################################################################

matrix.train.rfe.18 <- matrix.train.18[,rfe.18$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.18 <- createMultiFolds(labels.train.18, k=10, times = 5)  

outerctrl.GA.18 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.18,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.18<- gafs(matrix.train.rfe.18, labels.train.18, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.18,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.18   # yields 10 features
optVars.GA.18 <- Annotation[GA.18$optVariables,]  # export optimal variables
write.table(optVars.GA.18, file = "TrainingSet018_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.18 <- GA.18$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.18 <- arrange(performance.external.18, Iter)

performance.18 <- GA.18$ga$internal               # export internal accuracy (within the resample)
performance.18$AccuracyExternal <- aggregate(performance.external.18$Accuracy, by=list(Iter=performance.external.18$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.18 <- data.frame(Iter = c(performance.18$Iter,performance.18$Iter), Accuracy = c(performance.18$AccuracyExternal,performance.18$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.18 <- GA.18$averages[GA.18$optIter,2]
accuracy.external.18

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.18 <- subset(performance.external.18,performance.external.18$Iter == GA.18$optIter)
accuracy.external.opt.18 <- accuracy.external.opt.18$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.18, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 018 ########################
##############################################################################################################
matrix.train.18.bapred.GA <- t(matrix.train.18.batch[GA.18$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.18  <- train(matrix.train.18.bapred.GA,labels.train.18,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.18)
svmOpt.GA.18  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.18.qn  <- log2(t(qunormaddon(qunorm.train.18, t(RAW_test.18$E))))   # use qunormtrain object to normalize test set
matrix.test.18.qn  <- avereps(matrix.test.18.qn, ID= RAW_test.18$genes$ProbeName)
boxplot(cbind(matrix.train.18.qn,matrix.test.18.qn),col=c(pData.train.18$IVIM_Color,pData.test.18$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.18.bapred  <- t(combatbaaddon(combat.train.18, t(matrix.test.18.qn), batch = as.factor(pData.test.18$Batch)))

plot(Rtsne(t(cbind(matrix.train.18.bapred,matrix.test.18.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.18$Design_Color,pData.test.18$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.18.bapred.GA   <- t(matrix.test.18.bapred[GA.18$optVariables,])
matrix.test.18.bapred.full <- t(matrix.test.18.bapred[colnames(matrix.train.18),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.18.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.18$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.18$Design), col=unique(pData.train.18$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.18.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.18$Design_Color,pData.test.18$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.18), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.18 <- predict(svmOpt.GA.18,matrix.test.18.bapred.GA, type = "prob")
Prediction_GA.18$Prediction_GA.18 <- ifelse(Prediction_GA.18$transforming>0.50,"transforming","untransforming")
Prediction_GA.18 <- cbind(pData.test.18[,c(1:3)],TrueLabel=pData.test.18$Class,Prediction_GA.18)
write.table(Prediction_GA.18, file = paste("TestSet018_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.18 <- predict(svmFull.18, matrix.test.18.bapred.full, type = "prob")
Prediction_SVM_full.18$Prediction_SVM_full <- ifelse(Prediction_SVM_full.18$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.18 <- cbind(pData.test.18[,c(1:3)],TrueLabel=pData.test.18$Class,Prediction_SVM_full.18)
write.table(Prediction_SVM_full.18, file = paste("TestSet018_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 018 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet018_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.18$Prediction_GA.18), as.factor(Prediction_GA.18$TrueLabel))
sink()

sink("TestSet018_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.18$Prediction_SVM_full), as.factor(Prediction_SVM_full.18$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 018 ########################################################

Prediction_GA.18$Class <- as.factor(ifelse(Prediction_GA.18$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.18 <- roc(Prediction_GA.18$Class,                    # response vector (factor or character)
                 Prediction_GA.18$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.18$Class <- as.factor(ifelse(Prediction_SVM_full.18$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.18<- roc(Prediction_SVM_full.18$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.18$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 018 ####################################
Prediction_GA.18$Class_Code <- ifelse(Prediction_GA.18$TrueLabel == "transforming",1,0)
pr.GA.18 <- pr.curve(scores.class0 = Prediction_GA.18$transforming , weights.class0 = Prediction_GA.18$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.18, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.18$Class_Code <- ifelse(Prediction_SVM_full.18$TrueLabel == "transforming",1,0)
pr.full.18 <- pr.curve(scores.class0 = Prediction_SVM_full.18$transforming , weights.class0 = Prediction_SVM_full.18$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.18, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 019 = Batch 019 = IVIM #180822 #############################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 019 and TestSet 019 ############################################################
##############################################################################################################
pData.test.19  <- pData[pData$Batch==19,]    # set aside batch 19 / IVIM #180822
pData.train.19 <- pData[pData$Batch!=19,]    # use all remaining assays as training set
RAW_train.19   <- SAGA_RAW[,row.names(pData.train.19)]
RAW_test.19    <- SAGA_RAW[,row.names(pData.test.19)]

##############################################################################################################
#### 2.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(RAW_train.19$E), col=pData.train.19$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.19 quantile normalization using addon-quantile normalization ########################################### 
qunorm.train.19    <- qunormtrain(t(RAW_train.19$E))
matrix.train.19.qn <- log2(t(qunorm.train.19$xnorm))
matrix.train.19.qn <- avereps(matrix.train.19.qn, ID= RAW_train.19$genes$ProbeName)  
colnames(matrix.train.19.qn) <- row.names(pData.train.19)
boxplot(matrix.train.19.qn, col=pData.train.19$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.19.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.19$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.19         <- as.factor(ifelse(pData.train.19$Batch>19,pData.train.19$Batch-1,pData.train.19$Batch))                       
combat.train.19        <- combatba(t(matrix.train.19.qn), batch = batch.train.19)
matrix.train.19.batch  <- t(combat.train.19$xadj)
colnames(matrix.train.19.batch) <- row.names(pData.train.19)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.19.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.19$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 019 ######################################################
##############################################################################################################
fselect.19  <- genefilter(matrix.train.19.batch, filterfun(f1))
summary(fselect.19)
matrix.train.19 <-matrix.train.19.batch[fselect.19,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 019 #######################################################################
##############################################################################################################
matrix.train.19 <- (t(matrix.train.19))
labels.train.19 <- as.factor(pData.train.19$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.19 <- createMultiFolds(labels.train.19, k=10, times = 20)  

fullCtrl.19 <- trainControl(method = "repeatedcv",repeats = 20,
                            index = index.19,
                            summaryFunction = fiveStats,
                            classProbs = TRUE,
                            allowParallel = TRUE)

set.seed(721)
svmFull.19 <- train(matrix.train.19,labels.train.19,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.19)

svmFull.19  

##############################################################################################################
#### 3. SVM-RFE TrainingSet 019 ###############################################################################
##############################################################################################################

#### 3.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.19      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                index = index.19,
                                saveDetails = TRUE,
                                returnResamp="final", 
                                verbose = TRUE, 
                                seeds = seeds.rfe,
                                allowParallel = TRUE)

outerctrl.19$functions         <- caretFuncs
outerctrl.19$functions$summary <- fiveStats

#### 3.2 SVM-RFE TrainingSet 019 ##############################################################################
##############################################################################################################

system.time(rfe.19  <- rfe(matrix.train.19, labels.train.19, 
                           sizes=FeatureNumbers,
                           rfeControl=outerctrl.19,
                           metric = "Accuracy",
                           method="svmRadial",
                           tuneLength = 20,
                           trControl = innerctrl))

rfe.19   # 37 variables found by SVM-rfe
write.table(rfe.19$results, file = "TrainingSet019_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.19, type = c("g", "o"))
plot(rfe.19, type = c("g", "o"), xlim = c(0,61))

optFeatures.19 <- cbind(rfe.19$optVariables, Annotation[rfe.19$optVariables,])
write.table(optFeatures.19, file = "TrainingSet019_optVars_rfe.txt", sep="\t",col.names=NA)

#### 3.3 compare resampling performances TrainingSet 019 ######################################################
##############################################################################################################

### paired t-test according to Hothorn (2005) and Eugster (20019) since performances were measured on identically resampled sets

rfeResamples.19 <- resamples(list("SVM_full.19" = svmFull.19,"SVM_RFE.19" = rfe.19))
sink("TrainingSet019_resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.19)
sink()

modelDifferences.19 <- diff(rfeResamples.19)  # paired t-test for H0: difference = 0 between the different models. 
sink("TrainingSet019_Modeldifferences_rfe vs full.txt", append = TRUE)
summary(modelDifferences.19)
sink()


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 019 #########################
################################################################################################

matrix.train.rfe.19 <- matrix.train.19[,rfe.19$optVariables]   # subset fot the 20 optVars from rfe

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.19 <- createMultiFolds(labels.train.19, k=10, times = 5)  

outerctrl.GA.19 <- gafsControl(functions = svmGA,
                               method = "repeatedcv", repeats = 5,
                               index = index.GA.19,                                       
                               seeds = seeds.GA,                                      
                               returnResamp="all", 
                               verbose = TRUE,
                               maximize = c(internal = TRUE,
                                            external = TRUE),
                               allowParallel = TRUE)                                  


system.time(GA.19<- gafs(matrix.train.rfe.19, labels.train.19, 
                         iters = 40,
                         popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                         gafsControl = outerctrl.GA.19,
                         metric = "Accuracy",
                         method = "svmRadial",
                         # inner loop control for hyperparameter tuning
                         tuneLength = 12,
                         trControl = trainControl(method = "repeatedcv",
                                                  repeats = 2,
                                                  allowParallel = FALSE)))

### 4.2 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.19   # yields 15 features
optVars.GA.19 <- Annotation[GA.19$optVariables,]  # export optimal variables
write.table(optVars.GA.19, file = "TrainingSet019_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.19 <- GA.19$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.19 <- arrange(performance.external.19, Iter)

performance.19 <- GA.19$ga$internal               # export internal accuracy (within the resample)
performance.19$AccuracyExternal <- aggregate(performance.external.19$Accuracy, by=list(Iter=performance.external.19$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.19 <- data.frame(Iter = c(performance.19$Iter,performance.19$Iter), Accuracy = c(performance.19$AccuracyExternal,performance.19$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.19 <- GA.19$averages[GA.19$optIter,2]
accuracy.external.19

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.19 <- subset(performance.external.19,performance.external.19$Iter == GA.19$optIter)
accuracy.external.opt.19 <- accuracy.external.opt.19$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.19, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()

##############################################################################################################
##### 7 Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 019 ########################
##############################################################################################################
matrix.train.19.bapred.GA <- t(matrix.train.19.batch[GA.19$optVariables,]) # subset training matrix from step 2 to optimal predictors

#### 7.1 train SVM training set reduced to optVars  = SAGA classifier #######################################
#############################################################################################################
set.seed(721)
svmOpt.GA.19  <- train(matrix.train.19.bapred.GA,labels.train.19,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.19)
svmOpt.GA.19  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, training of classifier #############
##############################################################################################################

#### 7.2 Addon quantile normalization of test set ############################################################ 
##############################################################################################################
matrix.test.19.qn  <- log2(t(qunormaddon(qunorm.train.19, t(RAW_test.19$E))))   # use qunormtrain object to normalize test set
matrix.test.19.qn  <- avereps(matrix.test.19.qn, ID= RAW_test.19$genes$ProbeName)
boxplot(cbind(matrix.train.19.qn,matrix.test.19.qn),col=c(pData.train.19$IVIM_Color,pData.test.19$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

#### 7.2 Addon COMBAT batch correction of test set ########################################################### 
##############################################################################################################
matrix.test.19.bapred  <- t(combatbaaddon(combat.train.19, t(matrix.test.19.qn), batch = as.factor(pData.test.19$Batch)))

plot(Rtsne(t(cbind(matrix.train.19.bapred,matrix.test.19.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.19$Design_Color,pData.test.19$IVIM_Color), pch=16, cex=1.3)

#### 7.3 subset test matrix for the optimal variables ########################################################
##############################################################################################################
matrix.test.19.bapred.GA   <- t(matrix.test.19.bapred[GA.19$optVariables,])
matrix.test.19.bapred.full <- t(matrix.test.19.bapred[colnames(matrix.train.19),])

#### 7.4 PCA on best predictors found by genetic algorithm ###################################################
##############################################################################################################
pca.train <- prcomp(matrix.train.19.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.19$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.19$Design), col=unique(pData.train.19$Design_Color), pch=16, bty="n", cex=1.3)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.19.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.19$Design_Color,pData.test.19$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.19), cex= 0.4, pos=3, offset = 0.3) 

#### 7.5 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.19 <- predict(svmOpt.GA.19,matrix.test.19.bapred.GA, type = "prob")
Prediction_GA.19$Prediction_GA.19 <- ifelse(Prediction_GA.19$transforming>0.50,"transforming","untransforming")
Prediction_GA.19 <- cbind(pData.test.19[,c(1:3)],TrueLabel=pData.test.19$Class,Prediction_GA.19)
write.table(Prediction_GA.19, file = paste("TestSet019_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.19 <- predict(svmFull.19, matrix.test.19.bapred.full, type = "prob")
Prediction_SVM_full.19$Prediction_SVM_full <- ifelse(Prediction_SVM_full.19$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.19 <- cbind(pData.test.19[,c(1:3)],TrueLabel=pData.test.19$Class,Prediction_SVM_full.19)
write.table(Prediction_SVM_full.19, file = paste("TestSet019_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 8 Performance of optVars Classifier on  TestSet 019 #####################################################
#############################################################################################################

#### 8.1 Confusion matrix ###################################################################################  
sink("TestSet019_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.19$Prediction_GA.19), as.factor(Prediction_GA.19$TrueLabel))
sink()

sink("TestSet019_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.19$Prediction_SVM_full), as.factor(Prediction_SVM_full.19$TrueLabel))
sink()

#### 8.2 ROC on probability "transforming" TestSet 019 ########################################################

Prediction_GA.19$Class <- as.factor(ifelse(Prediction_GA.19$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.19 <- roc(Prediction_GA.19$Class,                    # response vector (factor or character)
                 Prediction_GA.19$transforming,             # predictor vector (numeric)
                 percent=TRUE, levels=c("nontransforming","transforming"),
                 plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                 print.auc=T)

Prediction_SVM_full.19$Class <- as.factor(ifelse(Prediction_SVM_full.19$TrueLabel == "transforming","transforming","nontransforming"))
roc.full.19<- roc(Prediction_SVM_full.19$Class,                    # response vector (factor or character)
                  Prediction_SVM_full.19$transforming,             # predictor vector (numeric)
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

#### 8.3 Precision-Recall curve on probability "transforming" TestSet 019 ####################################
Prediction_GA.19$Class_Code <- ifelse(Prediction_GA.19$TrueLabel == "transforming",1,0)
pr.GA.19 <- pr.curve(scores.class0 = Prediction_GA.19$transforming , weights.class0 = Prediction_GA.19$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.GA.19, rand.plot = TRUE, legend = F, color = 1,main = "")

Prediction_SVM_full.19$Class_Code <- ifelse(Prediction_SVM_full.19$TrueLabel == "transforming",1,0)
pr.full.19 <- pr.curve(scores.class0 = Prediction_SVM_full.19$transforming , weights.class0 = Prediction_SVM_full.19$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.full.19, rand.plot = TRUE, legend = F, color = 1,main = "")

#############################################################################################################################################
#############################################################################################################################################
#### FINAL MODEL WITHOUT INDEPENDENT TEST SET: Performance estimated only from resampling ###################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1.  quantile normalization & combat correction of training set only #####################################
##############################################################################################################
boxplot(log2(SAGA_RAW$E), col=pData$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.19 quantile normalization  ############################################################################# 
qunorm.train.FINAL    <- qunormtrain(t(SAGA_RAW$E))
matrix.train.FINAL.qn <- log2(t(qunorm.train.FINAL$xnorm))
matrix.train.FINAL.qn <- avereps(matrix.train.FINAL.qn, ID= SAGA_RAW$genes$ProbeName)  
colnames(matrix.train.FINAL.qn) <- row.names(pData)
boxplot(matrix.train.FINAL.qn, col=pData$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.19 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(matrix.train.FINAL.qn),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData$IVIM_Color, pch=16, cex=1.3) 

#### 2.19 COMBAT batch correction ############################################################################# 
batch.train.FINAL         <- as.factor(pData$Batch)                       
combat.train.FINAL        <- combatba(t(matrix.train.FINAL.qn), batch = batch.train.FINAL)
matrix.train.FINAL.batch  <- t(combat.train.FINAL$xadj)
colnames(matrix.train.FINAL.batch) <- row.names(pData)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.FINAL.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet FINAL ###################################################
##############################################################################################################
fselect.FINAL  <- genefilter(matrix.train.FINAL.batch, filterfun(f1))
summary(fselect.FINAL)
matrix.train.FINAL <-matrix.train.FINAL.batch[fselect.FINAL,]


##############################################################################################################
#### 2. SVM: FULL MODEL FINAL ################################################################################
##############################################################################################################
matrix.train.FINAL <- (t(matrix.train.FINAL))
labels.train.FINAL <- as.factor(pData$Class)

## create 200 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for comparison
set.seed(123)
index.FINAL <- createMultiFolds(labels.train.FINAL, k=10, times = 20)  

fullCtrl.FINAL <- trainControl(method = "repeatedcv", repeats = 20,index = index.FINAL,
                               summaryFunction = fiveStats,
                               returnData = TRUE,
                               returnResamp = "all",
                               savePredictions = "final",
                               classProbs = TRUE,
                               allowParallel = TRUE)

set.seed(721)
svmFull.FINAL <- train(matrix.train.FINAL,labels.train.FINAL,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.FINAL)

svmFull.FINAL  

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

system.time(rfe.FINAL  <- rfe(matrix.train.FINAL, labels.train.FINAL, 
                              sizes=FeatureNumbers,
                              rfeControl=outerctrl.FINAL,
                              metric = "Accuracy",
                              method="svmRadial",
                              tuneLength = 20,
                              trControl = innerctrl))

rfe.FINAL   # 26 optVars found by rfe
write.table(rfe.FINAL$results, file = "FinalSet154_Results_rfe.txt", sep="\t",col.names=NA)

# Figure 3h
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
sink("FinalSet152_Resamples_rfe vs full.txt", append = TRUE)
summary(rfeResamples.FINAL)
sink()

modelDifferences.FINAL <- diff(rfeResamples.FINAL)  # paired t-test for H0: difference = 0 between the different models. 
sink("FinalSet152_ModelDifferences_rfe vs full.txt", append = TRUE)
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
                            metric = "Accuracy",
                            method = "svmRadial",
                            # inner loop control for hyperparameter tuning
                            tuneLength = 12,
                            trControl = trainControl(method = "repeatedcv",
                                                     repeats = 2,
                                                     allowParallel = FALSE)))

### 4.3 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.FINAL  # 10 Features selected by GA
optVars.GA.FINAL <- Annotation[GA.FINAL$optVariables,]  # export optimal variables
write.table(optVars.GA.FINAL, file = "FinalSet152_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.FINAL <- GA.FINAL$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.FINAL <- arrange(performance.external.FINAL, Iter)

performance.FINAL <- GA.FINAL$ga$internal               # export internal accuracy (within the resample)
performance.FINAL$AccuracyExternal <- aggregate(performance.external.FINAL$Accuracy, by=list(Iter=performance.external.FINAL$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.FINAL <- data.frame(Iter = c(performance.FINAL$Iter,performance.FINAL$Iter), Accuracy = c(performance.FINAL$AccuracyExternal,performance.FINAL$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.FINAL <- GA.FINAL$averages[GA.FINAL$optIter,2]
accuracy.external.FINAL

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.FINAL <- subset(performance.external.FINAL,performance.external.FINAL$Iter == GA.FINAL$optIter)
accuracy.external.opt.FINAL <- accuracy.external.opt.FINAL$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.FINAL, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

### Figure 3i plot external accuracy over the iterations ################################################################
performance.aggr.FINAL <- subset(performance.long.FINAL,performance.long.FINAL$Group =="external")

ggplot(performance.aggr.FINAL, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T,colour = "#E8534F") +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())


##############################################################################################################
##### 5  PCA on optimal variables FINALSET152 ################################################################
##############################################################################################################

# PCA on 36,226 probes
pca.eset.batch <- prcomp(t(matrix.train.FINAL.batch),center = T, scale. = T)           
plot(pca.eset.batch$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(35,55, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

# PCA on 10 probes from SVM-GA
matrix.train.GA.FINAL  <- matrix.train.FINAL.batch[GA.FINAL$optVariables,]
pca.FINAL              <- prcomp(t(matrix.train.GA.FINAL),center = T, scale. = T)           
plot(pca.FINAL$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(2,4, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)



#############################################################################################################################################
#############################################################################################################################################
#### Analyze all results from the different Test/Training splits ############################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### 1. Summary of features, resampling- and test set accuracies for all 19 batchwise splits #################################################
#############################################################################################################################################
library(gmodels)

### 1. create data.frame for all 19 batchwise splits
Results <- data.frame(Testset = c(seq(1,19,by=1)), FeatureNumber_full = rep(NA,19), ResamplingAccuracy_full=rep(NA,19), 
                      FeatureNumber_rfe = rep(NA,19), ResamplingAccuracy_rfe = rep(NA,19), P_full_vs_rfe=rep(NA,19),
                      FeatureNumber_GA  = rep(NA,19), ResamplingAccuracy_GA = rep(NA,19),
                      Accuracy_TestSet_full = rep(NA,19), AUROC_full = rep(NA,19), AUPRC_full = rep(NA,19), 
                      Accuracy_TestSet_SAGA = rep(NA,19), AUROC_SAGA = rep(NA,19), AUPRC_SAGA = rep(NA,19),AUPRC_random = rep(NA,19))


### 2. number of features selected by IQR filter for each training/test split
for (i in 1:19) {Results[i,2] <- dim(eval(parse(text=paste("matrix.train.",i, sep=""))))[2] }

### 3. resampling accuracy for the full model / stored in svmFull.x$results dataframe
for (i in 1:19) { Results[i,3] <- round(eval(parse(text=paste("svmFull.",i, sep="")))$results[eval(parse(text=paste("svmFull.",i, sep="")))$results$C == eval(parse(text=paste("svmFull.",i, sep="")))$bestTune$C,"Accuracy"],3)}

### 4. number of features found by SVM-rfe for each training/test split
for (i in 1:19) {Results[i,4] <- eval(parse(text=paste("rfe.",i, sep="")))$optsize }

### 5. external resampling accuracy for the svm-rfe model
for (i in 1:19) { Results[i,5] <- round(eval(parse(text=paste("rfe.",i, sep="")))$results[eval(parse(text=paste("rfe.",i, sep="")))$results$Variables == eval(parse(text=paste("rfe.",i, sep="")))$optsize,5],3)}

### 6. p.value for H0 = resampling accuracy of full model = resampling accuracy of full model/ stored in rfeResamples.x$values dataframe
for (i in 1:16) {Results[i,6] <- eval(parse(text=paste("modelDifferences.",i, sep="")))$statistics$Accuracy[[1]][3]}
Results[17,10] <- NA
for (i in 18:19) {Results[i,6] <- eval(parse(text=paste("modelDifferences.",i, sep="")))$statistics$Accuracy[[1]][3]}

### 7. number of features found by SVM-GA 
for (i in c(1,3:19)) { Results[i,7] <- length(eval(parse(text=paste("GA.",i, sep="")))$optVariables)  }

### 8. Resampling accuracy for the GA model
for (i in c(1,3:19)) { Results[i,8] <- round(eval(parse(text=paste("accuracy.external.",i, sep=""))),3)}

### 9. TestSet accuracy for the full model
for (i in 1:19) {Results[i,9] <- round(confusionMatrix(as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$Prediction_SVM_full), as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$TrueLabel))$overall[1],3)}

### 10. TestSet AUROC for the full model
for (i in 1:19) {Results[i,10] <- round(eval(parse(text=paste("roc.full.",i, sep="")))$auc[1]/100,3) }

### 11. TestSet AUPRC for the full model
for (i in 1:19) {Results[i,11] <- round(eval(parse(text=paste("pr.full.",i, sep="")))$auc.integral,3)}

### 12. TestSet accuracy for SAGA = TestSet accuracy of rfe (TestSet 02; <10 features found by rfe) and GA otherwise (all other TestSets)
for (i in c(1,3:19)) { Results[i,12] <- round(confusionMatrix(as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))$TrueLabel))$overall[1],3)}
Results[2,12] <- round(confusionMatrix(as.factor(eval(parse(text=paste("Prediction_rfe.",2, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_rfe.",2, sep="")))$TrueLabel))$overall[1],3)

### 13. TestSet AUROC for SAGA
for (i in c(1,3:19)) { Results[i,13] <- round(eval(parse(text=paste("roc.GA.",i, sep="")))$auc[1]/100,3) }
Results[2,13] <- round(roc.rfe.2$auc[1]/100,3)

### 14. TestSet AUPRC for SAGA
for (i in c(1,3:19)) { Results[i,14] <- round(eval(parse(text=paste("pr.GA.",i, sep="")))$auc.integral,3)}
Results[2,14] <- round(pr.rfe.2$auc.integral,3)

### 15. TestSet AUPRC for a random classifier
for (i in c(1,3:19)) { Results[i,15] <- round(eval(parse(text=paste("pr.GA.",i, sep="")))$rand$auc.integral,3)}
Results[2,15] <- round(eval(parse(text=paste("pr.rfe.",2, sep="")))$rand$auc.integral,3)

write.table(Results, file = "Results_19 TestSets_batchwise_bapred upfront.txt", sep="\t",col.names=NA)


#############################################################################################################################################
### 2. aggregate Test set performances over the 19 batchwise splits #########################################################################
#############################################################################################################################################
colnames(Prediction_rfe.2)[7] <- "Prediction" 

colnames(Prediction_GA.1)[7] <- "Prediction" 
colnames(Prediction_GA.3)[7] <- "Prediction" 
colnames(Prediction_GA.4)[7] <- "Prediction" 
colnames(Prediction_GA.5)[7] <- "Prediction" 
colnames(Prediction_GA.6)[7] <- "Prediction" 
colnames(Prediction_GA.7)[7] <- "Prediction" 
colnames(Prediction_GA.8)[7] <- "Prediction" 
colnames(Prediction_GA.9)[7] <- "Prediction" 
colnames(Prediction_GA.10)[7] <- "Prediction" 
colnames(Prediction_GA.11)[7] <- "Prediction" 
colnames(Prediction_GA.12)[7] <- "Prediction" 
colnames(Prediction_GA.13)[7] <- "Prediction" 
colnames(Prediction_GA.14)[7] <- "Prediction" 
colnames(Prediction_GA.15)[7] <- "Prediction" 
colnames(Prediction_GA.16)[7] <- "Prediction" 
colnames(Prediction_GA.17)[7] <- "Prediction" 
colnames(Prediction_GA.18)[7] <- "Prediction" 
colnames(Prediction_GA.19)[7] <- "Prediction" 

#############################################################################################################################################
### 2.4 SAGA = Compound TestSet performances = 1xSVM-rfe (TestSet 02: RFE < 10 predictors) + 18 x SVM_GA models (RFE >10 predictors) ########
#############################################################################################################################################
Predictions_SAGA   <- rbind(Prediction_GA.1,Prediction_rfe.2, Prediction_GA.3,Prediction_GA.4,Prediction_GA.5,
                            Prediction_GA.6,Prediction_GA.7,Prediction_GA.8,Prediction_GA.9,Prediction_GA.10,
                            Prediction_GA.11,Prediction_GA.12,Prediction_GA.13,Prediction_GA.14,Prediction_GA.15,
                            Prediction_GA.16,Prediction_GA.17,Prediction_GA.18,Prediction_GA.19)


Predictions_SVM_full <- rbind(Prediction_SVM_full.1,Prediction_SVM_full.2, Prediction_SVM_full.3,Prediction_SVM_full.4,Prediction_SVM_full.5,
                                 Prediction_SVM_full.6,Prediction_SVM_full.7,Prediction_SVM_full.8,Prediction_SVM_full.9,Prediction_SVM_full.10,
                                 Prediction_SVM_full.11,Prediction_SVM_full.12,Prediction_SVM_full.13,Prediction_SVM_full.14,Prediction_SVM_full.15,
                                 Prediction_SVM_full.16,Prediction_SVM_full.17,Prediction_SVM_full.18,Prediction_SVM_full.19)


#############################################################################################################################################
### 2.5. Confusion matrices of aggregated TestSets ##########################################################################################
#############################################################################################################################################

sink("AllTestSets_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_SAGA$Prediction), as.factor(Predictions_SAGA$TrueLabel))
sink()

sink("AllTestSetsConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_SVM_full$Prediction_SVM_full), as.factor(Predictions_SVM_full$TrueLabel))
sink()

#############################################################################################################################################
### 2.6. AUROC and PRROC SAGA vs IVIM (all available IVIMs ): ###############################################################################
#############################################################################################################################################

### 2.6.1 read in all available IVIM data (from Fig. 1d, n=502) #############################################################################
#############################################################################################################################################
IVIM  <- read.delim("IVIM_MTT_ROC.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE)
IVIM$Classification_Code <- as.factor(IVIM$Classification_Code)
IVIM$Prediction <- ifelse(IVIM$MTT_Score>=3,"transforming","untransforming")  
IVIM$TrueLabel  <- ifelse(IVIM$Classification_Code==1,"transforming","untransforming")
confusionMatrix(as.factor(IVIM$Prediction), as.factor(IVIM$TrueLabel))

### 2.6.2. IVIM vs SAGA AUROC (for all vectors)  ############################################################################################
#############################################################################################################################################
IVIM.total <- roc(IVIM$Classification_Code,    # response vector (factor or character)
                  IVIM$MTT_Score,              # predictor vector (numeric)
                  percent=TRUE, smooth = F,
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                  col = "#8285BC", grid=F, lwd = 3, cex.lab=1.5, 
                  print.auc=T, print.thres = 3 )

roc.SAGA   <- roc( Predictions_SAGA$Class,                    
                   Predictions_SAGA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                   print.auc=T, add = T,print.thres="best")

### 2.6.3. compare ROC curves between SAGA and IVIM #########################################################################################
#############################################################################################################################################
roc.test(roc.SAGA, IVIM.total, alternative = "greater")   # p-value = 1.717e-05


### 2.6.4. IVIM vs SAGA on LTR.RV.SF only: Supplementary Figure 6b ###################################################################################################
#############################################################################################################################################
selected  <- c("MOCK","LTR.RV.SF")  # select Mock controls and LTR.SF in IVIM data
IVIM.sel1 <- subset(IVIM,IVIM$Vector %in% selected)
confusionMatrix(as.factor(IVIM.sel1$Prediction), as.factor(IVIM.sel1$TrueLabel))

IVIM.LTR.SF <- roc(IVIM.sel1$Classification_Code,    # response vector (factor or character)
                   IVIM.sel1$MTT_Score,              # predictor vector (numeric)
                   percent=TRUE, smooth = F,
                   plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                   col = "#8285BC", grid=F, lwd = 3, cex.lab=1.5, 
                   print.auc=T,print.thres=3)

# select Mock controls and LTR.SF in SAGA Predictions in the TestSets
Predictions_SAGA_LTR.SF <- Predictions_SAGA[(grepl("Mock",Predictions_SAGA$Name) | grepl("RV.SF",Predictions_SAGA$Name)),]

sink("ConfusionMatrix_SAGA_LTR.SF vs Mock.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_SAGA_LTR.SF$Prediction), as.factor(Predictions_SAGA_LTR.SF$TrueLabel))
sink()

ROC.SAGA.LTR.SF   <- roc(Predictions_SAGA_LTR.SF$Class,                    
                         Predictions_SAGA_LTR.SF$transforming,             
                         percent=TRUE, levels=c("nontransforming","transforming"),
                         plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                         print.auc=T, add = T,print.thres=0.5)

roc.test(ROC.SAGA.LTR.SF, IVIM.LTR.SF, alternative = "greater" )   #p-value = 6.029e-07

### 2.6.5. IVIM vs SAGA non-LTR vectors Supplementary Figure 6c #############################################################################
#############################################################################################################################################
IVIM.sel2 <- subset(IVIM, IVIM$Vector!= "LTR.RV.SF")
confusionMatrix(as.factor(IVIM.sel2$Prediction), as.factor(IVIM.sel2$TrueLabel))

ROC.IVIM.other <- roc(IVIM.sel2$Classification_Code,    # response vector (factor or character)
                      IVIM.sel2$MTT_Score,          # predictor vector (numeric)
                      percent=TRUE, smooth = F,
                      plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                      col = "#8285BC", grid=F, lwd = 3, cex.lab=1.5, 
                      print.auc=T,print.thres=3)

Predictions_SAGA_other <- Predictions_SAGA[(!grepl("RV.SF",Predictions_SAGA$Name)),]
confusionMatrix(as.factor(Predictions_SAGA_other$Prediction), as.factor(Predictions_SAGA_other$TrueLabel))

ROC.SAGA.other <- roc(Predictions_SAGA_other$Class,                    
                      Predictions_SAGA_other$transforming,             
                      percent=TRUE, levels=c("nontransforming","transforming"),
                      plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                      print.auc=T, add = T,print.thres=0.5)

roc.test(ROC.SAGA.other, ROC.IVIM.other, alternative = "greater" )  # p = 0.002314


### 2.6.6. IVIM vs SAGA PRROC ###############################################################################################################
#############################################################################################################################################

### 2.6.6.1 for all vectors: #######################################################################################################
AUPRC.SAGA <- pr.curve(scores.class0 = Predictions_SAGA$transforming , weights.class0 = Predictions_SAGA$Class_Code, curve = TRUE,rand.compute = T)
plot(AUPRC.SAGA, rand.plot = TRUE, legend = F, color = "#F9C35F", main = "",auc.main = T)  # AUPRC = 0.944  / AUC random = 0.422

IVIM$Class_Code <- ifelse(IVIM$TrueLabel == "transforming",1,0)
AUPRC.IVIM <- pr.curve(scores.class0 = IVIM$MTT_Score, weights.class0 = as.numeric(IVIM$Class_Code), curve = TRUE,rand.compute = T)
plot(AUPRC.IVIM, rand.plot = TRUE, legend = F, color = "#8285BC", add = TRUE) # AUPRC = 0.892  / AUC random = 0.582

### 2.6.2.2 for LTR.SF #########################################################################################################################
AUPRC.SAGA.LTR <- pr.curve(scores.class0 = Predictions_SAGA_LTR.SF$transforming , weights.class0 = Predictions_SAGA_LTR.SF$Class_Code, curve = TRUE,rand.compute = T)
plot(AUPRC.SAGA.LTR, rand.plot = TRUE, legend = F, color = "#F9C35F", main = "")  # AUPRC = 0.9993  / AUC random = 0.521

IVIM.sel1$Class_Code <- ifelse(IVIM.sel1$TrueLabel == "transforming",1,0)
AUPRC.IVIM.LTR <- pr.curve(scores.class0 = IVIM.sel1$MTT_Score, weights.class0 = as.numeric(IVIM.sel1$Class_Code), curve = TRUE,rand.compute = T)
plot(AUPRC.IVIM.LTR, rand.plot = TRUE, legend = F, color = "#8285BC", add = TRUE) # AUPRC = 0.94 / AUC random = 0.64

### 2.6.2.3 for non-LTR.SF  Supplementary Figure 6c #################################################################################################################
AUPRC.SAGA.other <- pr.curve(scores.class0 = Predictions_SAGA_other$transforming , weights.class0 = Predictions_SAGA_other$Class_Code, curve = TRUE,rand.compute = T)
plot(AUPRC.SAGA.other, rand.plot = TRUE, legend = F, color = "#F9C35F", main = "")  # AUPRC = 0.79  / AUC random = 0.24

IVIM.sel2$Class_Code <- ifelse(IVIM.sel2$TrueLabel == "transforming",1,0)
AUPRC.IVIM.other <- pr.curve(scores.class0 = IVIM.sel2$MTT_Score, weights.class0 = as.numeric(IVIM.sel2$Class_Code), curve = TRUE,rand.compute = T)
plot(AUPRC.IVIM.other, rand.plot = TRUE, legend = F, color = "#8285BC", add = TRUE) # AUPRC = 0.59 / AUC random = 0.28


#############################################################################################################################################
### 2.7. Toplist of predictors chosen during the 19 iterations / Table 1 ####################################################################
#############################################################################################################################################

all.optVars <- c(optVars.GA.1$GeneSymbol_FINAL,optFeatures.rfe.2$GeneSymbol_FINAL,optVars.GA.3$GeneSymbol_FINAL,optVars.GA.4$GeneSymbol_FINAL,optVars.GA.5$GeneSymbol_FINAL,
                 optVars.GA.6$GeneSymbol_FINAL,optVars.GA.7$GeneSymbol_FINAL,optVars.GA.8$GeneSymbol_FINAL,optVars.GA.9$GeneSymbol_FINAL,optVars.GA.10$GeneSymbol_FINAL,
                 optVars.GA.11$GeneSymbol_FINAL,optVars.GA.12$GeneSymbol_FINAL,optVars.GA.13$GeneSymbol_FINAL,optVars.GA.14$GeneSymbol_FINAL,optVars.GA.15$GeneSymbol_FINAL,
                 optVars.GA.16$GeneSymbol_FINAL,optVars.GA.17$GeneSymbol_FINAL,optVars.GA.18$GeneSymbol_FINAL,optVars.GA.19$GeneSymbol_FINAL)

b <- table(all.optVars)
b <- b[order(b,decreasing = T)]

write.table(b, file="Toplist_optimal predictors_batchwise_bapred upfront.txt", sep="\t")

