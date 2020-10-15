library(limma)
library(genefilter)
library(Rtsne)
library(RColorBrewer)
library(caret)
library(kernlab)
library(ggplot2)
library(dplyr)
library(pROC)
library(PRROC)
library(sva)
library(bapred)
library(doMC)
registerDoMC(cores = 10)    


#############################################################################################################################################
#############################################################################################################################################
#### A) create EListRAW of 154 arrays / 19 batches ##########################################################################################
#############################################################################################################################################
#############################################################################################################################################

# 1. all 169 arrays were read in and combined into a EListRaw object without further modification
# 2. 15 samples with unknown ground truth were subsequently removed from the dataset, resulting in 154 assays
# 3. IVIM ID 180523 was splitted due to severe class imbalance, the two mock controls were duplicated and used for #180523A and #180523B 
# 4. One separate IVIM assay is treated as one batch except for IVIM #160706 and IVIM #160525 which comprise a single batch (7) since RNA isolation + Microarray processing was performed together
# 5. quantile normalization, averaging, log2 and batch correction are performed on the training and test sets separately using the bapred packages to prevent any data leakage between test and training sets

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
SAGA_RAW <- SAGA_RAW.full[,row.names(pData)]             # EListRaw object with 154 samples in 19 batches /  157,712 probes (39,428 unique probes in quadruplicate)       

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

#### 2.1 quantile normalization  ############################################################################# 
RMA_train.1 <- normalizeBetweenArrays(RAW_train.1,method="quantile")      # limma: EList object with quantile normalized and log2 data
RMA_train.1 <- avereps(RMA_train.1, ID= RMA_train.1$genes$ProbeName)      # limma: EList object averaged over ProbeIDs    
boxplot(RMA_train.1$E, col=pData.train.1$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.1$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$IVIM_Color, pch=16, cex=1.3) 
legend(15,26, legend=unique(pData.train.1$IVIM_ID), col=unique(pData.train.1$IVIM_Color), pch=16, bty="n", cex=0.6)

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.1 <- pData.train.1$Batch-1                       
modcombat     <- model.matrix(~1, data=pData.train.1)         
matrix.train.1.batch <- ComBat(dat=RMA_train.1$E, batch=batch.train.1, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.1.batch <- matrix.train.1.batch[row.names(Annotation.known),] # filter for 36,226 annotated probes

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.1.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$IVIM_Color, pch=16, cex=1.3) 

set.seed(12)
plot(Rtsne(t(matrix.train.1.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$Design_Color, pch=16, cex=1.3) 
legend(-16,-6, legend=unique(pData.train.1$Design), col=unique(pData.train.1$Design_Color), pch=16, bty="n", cex=1.3)

##############################################################################################################
#### 3. nonspecific feature prefiltering of training set #####################################################
##############################################################################################################
f1       <- function(x) (IQR(x) > 0.8)    
fselect.1  <- genefilter(matrix.train.1.batch, filterfun(f1))
summary(fselect.1)                                  # 1226 predictors selected with interquartile range of log2-int >0.8
matrix.train.1 <-matrix.train.1.batch[fselect.1,]   # subset to 1226 predictors

##############################################################################################################
#### 4. SVM: FULL MODEL (1226 predictors) CV on training set #################################################
##############################################################################################################
matrix.train.1 <- (t(matrix.train.1))
labels.train.1 <- as.factor(pData.train.1$Class)

# calculate performance measures (accuracy, sensitivity, specificity, ROC) of external resamples
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))   

## create 200 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for direct model comparison
set.seed(123)
index.1 <- createMultiFolds(labels.train.1, k=10, times = 20)  

## define parameters of the train function  
fullCtrl.1 <- trainControl(method = "repeatedcv",repeats = 20,  # 10foldCVn20  
                           index = index.1,                     # fixed index to compare with the rfe model
                           summaryFunction = fiveStats,         # define summary function
                           classProbs = TRUE,                   # calculate class probabilities
                           allowParallel = TRUE)                # enable parallel computation 

## build full model on complete training data with all predictors    
set.seed(721)
system.time(svmFull.1 <- train(matrix.train.1,labels.train.1,  # define training set  
                               method = "svmRadial",           # support vector machine with radial kernel
                               metric = "Accuracy",            # use accuracy to select the best model
                               tuneLength = 20,                # number of cost values to test (Caret creates a range of values and uses a single value of sigma that is calculated internally with kernlab “sigest” function) 
                               trControl = fullCtrl.1))
svmFull.1  

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

rfe.1   # yields 14 features
write.table(rfe.1$results, file = "TrainingSet01_Results_rfe.txt", sep="\t",col.names=NA)
optFeatures.rfe.1 <- cbind(rfe.1$optVariables, Annotation[rfe.1$optVariables,])
write.table(optFeatures.rfe.1, file = "TrainingSet01_OptVars_rfe.txt", sep="\t",col.names=NA)

# Plot Accuracy over FeatureNumber
trellis.par.set(caretTheme())
plot(rfe.1, type = c("g", "o"))
plot(rfe.1, type = c("g", "o"), xlim = c(0,62))

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
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 01 #########################
################################################################################################

matrix.train.rfe.1 <- matrix.train.1[,rfe.1$optVariables]   # subset for the 14 optVars from rfe

#### 6.0 Set global SVM-GA Parameters:##########################################################
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



#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2334)
index.GA.1 <- createMultiFolds(labels.train.1, k=10, times = 5)  

# set all seeds for running the genetic algorithm in parallel over the 50 different resamples ##
set.seed(123)
seeds.GA <- vector(mode = "integer", length = length(index.GA.1)+1)    # B+1 elements where B is the number of resamples = 51
for(i in 1:length(index.GA.1)+1) seeds.GA[[i]] <- sample.int(10000, 1)

outerctrl.GA.1 <- gafsControl(functions = svmGA,
                              method = "repeatedcv", repeats = 5,
                              index = index.GA.1,                                       
                              seeds = seeds.GA,                                      
                              returnResamp="all", 
                              verbose = TRUE,
                              maximize = c(internal = TRUE,
                                           external = TRUE),
                              allowParallel = TRUE)                                  

#### 6.2 run GA  ###############################################################################
################################################################################################

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

### 6.3 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.1                                            # 8 features selected at iteration 26
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



##############################################################################################################
#### 7.0 reviewers approach w/o batch correction ############################################################# 
##############################################################################################################

#### 7.1 Subset normalized and batch-corrected training set from step 2 to 8 optimal predictors from step 3 ## 
##############################################################################################################
matrix.train.1.GA   <- t(matrix.train.1.batch[GA.1$optVariables,])

#### 7.2 PCA of training set from step 2 reduced to 8 optimal predictors from step 3 ######################### 
##############################################################################################################
pca.train <- prcomp(matrix.train.1.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.1$Design_Color), cex=1.8, asp=1)
legend(-4.2,4, legend=unique(pData.train.1$Design), col=unique(pData.train.1$Design_Color), pch=16, bty="n", cex=1.3)

summary(pca.train)  

#### 7.3 train SVM on training set from step 2 reduced to 8 optimal predictors from step 3  ##################
##############################################################################################################
set.seed(721)
SVM.8optPred.1  <- train(matrix.train.1.GA,labels.train.1,
                         method = "svmRadial",
                         metric = "Accuracy",
                         tuneLength = 20,         # number of cost values to test (Caret creates a range of values
                                                  # and uses a analytically determined value of sigma from kernlab's “sigest” 
                         trControl = fullCtrl.1)  # use same parameters as for the SVMfull model:  
                                                  # 20 times repeated 10-fold cross-validation for tuning
                                                  # fixed resamples (index.1)
SVM.8optPred.1  # best CV accuracy = 0.9557232

#### 7.4 prepare test set (quantile normalize and average within the batch)  ##################################
###############################################################################################################
RMA_test.1 <- normalizeBetweenArrays(RAW_test.1,method="quantile")     # limma: EList object with quantile normalized and log2 data
RMA_test.1 <- avereps(RMA_test.1, ID= RMA_test.1$genes$ProbeName)      # limma: EList object averaged over ProbeIDs    
RMA_test.1 <- RMA_test.1$E
RMA_test.1 <- RMA_test.1[row.names(Annotation.known),]                 # filter for 36,226 annotated probes

boxplot(log2(RAW_test.1$E), col=pData.test.1$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    
boxplot(RMA_test.1, col=pData.test.1$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 7.5 visualize test set and training set in global gene expression space of 36,226 predictors ##############
################################################################################################################
plot(Rtsne(t(cbind(matrix.train.1.batch,RMA_test.1)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
           col=c(pData.train.1$Design_Color,pData.test.1$IVIM_Color), pch=16, cex=1.3)

pca <- prcomp(t(cbind(matrix.train.1.batch,RMA_test.1)), center = T, scale. = T)           
plot(pca$x, pch=16, col=c(pData.train.1$Design_Color,pData.test.1$IVIM_Color), cex=1.8, asp=1)
summary(pca)    # 33.4 % of variance is batch effect, 8.1 % of variance is class-signal

#### 7.6 reduced test set to 8 optimal predictors and visualize test set and training set  #####################
################################################################################################################
matrix.test.1.reviewer.GA <- t(RMA_test.1[GA.1$optVariables,]) # reduce test set to 8 optimal probes
dim(matrix.test.1.reviewer.GA)			                           # test set: 7 samples x 8 predictors

## project test samples into PCA plot spanned by training set w/o altering training set PCA #################### 
## (prediction by multipling the test vectors with eigenvectors (loadings) from the training set PCA )
coord.pca.test <- predict(pca.train, newdata = matrix.test.1.reviewer.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.1$Design_Color,pData.test.1$IVIM_Color), cex=1.5, asp=1)
text(coord.pca.test, labels=row.names(pData.test.1), cex= 0.4, pos=3, offset = 0.3) 


#### 7.6 predict class labels of test set ######################################################################
################################################################################################################
Prediction_reviewer.1 <- predict(SVM.8optPred.1,matrix.test.1.reviewer.GA, type = "prob")
Prediction_reviewer.1$Prediction_GA.1 <- ifelse(Prediction_reviewer.1$transforming>0.50,"transforming","untransforming")
Prediction_reviewer.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_reviewer.1)
write.table(Prediction_GA.1, file = paste("TestSet01_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)
Prediction_reviewer.1

#### 8 Performance of optVars Classifier on  TestSet 01 #####################################################
#############################################################################################################

sink("TestSet01_ConfusionMatrix_GA_Reviewer.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_reviewer.1$Prediction_GA.1), as.factor(Prediction_reviewer.1$TrueLabel))
sink()

#### 7.2 ROC on probability "transforming" TestSet 01 ########################################################

Prediction_reviewer.1$Class <- as.factor(ifelse(Prediction_reviewer.1$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA.1 <- roc(Prediction_reviewer.1$Class,                    # response vector (factor or character)
                Prediction_reviewer.1$transforming,             # predictor vector (numeric)
                percent=TRUE, levels=c("nontransforming","transforming"),
                plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                print.auc=T)


#### 7.2 Precision-Recall curve on probability "transforming" TestSet 01 ####################################
Prediction_reviewer.1$Class_Code <- ifelse(Prediction_reviewer.1$TrueLabel == "transforming",1,0)
pr.rfe.1 <- pr.curve(scores.class0 = Prediction_reviewer.1$transforming , weights.class0 = Prediction_reviewer.1$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.rfe.1,rand.plot = TRUE, legend = F, color = 1, main = "")



#### 7.7 use different models on the training set to predict the test set ######################################
################################################################################################################

## Random Forest 
set.seed(721)
rfGA.1 <- train(matrix.train.1.bapred.GA,labels.train.1,
              method = "rf",
              metric = "Accuracy",
              ntree = 5000,
              trControl = fullCtrl.1)
rfGA.1 #best accuracy = 0.9437738 / mtry = 2 selected predictors

Prediction_rfGA <- predict(rfGA.1, matrix.test.1.reviewer.GA, type = "prob")
Prediction_rfGA$Prediction_rfGA <- ifelse(Prediction_rfGA$transforming>0.50,"transforming","untransforming")
Prediction_rfGA <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_rfGA)
Prediction_rfGA

## Generilized Linear model
lrGA <- train(matrix.train.1.bapred.GA,labels.train.1,
              method = "glm",
              metric = "Accuracy",
              trControl = fullCtrl.1)
lrGA  #best accuracy = 0.934

Prediction_lrGA <- predict(lrGA, matrix.test.1.reviewer.GA, type = "prob")
Prediction_lrGA$Prediction_lrGA <- ifelse(Prediction_lrGA$transforming>0.50,"transforming","untransforming")
Prediction_lrGA <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_lrGA)
Prediction_lrGA

## Linear discriminant analysis
ldaGA <- train(matrix.train.1.bapred.GA,labels.train.1,
              method = "lda",
              metric = "Accuracy",
              trControl = fullCtrl.1)
ldaGA  #best accuracy = 0.9106012

Prediction_ldaGA <- predict(ldaGA, matrix.test.1.reviewer.GA, type = "prob")
Prediction_ldaGA$Prediction_ldaGA <- ifelse(Prediction_ldaGA$transforming>0.50,"transforming","untransforming")
Prediction_ldaGA <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_ldaGA)
Prediction_ldaGA

## KNN
set.seed(721)
knnGA <- train(matrix.train.1.bapred.GA,labels.train.1,
               method = "knn",
               metric = "Accuracy",
               tuneLength = 20,
               trControl = fullCtrl.1)
knnGA #best accuracy = 0.9318254

Prediction_knnGA <- predict(knnGA, matrix.test.1.reviewer.GA, type = "prob")
Prediction_knnGA$Prediction_knnGA <- ifelse(Prediction_knnGA$transforming>0.50,"transforming","untransforming")
Prediction_knnGA <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_knnGA)
Prediction_knnGA

## qda
qdaGA <- train(matrix.train.1.bapred.GA,labels.train.1,
               method = "qda",
               metric = "Accuracy",
               trControl = fullCtrl.1)
qdaGA #best accuracy = 0.9187143

Prediction_qdaGA <- predict(qdaGA, matrix.test.1.reviewer.GA, type = "prob")
Prediction_qdaGA$Prediction_qdaGA <- ifelse(Prediction_qdaGA$transforming>0.50,"transforming","untransforming")
Prediction_qdaGA <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_qdaGA)
Prediction_qdaGA


#fdaTest <- roc(testing$Class, 
predict(fdaModel, testing, type = "prob")[,1], 
levels = rev(levels(testing$Class)))
#fdaTest



