#############################################################################################################################################
#################################### Construction of SAGA-SVM Classifier at IQR 0.8 #########################################################
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
#############################################################################################################################################
#### I. Split 01 / TestSet 01 ###############################################################################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. Divide into TrainingSet 01 and TestSet 01 ############################################################
##############################################################################################################
set.seed(123)
split.1        <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.1 <- eset.batch[, split.1] # training matrix 107 samples x 36,226 predictors 
matrix.test.1  <- eset.batch[,-split.1] # test matrix 45 samples x 36,226 predictors 
pData.train.1  <- pData[split.1,]
pData.test.1   <- pData[-split.1,]

table(pData.train.1$Design)             # sample distribution in the training set
table(pData.test.1$Design)              # sample distribution in the test set 

##############################################################################################################
#### 2. nonspecific feature prefiltering  ####################################################################
##############################################################################################################
f1       <- function(x) (IQR(x) > 0.8)    
fselect.1  <- genefilter(matrix.train.1, filterfun(f1))
summary(fselect.1)                            # 1151 genes selected with interquartile range of log2-int >0.8
matrix.train.1 <-matrix.train.1[fselect.1,]   # subset 

##############################################################################################################
#### 3. SVM: FULL MODEL (1151 predictors) ####################################################################
##############################################################################################################
matrix.train.1 <- (t(matrix.train.1))
labels.train.1 <- as.factor(pData.train.1$Class)

# calculate accuracy, sensitivity, specificity, ROC and  kappa of external resamples
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
svmFull.1 <- train(matrix.train.1,labels.train.1,              # define training set  
                   method = "svmRadial",                       # support vector machine with radial kernel
                   metric = "Accuracy",                        # use accuracy to select the best model
                   tuneLength = 20,                            # number of cost values to test (Caret creates a range of values and uses a single value of sigma that is calculated internally with kernlab “sigest” function) 
                   trControl = fullCtrl.1)

svmFull.1  
svmFull.1$results  

##############################################################################################################
#### 4. SVM-RFE on TrainingSet 01 ############################################################################
##############################################################################################################

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

##  set predictos subset sizes to test: 1,2,…,40,45,50,60,….,500,1151  = 52 subsets in total
FeatureNumbers <- c(seq(1,40,by=1),45,50,60,70,80,90,100,200,300,400,500)                 

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

#### 4.2 Parameters for inner (nested) cross-validation loop for hyperparameter tuning within each resample and for each subset size #######
##############################################################################################################
innerctrl <- trainControl(method = "repeatedcv",repeats = 3,           # 10CVn3 = 30 resamples
                          verboseIter = FALSE,
                          classProbs = TRUE,
                          allowParallel = FALSE)                      

#### 4.3.SVM-RFE #############################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning x 20 values for cost parameter = 6,240,000 models
# this will take around 3-4 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.1  <- rfe(matrix.train.1, labels.train.1, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.1,
                          metric = "Accuracy",
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.1
write.table(rfe.1$results, file = "TrainingSet01_Results_rfe.txt", sep="\t",col.names=NA)

optFeatures.rfe.1 <- cbind(rfe.1$optVariables, Annotation[rfe.1$optVariables,])
write.table(optFeatures.rfe.1, file = "TrainingSet01_OptVars_rfe.txt", sep="\t",col.names=NA)

# Plot Accuracy over FeatureNumber
trellis.par.set(caretTheme())
plot(rfe.1, type = c("g", "o"))
plot(rfe.1, type = c("g", "o"), xlim = c(0,61))

#### 4.4.compare resampling performances between full model and rfe ##########################################
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
legend(2,-1, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.1$Design_Color), pch=16, bty="n", cex=0.8)

##### compared to PCA on all variables ########################################################
matrix.full.1 <- rbind(matrix.train.1, t(matrix.test.1[colnames(matrix.train.1),]))
pca.full.1    <- prcomp(matrix.full.1)           
plot(pca.full.1$x, pch=16, col=pData.opt.1$Design_Color, cex=1.5, asp=1)
legend(25,-20, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.1$Design_Color), pch=16, bty="n", cex=0.8)


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

# note: these resampling values seem to be better than the external resampling results (above) due to selection bias, 
#  since the error rates are based on the SVM model after the optimal features for exactly this training set have been selected by rfe.
  

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

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################
sink("TestSet01_ConfusionMatrix_optVars.txt", append = TRUE)
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


################################################################################################
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION ########################################
################################################################################################
matrix.train.rfe.2 <- matrix.train.2[,rfe.2$optVariables]   # subset training matrix to optimal variables from SVM-RFE

#### 4.1 Set global SVM-GA Parameters:##########################################################
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

#### 4.2 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(1714)
index.GA.2 <- createMultiFolds(labels.train.2, k=10, times = 5)        # 10CVn5 crossvalidation = 50 resamples for external prediction 

# set all seeds for running the genetic algorithm in parallel over the 50 different resamples ##
set.seed(123)
seeds.GA <- vector(mode = "integer", length = length(index.GA.2)+1)    # B+1 elements where B is the number of resamples = 51
for(i in 1:length(index.GA.2)+1) seeds.GA[[i]] <- sample.int(10000, 1)

outerctrl.GA.2 <- gafsControl(functions = svmGA,                 # define helper functions
                            method = "repeatedcv", repeats = 5,  # define cv-scheme
                            index = index.GA.2,                  # set index for cv
                            seeds = seeds.GA,                    # set all seeds                    
                            returnResamp="all", 
                            verbose = TRUE,
                            maximize = c(internal = TRUE,        # maximize internal and external accuracy
                                         external = TRUE),
                            allowParallel = TRUE)                # allow feature selection in parallel in the outer loop                 

#### 4.3 run GA  ###############################################################################
################################################################################################

# This will take 4-5 hrs on AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM
system.time(GA.2<- gafs(matrix.train.rfe.2, labels.train.2, 
                        iters = 40,                           # run 40 iterations 
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  # define parameters of GA (default)
                        gafsControl = outerctrl.GA.2,         # parameters for outer layer
                        metric = "Accuracy",
                        method = "svmRadial",
                        # inner loop control for hyperparameter tuning 
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 4.4 analyze results from Genetic Algorithm #################################################
################################################################################################
GA.2
optVars.GA.2 <- Annotation[GA.2$optVariables,]  # export optimal variables
write.table(optVars.GA.2, file = "TrainingSet02_optVars_GeneticAlgorithm.txt", sep="\t",col.names=NA)

performance.external.2 <- GA.2$external         # export external accuracy (resamples prediction of held-outs of each resample)
performance.external.2 <- arrange(performance.external.2, Iter)
performance.2          <- GA.2$ga$internal               # export internal accuracy (within the resample)
performance.2$AccuracyExternal <- aggregate(performance.external.2$Accuracy, by=list(Iter=performance.external.2$Iter),mean)$x  # calculate average external accuracy during each iteration
performance.long.2    <- data.frame(Iter = c(performance.2$Iter,performance.2$Iter), Accuracy = c(performance.2$AccuracyExternal,performance.2$Accuracy), Group=c(rep("external",40),rep("internal",40)))

# extract average resampling accuracy at optimal iteration 
accuracy.external.2 <- GA.2$averages[GA.2$optIter,2]
accuracy.external.2

# extract all resampling accuracies at optimal iteration to compute confidence intervalls (below)
accuracy.external.opt.2 <- subset(performance.external.2,performance.external.2$Iter == GA.2$optIter)
accuracy.external.opt.2 <- accuracy.external.opt.2$Accuracy  

## plot internal and external accuracy over the iterations 
ggplot(performance.long.2, aes(Iter, Accuracy, col = Group)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw()


### plot external accuracy over the iterations ################################################################
performance.aggr.2 <- subset(performance.long.2,performance.long.2$Group =="external")

ggplot(performance.aggr.2, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())




##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.2     <- t(matrix.train.2[,GA.2$optVariables])
matrix.test.rfe.2      <- matrix.test.2[GA.2$optVariables,]
matrix.opt.2           <- cbind(matrix.train.rfe.2, matrix.test.rfe.2)
pData.opt.2            <- rbind(pData.train.2,pData.test.2)
pData.opt.2$Set        <- ifelse(row.names(pData.opt.2)%in% row.names(pData.train.2),"train","test")
pData.opt.2$Design_Color <- ifelse(pData.opt.2$Set == "train", pData.opt.2$Design_Color,"#000000")

pca.2         <- prcomp(t(matrix.opt.2))           
plot(pca.2$x, pch=16, col=pData.opt.2$Design_Color, cex=1, asp=1)
legend(4,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.2$Design_Color), pch=16, bty="n", cex=0.8)

# clean up
rm(matrix.train.rfe.2,matrix.test.rfe.2,matrix.opt.2,pData.opt.2)

##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet.2 #######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.test.rfe.2  <- t(matrix.test.2[rfe.2$optVariables,])
matrix.test.GA.2   <- t(matrix.test.2[GA.2$optVariables,])
matrix.train.GA.2  <- matrix.train.2[,GA.2$optVariables]

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

# train SVM and predict TestSet02 on 14 optVars found by rfe-GA ##############################################
set.seed(721)
svmOpt.GA.2  <- train(matrix.train.GA.2,labels.train.2,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.2)
svmOpt.GA.2  

Prediction_GA.2 <- predict(svmOpt.GA.2,matrix.test.GA.2, type = "prob")
Prediction_GA.2$Prediction_GA.2 <- ifelse(Prediction_GA.2$transforming>0.50,"transforming","untransforming")
Prediction_GA.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_GA.2)
write.table(Prediction_GA.2, file = paste("TestSet02_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


##### Predict TestSet02 with all predictors ##################################################################
matrix.test.2.full <- matrix.test.2[row.names(t(matrix.train.2)),]
Prediction_SVM_full.2 <- predict(svmFull.2,t(matrix.test.2.full), type = "prob")
Prediction_SVM_full.2$Prediction_SVM_full <- ifelse(Prediction_SVM_full.2$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_SVM_full.2)
write.table(Prediction_SVM_full.2, file = paste("TestSet02_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 02 #####################################################
#############################################################################################################

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
#############################################################################################################
sink("TestSet02_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.2$Prediction_rfe.2), as.factor(Prediction_rfe.2$TrueLabel))
sink()

sink("TestSet02_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.2$Prediction_GA.2), as.factor(Prediction_GA.2$TrueLabel))
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


Prediction_GA.2$Class <- as.factor(ifelse(Prediction_GA.2$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA <- roc(Prediction_GA.2$Class,                    # response vector (factor or character)
              Prediction_GA.2$transforming,             # predictor vector (numeric)
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.3    # 36 predictors found 
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

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ########################### 
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.4    # 8 predictors chosen
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

#### 6.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.5    # 3 optimal predictors found 
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

#### 6.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.6  # 40 optVars found 
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

### plot external accuracy over the iterations ################################################################
performance.aggr.6 <- subset(performance.long.6,performance.long.6$Group =="external")

ggplot(performance.aggr.6, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

#############################################################################################################
##### 5  PCA on optimal variables ###########################################################################
#############################################################################################################
matrix.train.rfe.6     <- t(matrix.train.6[,GA.6$optVariables])
matrix.test.rfe.6      <- matrix.test.6[GA.6$optVariables,]
matrix.opt.6           <- cbind(matrix.train.rfe.6, matrix.test.rfe.6)
pData.opt.6            <- rbind(pData.train.6,pData.test.6)
pData.opt.6$Set        <- ifelse(row.names(pData.opt.6)%in% row.names(pData.train.6),"train","test")
pData.opt.6$Design_Color <- ifelse(pData.opt.6$Set == "train", pData.opt.6$Design_Color,"#000000")

pca.6         <- prcomp(t(matrix.opt.6))           
plot(pca.6$x, pch=16, col=pData.opt.6$Design_Color, cex=1.5, asp=1)
legend(-5,3.1, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.6$Design_Color), pch=16, bty="n", cex=0.8)

# clean up
rm(matrix.train.rfe.6,matrix.test.rfe.6,matrix.opt.6,pData.opt.6)

##############################################################################################################
##### 6 Train SVM on optVars determined by rfe and GA and predict independent TestSet 06 #####################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.test.rfe.6  <- t(matrix.test.6[rfe.6$optVariables,])
matrix.test.GA.6   <- t(matrix.test.6[GA.6$optVariables,])
matrix.train.GA.6  <- matrix.train.6[,GA.6$optVariables]

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

# train SVM and predict TestSet06 on 14 optVars found by rfe-GA
set.seed(721)
svmOpt.GA.6  <- train(matrix.train.GA.6,labels.train.6,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.6)
svmOpt.GA.6  

Prediction_GA.6 <- predict(svmOpt.GA.6,matrix.test.GA.6, type = "prob")
Prediction_GA.6$Prediction_GA.6 <- ifelse(Prediction_GA.6$transforming>0.50,"transforming","untransforming")
Prediction_GA.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_GA.6)
write.table(Prediction_GA.6, file = paste("TestSet06_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


##### Predict TestSet06 with all predictors
matrix.test.6.full <- matrix.test.6[row.names(t(matrix.train.6)),]
Prediction_SVM_full.6 <- predict(svmFull.6,t(matrix.test.6.full), type = "prob")
Prediction_SVM_full.6$Prediction_SVM_full <- ifelse(Prediction_SVM_full.6$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_SVM_full.6)
write.table(Prediction_SVM_full.6, file = paste("TestSet06_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 06 #####################################################
#############################################################################################################

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
sink("TestSet06_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.6$Prediction_rfe.6), as.factor(Prediction_rfe.6$TrueLabel))
sink()

sink("TestSet06_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.6$Prediction_GA.6), as.factor(Prediction_GA.6$TrueLabel))
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


Prediction_GA.6$Class <- as.factor(ifelse(Prediction_GA.6$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA <- roc(Prediction_GA.6$Class,                    # response vector (factor or character)
              Prediction_GA.6$transforming,             # predictor vector (numeric)
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
#### Split 07 / TestSet 07 - Figure 3 b / c #################################################################################################
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.7    # 26 optVars found by rfe
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
GA.7    # 12 optVars found by GA
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

### Figure 3c plot external accuracy over the iterations ################################################################
performance.aggr.7 <- subset(performance.long.7,performance.long.7$Group =="external")

ggplot(performance.aggr.7, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T,colour = "#E8534F") +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())


##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.7     <- t(matrix.train.7[,GA.7$optVariables])
matrix.test.rfe.7      <- matrix.test.7[GA.7$optVariables,]
matrix.opt.7           <- cbind(matrix.train.rfe.7, matrix.test.rfe.7)
pData.opt.7            <- rbind(pData.train.7,pData.test.7)
pData.opt.7$Set        <- ifelse(row.names(pData.opt.7)%in% row.names(pData.train.7),"train","test")
pData.opt.7$Design_Color <- ifelse(pData.opt.7$Set == "train", pData.opt.7$Design_Color,"#000000")

pca.7         <- prcomp(t(matrix.opt.7))           
plot(pca.7$x, pch=16, col=pData.opt.7$Design_Color, cex=1.5, asp=1)
legend(3,3.7, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.7$Design_Color), pch=16, bty="n", cex=0.8)

# clean up
rm(matrix.train.rfe.7,matrix.test.rfe.7,matrix.opt.7,pData.opt.7)

##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet 07 ######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.test.rfe.7  <- t(matrix.test.7[rfe.7$optVariables,])
matrix.test.GA.7   <- t(matrix.test.7[GA.7$optVariables,])
matrix.train.GA.7  <- matrix.train.7[,GA.7$optVariables]

# train SVM and predict TestSet07 on 26 optVars found by rfe
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

# train SVM and predict TestSet07 on 12 optVars found by rfe-GA
set.seed(721)
svmOpt.GA.7  <- train(matrix.train.GA.7,labels.train.7,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.7)
svmOpt.GA.7  

Prediction_GA.7 <- predict(svmOpt.GA.7,matrix.test.GA.7, type = "prob")
Prediction_GA.7$Prediction_GA.7 <- ifelse(Prediction_GA.7$transforming>0.50,"transforming","untransforming")
Prediction_GA.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_GA.7)
write.table(Prediction_GA.7, file = paste("TestSet07_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


##### Predict TestSet07 with all predictors
matrix.test.7.full <- matrix.test.7[row.names(t(matrix.train.7)),]
Prediction_SVM_full.7 <- predict(svmFull.7,t(matrix.test.7.full), type = "prob")
Prediction_SVM_full.7$Prediction_SVM_full <- ifelse(Prediction_SVM_full.7$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_SVM_full.7)
write.table(Prediction_SVM_full.7, file = paste("TestSet07_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 7 Performance of optVars Classifier on  TestSet 07 #####################################################
#############################################################################################################

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
sink("TestSet07_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.7$Prediction_rfe.7), as.factor(Prediction_rfe.7$TrueLabel))
sink()

sink("TestSet07_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.7$Prediction_GA.7), as.factor(Prediction_GA.7$TrueLabel))
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


Prediction_GA.7$Class <- as.factor(ifelse(Prediction_GA.7$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA <- roc(Prediction_GA.7$Class,                    # response vector (factor or character)
              Prediction_GA.7$transforming,             # predictor vector (numeric)
              percent=TRUE, levels=c("nontransforming","transforming"),
              plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
              print.auc=T,print.thres=0.5)

Prediction_SVM_full.7$Class <- as.factor(ifelse(Prediction_SVM_full.7$TrueLabel == "transforming","transforming","nontransforming"))
roc.full<- roc(Prediction_SVM_full.7$Class,                    # response vector (factor or character)
               Prediction_SVM_full.7$transforming,             # predictor vector (numeric)
               percent=TRUE, levels=c("nontransforming","transforming"),
               plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
               print.auc=T,print.thres=0.5)


### ROC of 2140 cross-validation predictions 
pred.7 <- rfe.7$pred
pred.7 <- subset(pred.7, pred.7$Variables == 26)  # 2140 predictions with the best feature subset (26 variables)

pred.7$Class <- as.factor(ifelse(pred.7$obs == "transforming","transforming","nontransforming"))
roc.pred.7   <- roc(pred.7$Class,                    
                    pred.7$transforming,             
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.8   # 20 variables found by SVM-rfe
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

### plot external accuracy over the iterations ################################################################
performance.aggr.8 <- subset(performance.long.8,performance.long.8$Group =="external")

ggplot(performance.aggr.8, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())



##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.8     <- t(matrix.train.8[,GA.8$optVariables])
matrix.test.rfe.8      <- matrix.test.8[GA.8$optVariables,]
matrix.opt.8           <- cbind(matrix.train.rfe.8, matrix.test.rfe.8)
pData.opt.8            <- rbind(pData.train.8,pData.test.8)
pData.opt.8$Set        <- ifelse(row.names(pData.opt.8)%in% row.names(pData.train.8),"train","test")
pData.opt.8$Design_Color <- ifelse(pData.opt.8$Set == "train", pData.opt.8$Design_Color,"#000000")

pca.8         <- prcomp(t(matrix.opt.8))           
plot(pca.8$x, pch=16, col=pData.opt.8$Design_Color, cex=1.5, asp=1)
legend(2,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.8$Design_Color), pch=16, bty="n", cex=0.8)

# clean up
rm(matrix.train.rfe.8,matrix.test.rfe.8,matrix.opt.8,pData.opt.8)

##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet 08 ######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.test.rfe.8  <- t(matrix.test.8[rfe.8$optVariables,])
matrix.test.GA.8   <- t(matrix.test.8[GA.8$optVariables,])
matrix.train.GA.8  <- matrix.train.8[,GA.8$optVariables]

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

# train SVM and predict TestSet08 on 10 optVars found by rfe-GA
set.seed(721)
svmOpt.GA.8  <- train(matrix.train.GA.8,labels.train.8,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.8)
svmOpt.GA.8  

Prediction_GA.8 <- predict(svmOpt.GA.8,matrix.test.GA.8, type = "prob")
Prediction_GA.8$Prediction_GA.8 <- ifelse(Prediction_GA.8$transforming>0.50,"transforming","untransforming")
Prediction_GA.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_GA.8)
write.table(Prediction_GA.8, file = paste("TestSet08_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


##### Predict TestSet08 with all predictors
matrix.test.8.full <- matrix.test.8[row.names(t(matrix.train.8)),]
Prediction_SVM_full.8 <- predict(svmFull.8,t(matrix.test.8.full), type = "prob")
Prediction_SVM_full.8$Prediction_SVM_full <- ifelse(Prediction_SVM_full.8$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_SVM_full.8)
write.table(Prediction_SVM_full.8, file = paste("TestSet08_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 7 Performance of optVars Classifier on  TestSet 08 #####################################################
#############################################################################################################

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
sink("TestSet08_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.8$Prediction_rfe.8), as.factor(Prediction_rfe.8$TrueLabel))
sink()

sink("TestSet08_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.8$Prediction_GA.8), as.factor(Prediction_GA.8$TrueLabel))
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

Prediction_GA.8$Class <- as.factor(ifelse(Prediction_GA.8$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA <- roc(Prediction_GA.8$Class,                    # response vector (factor or character)
              Prediction_GA.8$transforming,             # predictor vector (numeric)
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.9 # 7 optVars found by rfe
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

#### 6.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
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
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.10   # 45 optVars selected
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

### plot external accuracy over the iterations ################################################################
performance.aggr.10 <- subset(performance.long.10,performance.long.10$Group =="external")

ggplot(performance.aggr.10, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.5,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())


##############################################################################################################
##### 5  PCA on optimal variables ############################################################################
##############################################################################################################
matrix.train.rfe.10     <- t(matrix.train.10[,GA.10$optVariables])
matrix.test.rfe.10      <- matrix.test.10[GA.10$optVariables,]
matrix.opt.10           <- cbind(matrix.train.rfe.10, matrix.test.rfe.10)
pData.opt.10            <- rbind(pData.train.10,pData.test.10)
pData.opt.10$Set        <- ifelse(row.names(pData.opt.10)%in% row.names(pData.train.10),"train","test")
pData.opt.10$Design_Color <- ifelse(pData.opt.10$Set == "train", pData.opt.10$Design_Color,"#000000")

pca.10         <- prcomp(t(matrix.opt.10))           
plot(pca.10$x, pch=16, col=pData.opt.10$Design_Color, cex=1.5, asp=1)
legend(3.8,3, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.10$Design_Color), pch=16, bty="n", cex=0.8)

# clean up
rm(matrix.train.rfe.10,matrix.test.rfe.10,matrix.opt.10,pData.opt.10)

##############################################################################################################
##### 6 Train SVM on optVar determined by rfe and GA and predict independent TestSet 10 ######################
##############################################################################################################

# subset matrices on the different optimal predictors 
matrix.test.rfe.10  <- t(matrix.test.10[rfe.10$optVariables,])
matrix.test.GA.10   <- t(matrix.test.10[GA.10$optVariables,])
matrix.train.GA.10  <- matrix.train.10[,GA.10$optVariables]

# train SVM and predict TestSet10 on 45 optVars found by rfe
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

# train SVM and predict TestSet10 on 19 optVars found by rfe-GA
set.seed(721)
svmOpt.GA.10  <- train(matrix.train.GA.10,labels.train.10,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 20,
                      trControl = fullCtrl.10)
svmOpt.GA.10  

Prediction_GA.10 <- predict(svmOpt.GA.10,matrix.test.GA.10, type = "prob")
Prediction_GA.10$Prediction_GA.10 <- ifelse(Prediction_GA.10$transforming>0.50,"transforming","untransforming")
Prediction_GA.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_GA.10)
write.table(Prediction_GA.10, file = paste("TestSet10_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet10 with all predictors
matrix.test.10.full <- matrix.test.10[row.names(t(matrix.train.10)),]
Prediction_SVM_full.10 <- predict(svmFull.10,t(matrix.test.10.full), type = "prob")
Prediction_SVM_full.10$Prediction_SVM_full <- ifelse(Prediction_SVM_full.10$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_SVM_full.10)
write.table(Prediction_SVM_full.10, file = paste("TestSet10_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#############################################################################################################
#### 7 Performance of optVars Classifier on  TestSet 10 #####################################################
#############################################################################################################

#### 7.1 Confusion matrix /Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8  ###########################  
sink("TestSet10_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.10$Prediction_rfe.10), as.factor(Prediction_rfe.10$TrueLabel))
sink()

sink("TestSet10_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.10$Prediction_GA.10), as.factor(Prediction_GA.10$TrueLabel))
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


Prediction_GA.10$Class <- as.factor(ifelse(Prediction_GA.10$TrueLabel == "transforming","transforming","nontransforming"))
roc.GA <- roc(Prediction_GA.10$Class,                    # response vector (factor or character)
              Prediction_GA.10$transforming,             # predictor vector (numeric)
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
summary(fselect.FINAL)                              # 1243 features selected by nonspecific filtering

Annotation.matrix.train.FINAL <- Annotation[row.names(matrix.train.FINAL),]
write.table(Annotation.matrix.train.FINAL, file = "Annotation.matrix.train.FINAL.txt", sep="\t",col.names=NA)

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

## export Variable importances from the final set / Supplementary Table 4/ Variable importance
varImp_FINAL <- varImp(svmFull.FINAL)$importance
varImp_FINAL <- varImp_FINAL[order(varImp_FINAL$transforming,decreasing = T),]
varImp_FINAL <- cbind(varImp_FINAL,Annotation[row.names(varImp_FINAL),])
write.table(varImp_FINAL, file = "FINAL_VariableImportance_ROC_IQR0.8.txt", sep="\t",col.names=NA)

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
                              method="svmRadial",
                              tuneLength = 20,
                              trControl = innerctrl))

rfe.FINAL   # 20 optVars found by rfe
write.table(rfe.FINAL$results, file = "FinalSet152_Results_rfe.txt", sep="\t",col.names=NA)

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
GA.FINAL  # 11 Features selected by GA
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

### Figure 3j: PCA on 36,226 probes
pca.eset.batch <- prcomp(t(eset.batch))           
plot(pca.eset.batch$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(35,55, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

### Figure 3k: PCA on 11 probes from SVM-GA
matrix.train.GA.FINAL  <- eset.batch[GA.FINAL$optVariables,]
pca.FINAL              <- prcomp(t(matrix.train.GA.FINAL))           
plot(pca.FINAL$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(2,4, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

### PCA on 20 probes from SVM-RFE
matrix.train.rfe.FINAL <- eset.batch[rfe.FINAL$optVariables,]
pca.rfe.FINAL          <- prcomp(t(matrix.train.rfe.FINAL))           
plot(pca.rfe.FINAL$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(2,5, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

################################################################################################################
### 6. ROC curve for outer resamples of SVM.FINAL / SVM.rfe.FINAL: #############################################
################################################################################################################

### 6.1 ROC curve SVM.rfe resamples 
a <- rfe.FINAL$pred
b <- subset(a, a$Variables == 20)  # only for the best feature subset

b$Class <- as.factor(ifelse(b$obs == "transforming","transforming","nontransforming"))
roc.b <- roc(b$Class,                    
             b$transforming,             
             percent=TRUE, levels=c("nontransforming","transforming"),
             plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
             print.auc=T,print.thres=0.5)





#############################################################################################################################################
#############################################################################################################################################
#### Analyze all results from the different Test/Training splits ############################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### 1. proportion of samples that were predicted at least once:  ############################################################################
#############################################################################################################################################
a <- unique(c(row.names(pData.test.1),row.names(pData.test.2),row.names(pData.test.3),row.names(pData.test.4),row.names(pData.test.5),
              row.names(pData.test.6),row.names(pData.test.7),row.names(pData.test.8),row.names(pData.test.9),row.names(pData.test.10)))
length(a)  # 141/152 = 93 % of the samples have been predicted at least once

#############################################################################################################################################
### 2. Supplementary Table 4: summary of features, resampling- and test set accuracies ######################################################
#############################################################################################################################################
library(gmodels)

### 1. create data.frame 
Results <- data.frame(Testset = c(seq(1,10,by=1),"FinalSet152"), FeatureNumber_full = rep(NA,11), ResamplingAccuracy_full=rep(NA,11),CI_full_lower=rep(NA,11),CI_full_upper=rep(NA,11), 
                      FeatureNumber_rfe = rep(NA,11), ResamplingAccuracy_rfe = rep(NA,11), CI_rfe_lower=rep(NA,11),CI_rfe_upper=rep(NA,11), P_full_vs_rfe=rep(NA,11),
                      FeatureNumber_GA  = rep(NA,11), ResamplingAccuracy_GA = rep(NA,11),CI_GA_lower=rep(NA,11),CI_GA_upper=rep(NA,11),
                      Accuracy_TestSet_full = rep(NA,11), CI_Test_full_lower=rep(NA,11),CI_Test_full_upper=rep(NA,11),
                      Accuracy_TestSet_rfe = rep(NA,11),  CI_Test_rfe_lower=rep(NA,11),CI_Test_rfe_upper=rep(NA,11),
                      Accuracy_TestSet_GA = rep(NA,11),   CI_Test_GA_lower=rep(NA,11),CI_Test_GA_upper=rep(NA,11))
                      
### 2. number of features selected by IQR filter for each training/test split
for (i in 1:10) {Results[i,2] <- dim(eval(parse(text=paste("matrix.train.",i, sep=""))))[2] }
Results[11,2] <- dim(matrix.train.FINAL)[2]

### 3. resampling accuracy for the full model / stored in svmFull.x$results dataframe
for (i in 1:10) { Results[i,3] <- eval(parse(text=paste("svmFull.",i, sep="")))$results[eval(parse(text=paste("svmFull.",i, sep="")))$results$C == eval(parse(text=paste("svmFull.",i, sep="")))$bestTune$C,"Accuracy"] }
Results[11,3] <- svmFull.FINAL$results[svmFull.FINAL$results$C == svmFull.FINAL$bestTune$C,"Accuracy"]

### 4. lower bound of confidence intervall for the resampling accuracy of full model / stored in rfeResamples.x$values dataframe
for (i in 1:10) {Results[i,4] <- ci(eval(parse(text=paste("rfeResamples.",i, sep="")))$values[,2])[2]}
Results[11,4] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[2] 

### 5. upper bound of confidence intervall for the resampling accuracy of full model / stored in rfeResamples.x$values dataframe
for (i in 1:10) {Results[i,5] <- ci(eval(parse(text=paste("rfeResamples.",i, sep="")))$values[,2])[3]}
Results[11,5] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[3] 

### 6. number of features found by SVM-rfe for each training/test split
for (i in 1:10) {Results[i,6] <- eval(parse(text=paste("rfe.",i, sep="")))$optsize }
Results[11,6] <- rfe.FINAL$optsize

### 7. external resampling accuracy for the svm-rfe model
for (i in 1:10) { Results[i,7] <- eval(parse(text=paste("rfe.",i, sep="")))$results[eval(parse(text=paste("rfe.",i, sep="")))$results$Variables == eval(parse(text=paste("rfe.",i, sep="")))$optsize,5] }
Results[11,7] <- rfe.FINAL$results[rfe.FINAL$results$Variables == rfe.FINAL$optsize,5]

### 8. fill in lower bound confidence intervall for the resampling accuracy of rfe model / stored in rfeResamples.x$values dataframe
for (i in 1:10) {Results[i,8] <- ci(eval(parse(text=paste("rfeResamples.",i, sep="")))$values[,7])[2]}
Results[11,8] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[2] 

### 9. upper bound confidence intervall for the resampling accuracy of rfe model / stored in rfeResamples.x$values dataframe
for (i in 1:10) {Results[i,9] <- ci(eval(parse(text=paste("rfeResamples.",i, sep="")))$values[,7])[3]}
Results[11,9] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[3] 

### 10. p.value for H0 = resampling accuracy of full model = resampling accuracy of full model/ stored in rfeResamples.x$values dataframe
for (i in 1:10) {Results[i,10] <- eval(parse(text=paste("modelDifferences.",i, sep="")))$statistics$Accuracy[[1]][3]}
Results[11,10] <- modelDifferences.FINAL$statistics$Accuracy[[1]][3]

### 11. number of features found by SVM-GA 
for (i in c(2,3,6,7,8,10)) { Results[i,11] <- length(eval(parse(text=paste("GA.",i, sep="")))$optVariables)  }
Results[11,11] <- length(GA.FINAL$optVariables)

### 12. Resampling accuracy for the GA model
for (i in c(2,3,6,7,8,10)) { Results[i,12] <- eval(parse(text=paste("accuracy.external.",i, sep=""))) }
Results[11,12] <- accuracy.external.FINAL

### 13. calculate lower bound confidence intervall for the resampling accuracy of GA model / stored in rfeResamples.x$values dataframe
for (i in c(2,3,6,7,8,10)) {Results[i,13] <- ci(eval(parse(text=paste("accuracy.external.opt.",i, sep=""))))[2]}
Results[11,13] <- ci(accuracy.external.opt.FINAL)[2] 

### 14. calculate upper bound confidence intervall for the resampling accuracy of GA model / stored in rfeResamples.x$values dataframe
for (i in c(2,3,6,7,8,10)) {Results[i,14] <- ci(eval(parse(text=paste("accuracy.external.opt.",i, sep=""))))[3]}
Results[11,14] <- ci(accuracy.external.opt.FINAL)[3] 

### 15. TestSet accuracy for the full model
for (i in 1:10) {Results[i,15] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$Prediction_SVM_full), as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$TrueLabel))$overall[1] }

### 16. lower bound of confidence intervall of the TestSet accuracy for the full model
for (i in 1:10) {Results[i,16] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$Prediction_SVM_full), as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$TrueLabel))$overall[3] }

### 17. upper bound confidence intervall of the TestSet accuracy for the full model
for (i in 1:10) {Results[i,17] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$Prediction_SVM_full), as.factor(eval(parse(text=paste("Prediction_SVM_full.",i, sep="")))$TrueLabel))$overall[4] }

### 18. TestSet accuracy for the SVM-rfe model
for (i in 1:10) { Results[i,18] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))$TrueLabel))$overall[1] }

### 19. lower bound of confidence intervall of the TestSet accuracy for the SVM-rfe model
for (i in 1:10) { Results[i,19] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))$TrueLabel))$overall[3] }

### 20. upper bound confidence intervall of the TestSet accuracy for the SVM-rfe model
for (i in 1:10) { Results[i,20] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_rfe.",i, sep="")))$TrueLabel))$overall[4] }

### 21. TestSet accuracy for the SVM-GA model
for (i in c(2,3,6,7,8,10)) { Results[i,21] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))$TrueLabel))$overall[1] }

### 22. lower bound of confidence intervall TestSet accuracy for the SVM-GA model
for (i in c(2,3,6,7,8,10)) { Results[i,22] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))$TrueLabel))$overall[3] }

### 23. upper bound of confidence intervall TestSet accuracy for the SVM-GA model
for (i in c(2,3,6,7,8,10)) { Results[i,23] <- confusionMatrix(as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))[,7]), as.factor(eval(parse(text=paste("Prediction_GA.",i, sep="")))$TrueLabel))$overall[4] }

write.table(Results, file = "AllResults_rfe_GA40_10TestSets_FinalSet152_IQR0.8.txt", sep="\t",col.names=NA)

#############################################################################################################################################
### 3. Plot resampling results vs TestSet results: Figure 3 d-e #############################################################################
#############################################################################################################################################

df.results <- Results[c(1:10),]   # export data for the 10 Training / Test Sets w/o FinalSet

### for full model / Figure 3d: #############################################################################
Accuracy_TestSet_full.median   <- median(df.results$Accuracy_TestSet_full)
ResamplingAccuracy_full.median <- median(df.results$ResamplingAccuracy_full)

ggplot(df.results, aes(Accuracy_TestSet_full,ResamplingAccuracy_full, label = Testset)) +
  geom_point(size=7) +
  geom_errorbarh(aes(xmin = CI_Test_full_lower,  xmax= CI_Test_full_upper, height = 0.003),colour = "black",size = 0.2) +
  geom_errorbar (aes(ymin = CI_full_lower, ymax = CI_full_upper, width = 0.003),colour = "black",size = 0.3) +
  geom_point(aes(x=Accuracy_TestSet_full.median, y=ResamplingAccuracy_full.median), colour="red", size = 7)+
  coord_fixed(ratio = 1, xlim=c(0.65,1),ylim=c(0.65,1)) +
  geom_abline(intercept = 0, slope = 1, color = "red",linetype = "dashed") + 
  geom_text(aes(label=Testset),hjust=0.5, vjust=0.5, colour = "white", size = 5,fontface="bold") +
  scale_y_continuous(breaks=seq(0.65,1, 0.05)) + 
  scale_x_continuous(breaks=seq(0.65,1, 0.05)) + 
  xlab("Test set accuracy" ) +
  ylab("Cross-validation accuracy") +
  theme_bw() +
  theme(axis.title.x = element_text(size=17),axis.title.y = element_text(size=17),
        axis.text = element_text(size=14, color ="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())


### for SVM-rfe / Figure 3e:  #############################################################################
Accuracy_TestSet_rfe.median   <- median(df.results$Accuracy_TestSet_rfe)
ResamplingAccuracy_rfe.median <- median(df.results$ResamplingAccuracy_rfe)

ggplot(df.results, aes(Accuracy_TestSet_rfe,ResamplingAccuracy_rfe, label = Testset)) +
  geom_point(size=7) +
  geom_errorbarh(aes(xmin = CI_Test_rfe_lower,  xmax= CI_Test_rfe_upper, height = 0.003),colour = "black",size = 0.2) +
  geom_errorbar (aes(ymin = CI_rfe_lower, ymax = CI_rfe_upper, width = 0.003),colour = "black",size = 0.3) +
  geom_point(aes(x=Accuracy_TestSet_rfe.median, y=ResamplingAccuracy_rfe.median), colour="red", size = 7)+
  coord_fixed(ratio = 1, xlim=c(0.65,1),ylim=c(0.65,1)) +
  geom_abline(intercept = 0, slope = 1, color = "red",linetype = "dashed") + 
  geom_text(aes(label=Testset),hjust=0.5, vjust=0.5, colour = "white", size = 5,fontface="bold") +
  scale_y_continuous(breaks=seq(0.65,1, 0.05)) + 
  scale_x_continuous(breaks=seq(0.65,1, 0.05)) + 
  xlab("Test set accuracy" ) +
  ylab("Cross-validation accuracy") +
  theme_bw() +
  theme(axis.title.x = element_text(size=17),axis.title.y = element_text(size=17),
        axis.text = element_text(size=14, color ="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())


### for SVM-GA / Figure 3f  #############################################################################
df.results.GA              <- df.results[c(2,3,6,7,8,10),]    # subset for the Training / TestSet splits in which GA was performed 
Accuracy_TestSet_GA.median   <- median(df.results.GA$Accuracy_TestSet_GA)
ResamplingAccuracy_GA.median <- median(df.results.GA$ResamplingAccuracy_GA)

ggplot(df.results.GA, aes(Accuracy_TestSet_GA,ResamplingAccuracy_GA, label = Testset)) +
  geom_point(size=7) +
  geom_errorbarh(aes(xmin = CI_Test_GA_lower,  xmax= CI_Test_GA_upper, height = 0.003),colour = "black",size = 0.2) +
  geom_errorbar (aes(ymin = CI_GA_lower, ymax = CI_GA_upper, width = 0.003),colour = "black",size = 0.3) +
  geom_point(aes(x=Accuracy_TestSet_GA.median, y=ResamplingAccuracy_GA.median), colour="red", size = 7)+
  coord_fixed(ratio = 1, xlim=c(0.65,1),ylim=c(0.65,1)) +
  geom_abline(intercept = 0, slope = 1, color = "red",linetype = "dashed") + 
  geom_text(aes(label=Testset),hjust=0.5, vjust=0.5, colour = "white", size = 5,fontface="bold") +
  scale_y_continuous(breaks=seq(0.65,1, 0.05)) + 
  scale_x_continuous(breaks=seq(0.65,1, 0.05)) + 
  xlab("Test set accuracy" ) +
  ylab("Cross-validation accuracy") +
  theme_bw() +
  theme(axis.title.x = element_text(size=17),axis.title.y = element_text(size=17),
        axis.text = element_text(size=14, color ="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())


 
#############################################################################################################################################
### 4. aggregate Test set performances over the different splits  ###########################################################################
#############################################################################################################################################

### 4.1 TestSet performances for all 10 SVM_full models   ###################################################################################
Predictions_SVM_full    <- rbind(Prediction_SVM_full.1,Prediction_SVM_full.2, Prediction_SVM_full.3,Prediction_SVM_full.4,Prediction_SVM_full.5,
                                 Prediction_SVM_full.6,Prediction_SVM_full.7,Prediction_SVM_full.8,Prediction_SVM_full.9,Prediction_SVM_full.10)

### 4.2.TestSet performances for the 6 SVM_full models for which rfe and GA was performed ###################################################
Predictions_SVM_full_GA <- rbind(Prediction_SVM_full.2, Prediction_SVM_full.3,
                                 Prediction_SVM_full.6,Prediction_SVM_full.7,Prediction_SVM_full.8,Prediction_SVM_full.10)

### 4.3 TestSet performances for all 10 SVM_rfe models   ####################################################################################
colnames(Prediction_rfe.1)[7] <- "Prediction" 
colnames(Prediction_rfe.2)[7] <- "Prediction" 
colnames(Prediction_rfe.3)[7] <- "Prediction" 
colnames(Prediction_rfe.4)[7] <- "Prediction" 
colnames(Prediction_rfe.5)[7] <- "Prediction" 
colnames(Prediction_rfe.6)[7] <- "Prediction" 
colnames(Prediction_rfe.7)[7] <- "Prediction" 
colnames(Prediction_rfe.8)[7] <- "Prediction" 
colnames(Prediction_rfe.9)[7] <- "Prediction" 
colnames(Prediction_rfe.10)[7] <- "Prediction" 

Predictions_RFE   <- rbind(Prediction_rfe.1,Prediction_rfe.2, Prediction_rfe.3,Prediction_rfe.4,Prediction_rfe.5,
                           Prediction_rfe.6,Prediction_rfe.7,Prediction_rfe.8,Prediction_rfe.9,Prediction_rfe.10)

### 4.4.TestSet performances for the 6 SVM_rfe models for which GA was performed ###########################################################
Predictions_RFE_GA   <- rbind(Prediction_rfe.2, Prediction_rfe.3,
                              Prediction_rfe.6,Prediction_rfe.7,Prediction_rfe.8,Prediction_rfe.10)

### 4.5 TestSet performances for all 6 SVM_GA models   #####################################################################################
colnames(Prediction_GA.2)[7] <- "Prediction" 
colnames(Prediction_GA.3)[7] <- "Prediction" 
colnames(Prediction_GA.6)[7] <- "Prediction" 
colnames(Prediction_GA.7)[7] <- "Prediction" 
colnames(Prediction_GA.8)[7] <- "Prediction" 
colnames(Prediction_GA.10)[7] <- "Prediction" 
Predictions_GA <- rbind(Prediction_GA.2, Prediction_GA.3,Prediction_GA.6,Prediction_GA.7,Prediction_GA.8,Prediction_GA.10)

#############################################################################################################################################
### 4.6 SAGA = Compound TestSet performances = 4xSVM-rfe (RFE < 10 predictors) + 6xSVM_GA models (RFE >10 predictors) #######################
#############################################################################################################################################
Predictions_SAGA <- rbind(Prediction_rfe.1,Prediction_GA.2, Prediction_GA.3,Prediction_rfe.4,Prediction_rfe.5,
                          Prediction_GA.6,Prediction_GA.7,Prediction_GA.8,Prediction_rfe.9,Prediction_GA.10)

#############################################################################################################################################
### 5. Confusion matrices of aggregated TestSets / Supplementary Table 4/ 4b_Performance_TestSets_IQR0.8   ##################################
#############################################################################################################################################

sink("AllTestSetsConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_SVM_full$Prediction_SVM_full), as.factor(Predictions_SVM_full$TrueLabel))
sink()

sink("AllTestSetsConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_RFE$Prediction), as.factor(Predictions_RFE$TrueLabel))
sink()

sink("AllTestSets_ConfusionMatrix_SAGA.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_SAGA$Prediction), as.factor(Predictions_SAGA$TrueLabel))
sink()

sink("TestSets2367810_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_GA$Prediction), as.factor(Predictions_GA$TrueLabel))
sink()

sink("TestSets2367810_ConfusionMatrix_Full.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_SVM_full_GA$Prediction), as.factor(Predictions_SVM_full_GA$TrueLabel))
sink()

sink("TestSets2367810_ConfusionMatrix_RFE.txt", append = TRUE)
confusionMatrix(as.factor(Predictions_RFE_GA$Prediction), as.factor(Predictions_RFEGA$TrueLabel))
sink()


#############################################################################################################################################
### 6. ROC SAGA vs IVIM (all available IVIMs ): #############################################################################################
#############################################################################################################################################

### 6.1 read in all available IVIM data (from Fig. 1d, n=502) ###############################################################################
#############################################################################################################################################
IVIM  <- read.delim("IVIM_MTT_ROC.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE)
IVIM$Classification_Code <- as.factor(IVIM$Classification_Code)
IVIM$Prediction <- ifelse(IVIM$MTT_Score>=3,"transforming","untransforming")  
IVIM$TrueLabel  <- ifelse(IVIM$Classification_Code==1,"transforming","untransforming")
confusionMatrix(as.factor(IVIM$Prediction), as.factor(IVIM$TrueLabel))

### 6.2.Figure 3g / S4a: IVIM vs SAGA performance (for all vectors) #########################################################################
#############################################################################################################################################
IVIM.total <- roc(IVIM$Classification_Code,    # response vector (factor or character)
                  IVIM$MTT_Score,          # predictor vector (numeric)
                  percent=TRUE, smooth = F,
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                  col = "#8285BC", grid=F, lwd = 3, cex.lab=1.5, 
                  print.auc=T,print.thres=3 )

roc.rfe.all <- roc(Predictions_RFE$Class,                    
                   Predictions_RFE$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col = "#4A6893", grid=F,
                   print.auc=T,print.thres=0.5, add = T)

roc.GA.all <- roc( Predictions_GA$Class,                    
                   Predictions_GA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#E8534F", grid=F,
                   print.auc=T, add = T,print.thres=0.5)

roc.SAGA   <- roc( Predictions_SAGA$Class,                    
                   Predictions_SAGA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                   print.auc=T, add = T,print.thres=0.5)

### 6.3. compare ROC curves between SAGA and IVIM ###########################################################################################
#############################################################################################################################################
roc.test(roc.SAGA, IVIM.total, alternative = "greater")   # p-value = 9.173e-09


### 6.4.Supplementary Figure 4b IVIM vs SAGA on LTR.RV.SF only ##############################################################################
#############################################################################################################################################
selected  <- c("MOCK","LTR.RV.SF")  # select Mock controls and LTR.SF in IVIM data
IVIM.sel1 <- subset(IVIM,IVIM$Vector %in% selected)
confusionMatrix(as.factor(IVIM.sel1$Prediction), as.factor(IVIM.sel1$TrueLabel))

IVIM.LTR.SF <- roc(IVIM.sel1$Classification_Code,    # response vector (factor or character)
                   IVIM.sel1$MTT_Score,          # predictor vector (numeric)
                   percent=TRUE, smooth = F,
                   plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                   col = "#8285BC", grid=F, lwd = 3, cex.lab=1.5, 
                   print.auc=T,print.thres=3)

# select Mock controls and LTR.SF in SAGA Predictions in the TestSets
Predictions_SAGA_LTR.SF <- Predictions_SAGA[(grepl("Mock",Predictions_SAGA$Name) | grepl("RV.SF",Predictions_SAGA$Name)),]
confusionMatrix(as.factor(Predictions_SAGA_LTR.SF$Prediction), as.factor(Predictions_SAGA_LTR.SF$TrueLabel))

ROC.SAGA.LTR.SF   <- roc(Predictions_SAGA_LTR.SF$Class,                    
                         Predictions_SAGA_LTR.SF$transforming,             
                         percent=TRUE, levels=c("nontransforming","transforming"),
                         plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                         print.auc=T, add = T,print.thres=0.5)

roc.test(ROC.SAGA.LTR.SF, IVIM.LTR.SF, alternative = "greater" )   #p-value = 5.696e-07

### 6.5.Supplementary Figure 4c IVIM vs SAGA non-LTR vectors ################################################################################
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

roc.test(ROC.SAGA.other, ROC.IVIM.other, alternative = "greater" )

#############################################################################################################################################
### 7. ROC SAGA vs IVIM (only the 152 samples for which SAGA and IVIM are available ): ######################################################
#############################################################################################################################################
pData.152           <- pData
pData.152$TrueLabel <- pData$Class  
pData.152$Class     <- ifelse(pData.152$TrueLabel=="transforming","transforming","nontransforming")

pData.152$IVIM_Prediction <- ifelse(pData.152$MTT_4>=3,"transforming","untransforming")
confusionMatrix(as.factor(pData.152$IVIM_Prediction), as.factor(pData.152$TrueLabel))

### 7.1.Supplementary Figure 4d: IVIM vs SAGA performance (for all vectors) #################################################################
#############################################################################################################################################
IVIM.152   <- roc(pData.152$Class,    
                  pData.152$MTT_4,          
                  percent=TRUE, smooth = F,
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                  col = "#8285BC", grid=F, lwd = 3, cex.lab=1.5, 
                  print.auc=T,print.thres=3 )

roc.SAGA   <- roc( Predictions_SAGA$Class,                    
                   Predictions_SAGA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F7C35E", grid=F,
                   print.auc=T, add = T,print.thres=0.5)

roc.test(roc.SAGA, IVIM.152, alternative = "greater")


### 7.2.Supplementary Figure 4e IVIM vs SAGA on LTR.RV.SF only ##############################################################################
#############################################################################################################################################
selected.2  <- c("A2_MOCK","A1_LTR.RV.SF.eGFP")  
pData.LTR   <- subset(pData.152,pData.152$Vector %in% selected.2)
confusionMatrix(as.factor(pData.LTR$IVIM_Prediction), as.factor(pData.LTR$TrueLabel))

IVIM.LTR.152 <- roc(pData.LTR$Class,    
                    pData.LTR$MTT_4,          
                    percent=TRUE, smooth = F,
                    plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                    col = "#C1C2E0", grid=F, lwd = 3, cex.lab=1.5, 
                    print.auc=T,print.thres=3 )

ROC.SAGA.LTR.SF   <- roc(Predictions_SAGA_LTR.SF$Class,                    
                         Predictions_SAGA_LTR.SF$transforming,             
                         percent=TRUE, levels=c("nontransforming","transforming"),
                         plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                         print.auc=T, add = T,print.thres=0.5)


roc.test(ROC.SAGA.LTR.SF, IVIM.LTR.152,alternative = "greater")

### 7.2.Supplementary Figure 4f IVIM vs SAGA other vectors ##################################################################################
#############################################################################################################################################
pData.noLTR <- subset(pData.152,pData.152$Vector!= "A1_LTR.RV.SF.eGFP")
confusionMatrix(as.factor(pData.noLTR$IVIM_Prediction), as.factor(pData.noLTR$TrueLabel))

IVIM.noLTR.152 <- roc(pData.noLTR$Class,    
                      pData.noLTR$MTT_4,          
                      percent=TRUE, smooth = F,
                      plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                      col = "#C1C2E0", grid=F, lwd = 3, cex.lab=1.5, 
                      print.auc=T,print.thres=3 )

ROC.SAGA.other <- roc(Predictions_SAGA_other$Class,                    
                      Predictions_SAGA_other$transforming,             
                      percent=TRUE, levels=c("nontransforming","transforming"),
                      plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                      print.auc=T, add = T,print.thres=0.5)


roc.test(ROC.SAGA.other, IVIM.noLTR.152,alternative = "greater")





### 6.2. Supplementary Figure 3l: 
roc.full <-  roc(Predictions_SVM_full$Class,                    
                    Predictions_SVM_full$transforming,             
                    percent=TRUE, levels=c("nontransforming","transforming"),
                    plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                    print.auc=T,print.thres=0.5)


roc.rfe.all <- roc(Predictions_RFE$Class,                    
                   Predictions_RFE$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col = "#4A6893", grid=F,
                   print.auc=T,print.thres=0.5, add = T)

roc.GA.all <- roc( Predictions_GA$Class,                    
                   Predictions_GA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#E8534F", grid=F,
                   print.auc=T, add = T,print.thres=0.5)

roc.comp   <- roc( Predictions_comp$Class,                    
                   Predictions_comp$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#F9C35F", grid=F,
                   print.auc=T, add = T,print.thres=0.5)


### 6.2. Supplementary Figure 3m: 
## comparisons for the Test/Training Splits 2,3,6,7,8,10 for which full, rfe and GA is available = 255 predictions each 

roc.full.GA <-  roc(Predictions_SVM_full_GA$Class,                    
                    Predictions_SVM_full_GA$transforming,             
                    percent=TRUE, levels=c("nontransforming","transforming"),
                    plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                    print.auc=T,print.thres=0.5)

roc.rfe.GA <- roc(Predictions_RFE_GA$Class,                    
                  Predictions_RFE_GA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col = "#A50F15", grid=F,
                   print.auc=T,print.thres=0.5, add = T)

roc.GA.all <- roc( Predictions_GA$Class,                    
                   Predictions_GA$transforming,             
                   percent=TRUE, levels=c("nontransforming","transforming"),
                   plot=T, auc.polygon=F, max.auc.polygon=F, col ="#08306B", grid=F,
                   print.auc=T, add = T,print.thres=0.5)

roc.test(roc.full.GA, roc.GA.all, reuse.auc=FALSE, partial.auc=c(70,100), partial.auc.focus="se")
roc.test(roc.full.GA, roc.rfe.GA, reuse.auc=FALSE, partial.auc=c(70,100), partial.auc.focus="se")


#############################################################################################################################################
### 8. compare resampling performances between full and rfe models ##########################################################################
#############################################################################################################################################

# Full models vs RFE-models
mean(df.results$ResamplingAccuracy_full)
mean(df.results$ResamplingAccuracy_rfe)
mean(df.results$P_full_vs_rfe)
t.resample <- t.test(df.results$ResamplingAccuracy_full, df.results$ResamplingAccuracy_rfe, paired = TRUE)
t.resample   

# Full models vs RFE-models
t.TestSet <- t.test(df.results$Accuracy_TestSet_full, df.results$Accuracy_TestSet_rfe, paired = TRUE)
t.TestSet    

# Resampling Full models vs GA-models
mean(df.results.GA$ResamplingAccuracy_full)
mean(df.results.GA$ResamplingAccuracy_GA)
t.resample.GA <- t.test(df.results.GA$ResamplingAccuracy_full, df.results.GA$ResamplingAccuracy_GA, paired = TRUE)
t.resample.GA

# Resampling Full models vs RFE-models for Split 1/4/5/9
df.results.rfe <- df.results[c(1,4,5,9),]
median(df.results.rfe$ResamplingAccuracy_full)
median(df.results.rfe$ResamplingAccuracy_rfe)
median(df.results.rfe$P_full_vs_rfe)
median(df.results.rfe$Accuracy_TestSet_full)
median(df.results.rfe$Accuracy_TestSet_rfe)

t.resample.rfe <- t.test(df.results.rfe$ResamplingAccuracy_full, df.results.rfe$ResamplingAccuracy_rfe, paired = TRUE)
t.resample.rfe








