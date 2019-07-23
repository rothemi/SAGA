#############################################################################################################################################
#################################### Construction of SAGA-SVM Classifier ####################################################################
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
registerDoMC(cores = 60)    


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
#### I. TestSet 1 ###########################################################################################################################
#############################################################################################################################################
#############################################################################################################################################


##############################################################################################################
#### 1. Divide into Training and Test set1 ###################################################################
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
summary(fselect.1)
matrix.train.1 <-matrix.train.1[fselect.1,]

##############################################################################################################
#### 3. SVM: FULL MODEL (1151 predictors) ####################################################################
##############################################################################################################

matrix.train.1 <- (t(matrix.train.1))
labels.train.1 <- as.factor(pData.train.1$Class)


#### 3.1 SetUp SVMrad ########################################################################################
##############################################################################################################
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...)) # calculates Accuracy, Sens, Spec, ROC, kappa of external resamples

## create 200 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for comparison
set.seed(123)
index.1 <- createMultiFolds(labels.train.1, k=10, times = 20)  

fullCtrl.1 <- trainControl(method = "repeatedcv",
                         repeats = 20,
                         summaryFunction = fiveStats,
                         classProbs = TRUE,
                         index = index.1,
                         allowParallel = TRUE)

set.seed(721)
svmFull.1 <- train(matrix.train.1,labels.train.1,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.1)

svmFull.1  

# Resampling results for the FullModel with 107 samples 1151 predictors: 
# final tuning parameters: sigma = 0.0006816859, C = 0.5 
#    ROC       Sens   Spec       Accuracy   Kappa
# 0.9508452  0.84575  0.9377381  0.8978712  0.7885724

# ROC-based predictor importance
varImp.1 <- varImp(svmFull.1)$importance
varImp.1 <- varImp.1[order(varImp.1$transforming,decreasing = T),]
varImp.1 <- cbind(varImp.1,Annotation[row.names(varImp.1),])
write.table(varImp.1, file = "TestSet1_VariableImportance_ROC.txt", sep="\t",col.names=NA)

##############################################################################################################
#### 4. SVM-RFE ##############################################################################################
##############################################################################################################

#### 4.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

## set number of features to test (subsetSizes
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
                                 allowParallel = TRUE)                 # parallel on 60/90 cores 

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

rfe.1
write.table(rfe.1$results, file = "TestSet1_Results_Split7S123_IQR8_CVn20S123.txt", sep="\t",col.names=NA)

optFeatures.1 <- cbind(rfe.1$optVariables, Annotation[rfe.1$optVariables,])
write.table(optFeatures.1, file = "TestSet1_OptVars_Split7S123_IQR8_CVn20S123.txt", sep="\t",col.names=NA)

# Plot Accuracy over FeatureNumber
trellis.par.set(caretTheme())
plot(rfe.1, type = c("g", "o"))
plot(rfe.1, type = c("g", "o"), xlim = c(0,61))

# Plot ROC over FeatureNumber
plot(rfe.1$results$Variables,rfe.1$results$ROC, xlim = c(0,201),type = "o",panel.first = grid())

                         
#### 4.4.compare resampling performances #####################################################################
##############################################################################################################

rfeResamples.1 <- resamples(list("SVM_full.1" = svmFull.1,"SVM_RFE.1" = rfe.1))
sink("TestSet1_Resamples_Split7S123_IQR8_CVn20S123.txt", append = TRUE)
summary(rfeResamples.1)
sink()

modelDifferences.1 <- diff(rfeResamples.1)  # paired t-test for H0: difference = 0 between the different models. 
sink("TestSet1_ModelDifferences_Split7S123_IQR8_CVn20S123.txt", append = TRUE)
summary(modelDifferences.1)
sink()


##### 5  PCA on optimal variables ##############################################################
################################################################################################
matrix.train.opt.1     <- t(matrix.train.1[,row.names(optFeatures.1)])   # subset to the optVars
matrix.test.opt.1      <- matrix.test.1[row.names(optFeatures.1),]       # subset to the optVars
matrix.opt.1           <- cbind(matrix.train.opt.1, matrix.test.opt.1)   # combine training and test set
pData.opt.1            <- rbind(pData.train.1,pData.test.1)
pData.opt.1$Set        <- ifelse(row.names(pData.opt.1)%in% row.names(pData.train.1),"train","test")
pData.opt.1$Design_Color <- ifelse(pData.opt.1$Set == "train", pData.opt.1$Design_Color,"#000000")   # TestSet samples are black

pca.1         <- prcomp(t(matrix.opt.1))           
plot(pca.1$x, pch=16, col=pData.opt.1$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.1$Design_Color), pch=16, bty="n", cex=0.8)

##### compared to PCA on all variables ########################################################
matrix.full.1 <- rbind(matrix.train.1, t(matrix.test.1[colnames(matrix.train.1),]))
pca.full.1    <- prcomp(matrix.full.1)           
plot(pca.full.1$x, pch=16, col=pData.opt.1$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.1$Design_Color), pch=16, bty="n", cex=0.8)


##### 6  Train SVM on optVar FeatureSet and predict independent TestSet.1 ######################
################################################################################################

set.seed(721)
svmOpt.1  <- train(t(matrix.train.opt.1),labels.train.1,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.1)   # use same parameters as for the SVMfull model (10CVn20/index.1)
svmOpt.1  

# Resampling results for the optVars Model with 107 samples and 7 predictors 
# final tuning parameters: sigma = 0.3736986, C = 1073741824 
#    ROC       Sens   Spec       Accuracy   Kappa
# 0.9797976  0.95425  0.9582143  0.9568106  0.9118657  
# note: these values are better than the external resampling results (above) due to 
# positive bias since the optimal feature were determined for this set

##### Predict TestSet 1 with optimal predictors
Prediction_optVars.1 <- predict(svmOpt.1,t(matrix.test.opt.1), type = "prob")
Prediction_optVars.1$Prediction_optVars.1 <- ifelse(Prediction_optVars.1$transforming>0.50,"transforming","untransforming")
Prediction_optVars.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_optVars.1)
write.table(Prediction_optVars.1, file = paste("TestSet1_Predictions_optVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet 1 with all 1151 predictors
matrix.test.1.full <- matrix.test.1[row.names(t(matrix.train.1)),]
Prediction_SVM_full.1 <- predict(svmFull.1,t(matrix.test.1.full), type = "prob")
Prediction_SVM_full.1$Prediction_SVM_full <- ifelse(Prediction_SVM_full.1$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.1 <- cbind(pData.test.1[,c(1:3)],TrueLabel=pData.test.1$Class,Prediction_SVM_full.1)
write.table(Prediction_SVM_full.1, file = paste("TestSet1_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 7 Performance of opVars Classifier on TestSet 1 ########################################################
#############################################################################################################

#### 7.1 Confusion matrix  ##################################################################################  
sink("TestSet1_ConfusionMatrix_optVars.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.1$Prediction_optVars.1), as.factor(Prediction_optVars.1$TrueLabel))
sink()

sink("TestSet1_ConfusionMatrix_allVars.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.1$Prediction_SVM_full), as.factor(Prediction_SVM_full.1$TrueLabel))
sink()


#### 7.2.ROC on probability "transforming" on TestSet 1 #####################################################

Prediction_optVars.1$Class <- as.factor(ifelse(Prediction_optVars.1$TrueLabel == "transforming","transforming","nontransforming"))

roc1 <- roc(Prediction_optVars.1$Class,                    # response vector (factor or character)
            Prediction_optVars.1$transforming,             # predictor vector (numeric)
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
#### II. TestSet 02 #########################################################################################################################
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

intersect(row.names(pData.test.1),row.names(pData.test.2))

unique(c(row.names(pData.test.1),row.names(pData.test.2),row.names(pData.test.3),row.names(pData.test.4),row.names(pData.test.5),
      row.names(pData.test.6),row.names(pData.test.7),row.names(pData.test.8),row.names(pData.test.9),row.names(pData.test.10)))


##############################################################################################################
#### 3. nonspecific feature prefiltering TestSet 02  #########################################################
##############################################################################################################

fselect.2  <- genefilter(matrix.train.2, filterfun(f1))
summary(fselect.2)
matrix.train.2 <-matrix.train.2[fselect.2,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 02 ###########################################################################
##############################################################################################################

matrix.train.2 <- (t(matrix.train.2))
labels.train.2 <- as.factor(pData.train.2$Class)


#### 4.1 SetUp SVMrad for full model TestSet 02 ##############################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
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
                   tuneLength = 50,
                   trControl = fullCtrl.2)

svmFull.2  

# Resampling results for the FullModel with 107 samples 1201 predictors: 
# final tuning parameters: sigma = 0.000713516, C = 64 
#    ROC       Sens   Spec       Accuracy   Kappa
# 0.9688571  0.86125  0.9110714  0.8902045  0.7739778


##############################################################################################################
#### 5. SVM-RFE TestSet 02 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.2      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.2,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.2$functions         <- caretFuncs
outerctrl.2$functions$summary <- fiveStats


#### 5.3.SVM-RFE TestSet 02 #####################################################################################
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

rfe.2
write.table(rfe.2$results, file = "Results_TestSet02.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.2, type = c("g", "o"))
plot(rfe.2, type = c("g", "o"), xlim = c(0,61))

optFeatures.2 <- cbind(rfe.2$optVariables, Annotation[rfe.2$optVariables,])
write.table(optFeatures.2, file = "TestSet2_OptVars_Split7S1440_IQR8_CVn20S1234.txt", sep="\t",col.names=NA)

#### 5.4.compare resampling performances TestSet 02 #############################################################
##############################################################################################################

rfeResamples.2 <- resamples(list("SVM_full.2" = svmFull.2,"SVM_RFE.2" = rfe.2))
sink("TestSet02_Resamples.txt", append = TRUE)
summary(rfeResamples.2)
sink()

modelDifferences.2 <- diff(rfeResamples.2)  # paired t-test for H0: difference = 0 between the different models. 
sink("TestSet02_ModelDifferences.txt", append = TRUE)
summary(modelDifferences.2)
sink()


################################################################################################
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TestSet02 ##############################
################################################################################################

matrix.train.GA.2 <- matrix.train.2[,rfe.2$optVariables]   # subset fot the 35 optVars

#### 6.0 Global SVM-GA Paramneters: starting population   ######################################
################################################################################################

svmGA <- caretGA

initial33 <- function (vars, popSize, ...) {
                                            x <- matrix(NA, nrow = popSize, ncol = vars)
                                            probs <- rep(0.6666,length = popSize)
                                            for (i in 1:popSize) {
                                                x[i, ] <- sample(0:1, replace = TRUE, size = vars, prob = c(probs[i], 
                                                                1 - probs[i]))
                                                                  }
                                            var_count <- apply(x, 1, sum)
                                            if (any(var_count == 0)) {
                                            for (i in which(var_count == 0)) {
                                            x[i, ] <- sample(0:1, replace = TRUE, size = vars)
                                                                              }
                                                                      }
                                             x
                                            }

environment(initial33) <- asNamespace('caret')
svmGA$initial <- initial33


#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(1714)
index.GA.2 <- createMultiFolds(labels.train.2, k=10, times = 5)  

set.seed(123)
seeds.GA <- vector(mode = "integer", length = length(index.GA.2)+1)    # B+1 elements where B is the number of resamples = 51
for(i in 1:length(index.GA.2)+1) seeds.GA[[i]] <- sample.int(10000, 1)


outerctrl.GA.2 <- gafsControl(functions = svmGA,
                            method = "repeatedcv", repeats = 5,
                            index = index.GA.2,                                       
                            seeds = seeds.GA,                                      
                            returnResamp="all", 
                            verbose = TRUE,
                            maximize = c(internal = TRUE,
                                         external = TRUE),
                            allowParallel = TRUE)                                  




#### 6.3 run GA  ###############################################################################
################################################################################################

system.time(GA.2<- gafs(matrix.train.GA.2, labels.train.2, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.2,
                        ## Now we pass options to `train` via "svmGA":               
                        metric = "Accuracy",
                        method = "svmRadial",
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 6.1.1 analyze results  ##################################################################
GA.2
GA.2$optIter
GA.2$optVariables
optVars.GA <- Annotation[GA.2$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.2_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)


##### 6  PCA on optimal variables TestSet 02 ####################################################################
##############################################################################################################
matrix.train.opt.2     <- t(matrix.train.2[,row.names(optFeatures.2)])
matrix.test.opt.2      <- matrix.test.2[row.names(optFeatures.2),]
matrix.opt.2           <- cbind(matrix.train.opt.2, matrix.test.opt.2)
pData.opt.2            <- rbind(pData.train.2,pData.test.2)
pData.opt.2$Set        <- ifelse(row.names(pData.opt.2)%in% row.names(pData.train.2),"train","test")
pData.opt.2$Design_Color <- ifelse(pData.opt.2$Set == "train", pData.opt.2$Design_Color,"#000000")

pca.2         <- prcomp(t(matrix.opt.2))           
plot(pca.2$x, pch=16, col=pData.opt.2$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.2$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.2 and predict independent TestSet.2 ##################################
##############################################################################################################

optCtrl.2 <- trainControl(method = "repeatedcv",repeats = 5,
                          summaryFunction = fiveStats,
                          classProbs = TRUE,
                          index = index.2,
                          allowParallel = TRUE)

set.seed(721)
svmOpt.2  <- train(t(matrix.train.opt.2),labels.train.2,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = optCtrl.2)
svmOpt.2  

Prediction_optVars.2 <- predict(svmOpt.2,t(matrix.test.opt.2), type = "prob")
Prediction_optVars.2$Prediction_optVars.2 <- ifelse(Prediction_optVars.2$transforming>0.50,"transforming","untransforming")
Prediction_optVars.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_optVars.2)
write.table(Prediction_optVars.2, file = paste("Predictions_TestSet 02.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

##### Predict TestSet 2 with all predictors
matrix.test.2.full <- matrix.test.2[row.names(t(matrix.train.2)),]
Prediction_SVM_full.2 <- predict(svmFull.2,t(matrix.test.2.full), type = "prob")
Prediction_SVM_full.2$Prediction_SVM_full <- ifelse(Prediction_SVM_full.2$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.2 <- cbind(pData.test.2[,c(1:3)],TrueLabel=pData.test.2$Class,Prediction_SVM_full.2)
write.table(Prediction_SVM_full.2, file = paste("Predictions_TestSet1_SVMfull_1151Vars_Accuracy.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 02 ##########################################
#############################################################################################################

#### 8.2 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet02.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.2$Prediction_optVars.2), as.factor(Prediction_optVars.2$TrueLabel))
sink()

sink("ConfusionMatrix_TestSet1_SVMfull_1151Vars_Accuracy.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.2$Prediction_SVM_full), as.factor(Prediction_SVM_full.2$TrueLabel))
sink()

#### 6.2.2 ROC on probability "transforming" TestSet 02 ########################################################

Prediction_optVars.2$Class <- as.factor(ifelse(Prediction_optVars.2$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet02.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.2$Class,                    # response vector (factor or character)
            Prediction_optVars.2$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=0.5)
dev.off()

Prediction_SVM_full.1$Class <- as.factor(ifelse(Prediction_SVM_full.1$TrueLabel == "transforming","transforming","nontransforming"))
roc2 <- roc(Prediction_SVM_full.1$Class,                    # response vector (factor or character)
            Prediction_SVM_full.1$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=0.5)



#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 03 #########################################################################################################################
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

intersect(row.names(pData.test.1),row.names(pData.test.3))

##############################################################################################################
#### 3. nonspecific feature prefiltering TestSet 03  ############################################################
##############################################################################################################
fselect.3  <- genefilter(matrix.train.3, filterfun(f1))
summary(fselect.3)
matrix.train.3 <-matrix.train.3[fselect.3,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 03 ##############################################################################
##############################################################################################################

matrix.train.3 <- (t(matrix.train.3))
labels.train.3 <- as.factor(pData.train.3$Class)


#### 4.1 SetUp SVMrad for full model TestSet 03 #################################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(12345)
index.3 <- createMultiFolds(labels.train.3, k=10, times = 20)  

fullCtrl.3 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.3,
                           allowParallel = TRUE)

set.seed(721)
svmFull.3 <- train(matrix.train.3,labels.train.3,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.3)

svmFull.3  

##############################################################################################################
#### 5. SVM-RFE TestSet 03 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.3      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.3,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.3$functions         <- caretFuncs
outerctrl.3$functions$summary <- fiveStats


#### 5.3.SVM-RFE TestSet 03 #####################################################################################
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

rfe.3
write.table(rfe.3$results, file = "Results_TestSet03.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.3, type = c("g", "o"))
plot(rfe.3, type = c("g", "o"), xlim = c(0,61))

optFeatures.3 <- cbind(rfe.3$optVariables, Annotation[rfe.3$optVariables,])
write.table(optFeatures.3, file = "OptVars_rfe152_TestSet03.txt", sep="\t",col.names=NA)

#### 5.4.compare resampling performances TestSet 03 #############################################################
##############################################################################################################

rfeResamples.3 <- resamples(list("SVM_full.3" = svmFull.3,"SVM_RFE.3" = rfe.3))
sink("Resamples_rfe152_TestSet03.txt", append = TRUE)
summary(rfeResamples.3)
sink()

modelDifferences.3 <- diff(rfeResamples.3)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet03.txt", append = TRUE)
summary(modelDifferences.3)
sink()


################################################################################################
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TestSet03 ##############################
################################################################################################
matrix.train.GA.3 <- matrix.train.3[,rfe.3$optVariables]   # subset fot the 36 optVars

#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
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

#### 6.2 run GA  ###############################################################################
################################################################################################

system.time(GA.3<- gafs(matrix.train.GA.3, labels.train.3, 
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

### 6.1.1 analyze results  ##################################################################
GA.3
GA.3$optIter
GA.3$optVariables
optVars.GA <- Annotation[GA.3$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.3_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)



##### 6  PCA on optimal variables TestSet 03 ####################################################################
##############################################################################################################
matrix.train.opt.3     <- t(matrix.train.3[,row.names(optFeatures.3)])
matrix.test.opt.3      <- matrix.test.3[row.names(optFeatures.3),]
matrix.opt.3           <- cbind(matrix.train.opt.3, matrix.test.opt.3)
pData.opt.3            <- rbind(pData.train.3,pData.test.3)
pData.opt.3$Set        <- ifelse(row.names(pData.opt.3)%in% row.names(pData.train.3),"train","test")
pData.opt.3$Design_Color <- ifelse(pData.opt.3$Set == "train", pData.opt.3$Design_Color,"#000000")

pca.3         <- prcomp(t(matrix.opt.3))           
plot(pca.3$x, pch=16, col=pData.opt.3$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.3$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.3 and predict independent TestSet.3 ##################################
##############################################################################################################

set.seed(721)
svmOpt.3  <- train(t(matrix.train.opt.3),labels.train.3,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.3)
svmOpt.3  

Prediction_optVars.3 <- predict(svmOpt.3,t(matrix.test.opt.3), type = "prob")
Prediction_optVars.3$Prediction_optVars.3 <- ifelse(Prediction_optVars.3$transforming>0.50,"transforming","untransforming")
Prediction_optVars.3 <- cbind(pData.test.3[,c(1:3)],TrueLabel=pData.test.3$Class,Prediction_optVars.3)
write.table(Prediction_optVars.3, file = paste("Predictions_TestSet 03.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 03 #############################################
#############################################################################################################

#### 8.3 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet03.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.3$Prediction_optVars.3), as.factor(Prediction_optVars.3$TrueLabel))
sink()

#### 6.3.3 ROC on probability "transforming" TestSet 03 ########################################################

Prediction_optVars.3$Class <- as.factor(ifelse(Prediction_optVars.3$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet03.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.3$Class,                    # response vector (factor or character)
            Prediction_optVars.3$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()


#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 04 #########################################################################################################################
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
#### 3. nonspecific feature prefiltering TestSet 04 ##########################################################
##############################################################################################################
fselect.4  <- genefilter(matrix.train.4, filterfun(f1))
summary(fselect.4)
matrix.train.4 <-matrix.train.4[fselect.4,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 04 ###########################################################################
##############################################################################################################

matrix.train.4 <- (t(matrix.train.4))
labels.train.4 <- as.factor(pData.train.4$Class)


#### 4.1 SetUp SVMrad for full model TestSet 04 ##############################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
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
#### 5. SVM-RFE TestSet 04 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
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


#### 5.4.SVM-RFE TestSet 04 ##################################################################################
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

rfe.4
write.table(rfe.4$results, file = "Results_TestSet04.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.4, type = c("g", "o"))
plot(rfe.4, type = c("g", "o"), xlim = c(0,61))

optFeatures.4 <- cbind(rfe.4$optVariables, Annotation[rfe.4$optVariables,])
write.table(optFeatures.4, file = "OptVars_rfe152_TestSet03.txt", sep="\t",col.names=NA)

#### 5.4.compare resampling performances TestSet 04 ##########################################################
##############################################################################################################

rfeResamples.4 <- resamples(list("SVM_full.4" = svmFull.4,"SVM_RFE.4" = rfe.4))
sink("Resamples_rfe152_TestSet03.txt", append = TRUE)
summary(rfeResamples.4)
sink()

modelDifferences.4 <- diff(rfeResamples.4)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet03.txt", append = TRUE)
summary(modelDifferences.4)
sink()


##### 6  PCA on optimal variables TestSet 04 ####################################################################
##############################################################################################################
matrix.train.opt.4     <- t(matrix.train.4[,row.names(optFeatures.4)])
matrix.test.opt.4      <- matrix.test.4[row.names(optFeatures.4),]
matrix.opt.4           <- cbind(matrix.train.opt.4, matrix.test.opt.4)
pData.opt.4            <- rbind(pData.train.4,pData.test.4)
pData.opt.4$Set        <- ifelse(row.names(pData.opt.4)%in% row.names(pData.train.4),"train","test")
pData.opt.4$Design_Color <- ifelse(pData.opt.4$Set == "train", pData.opt.4$Design_Color,"#000000")

pca.4         <- prcomp(t(matrix.opt.4))           
plot(pca.4$x, pch=16, col=pData.opt.4$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.4$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.4 and predict independent TestSet.4 ##################################
##############################################################################################################

optCtrl.4 <- trainControl(method = "repeatedcv",repeats = 5,
                          summaryFunction = fiveStats,
                          classProbs = TRUE,
                          index = index.4,
                          allowParallel = TRUE)

set.seed(721)
svmOpt.4  <- train(t(matrix.train.opt.4),labels.train.4,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = optCtrl.4)
svmOpt.4  

Prediction_optVars.4 <- predict(svmOpt.4,t(matrix.test.opt.4), type = "prob")
Prediction_optVars.4$Prediction_optVars.4 <- ifelse(Prediction_optVars.4$transforming>0.50,"transforming","untransforming")
Prediction_optVars.4 <- cbind(pData.test.4[,c(1:3)],TrueLabel=pData.test.4$Class,Prediction_optVars.4)
write.table(Prediction_optVars.4, file = paste("Predictions_TestSet 03.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 04 ##########################################
#############################################################################################################

#### 8.4 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet03.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.4$Prediction_optVars.4), as.factor(Prediction_optVars.4$TrueLabel))
sink()

#### 6.4.4 ROC on probability "transforming" TestSet 04 ########################################################

Prediction_optVars.4$Class <- as.factor(ifelse(Prediction_optVars.4$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet03.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.4$Class,                    # response vector (factor or character)
            Prediction_optVars.4$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()




#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 05 #########################################################################################################################
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
#### 3. nonspecific feature prefiltering TestSet 05 ##########################################################
##############################################################################################################
fselect.5  <- genefilter(matrix.train.5, filterfun(f1))
summary(fselect.5)
matrix.train.5 <-matrix.train.5[fselect.5,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 05 ###########################################################################
##############################################################################################################

matrix.train.5 <- (t(matrix.train.5))
labels.train.5 <- as.factor(pData.train.5$Class)


#### 4.1 SetUp SVMrad for full model TestSet 05 ##############################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.5 <- createMultiFolds(labels.train.5, k=10, times = 20)  

fullCtrl.5 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.5,
                           allowParallel = TRUE)

set.seed(721)
svmFull.5 <- train(matrix.train.5,labels.train.5,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.5)

svmFull.5  

##############################################################################################################
#### 5. SVM-RFE TestSet 05 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
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


#### 5.5.SVM-RFE TestSet 05 ##################################################################################
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

rfe.5
write.table(rfe.5$results, file = "Results_TestSet05.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.5, type = c("g", "o"))
plot(rfe.5, type = c("g", "o"), xlim = c(0,61))

optFeatures.5 <- cbind(rfe.5$optVariables, Annotation[rfe.5$optVariables,])
write.table(optFeatures.5, file = "OptVars_rfe152_TestSet03.txt", sep="\t",col.names=NA)

#### 5.5.compare resampling performances TestSet 05 ##########################################################
##############################################################################################################

rfeResamples.5 <- resamples(list("SVM_full.5" = svmFull.5,"SVM_RFE.5" = rfe.5))
sink("Resamples_rfe152_TestSet03.txt", append = TRUE)
summary(rfeResamples.5)
sink()

modelDifferences.5 <- diff(rfeResamples.5)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet03.txt", append = TRUE)
summary(modelDifferences.5)
sink()


##### 6  PCA on optimal variables TestSet 05 ####################################################################
##############################################################################################################
matrix.train.opt.5     <- t(matrix.train.5[,row.names(optFeatures.5)])
matrix.test.opt.5      <- matrix.test.5[row.names(optFeatures.5),]
matrix.opt.5           <- cbind(matrix.train.opt.5, matrix.test.opt.5)
pData.opt.5            <- rbind(pData.train.5,pData.test.5)
pData.opt.5$Set        <- ifelse(row.names(pData.opt.5)%in% row.names(pData.train.5),"train","test")
pData.opt.5$Design_Color <- ifelse(pData.opt.5$Set == "train", pData.opt.5$Design_Color,"#000000")

pca.5         <- prcomp(t(matrix.opt.5))           
plot(pca.5$x, pch=16, col=pData.opt.5$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.5$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.5 and predict independent TestSet.5 ##################################
##############################################################################################################

optCtrl.5 <- trainControl(method = "repeatedcv",repeats = 5,
                          summaryFunction = fiveStats,
                          classProbs = TRUE,
                          index = index.5,
                          allowParallel = TRUE)

set.seed(721)
svmOpt.5  <- train(t(matrix.train.opt.5),labels.train.5,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = optCtrl.5)
svmOpt.5  

Prediction_optVars.5 <- predict(svmOpt.5,t(matrix.test.opt.5), type = "prob")
Prediction_optVars.5$Prediction_optVars.5 <- ifelse(Prediction_optVars.5$transforming>0.50,"transforming","untransforming")
Prediction_optVars.5 <- cbind(pData.test.5[,c(1:3)],TrueLabel=pData.test.5$Class,Prediction_optVars.5)
write.table(Prediction_optVars.5, file = paste("Predictions_TestSet 03.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 05 ##########################################
#############################################################################################################

#### 8.5 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet03.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.5$Prediction_optVars.5), as.factor(Prediction_optVars.5$TrueLabel))
sink()

#### 6.5.5 ROC on probability "transforming" TestSet 05 ########################################################

Prediction_optVars.5$Class <- as.factor(ifelse(Prediction_optVars.5$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet03.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.5$Class,                    # response vector (factor or character)
            Prediction_optVars.5$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()




#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 06 #########################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
set.seed(678)
split.6 <- createDataPartition(as.factor(pData$Class), p = .7, list = FALSE) # group-stratified sampling: 107 Training, 45 TestSamples

matrix.train.6 <- eset.batch[, split.6] 
matrix.test.6  <- eset.batch[,-split.6] 
pData.train.6  <- pData[split.6,]
pData.test.6   <- pData[-split.6,]

table(pData.train.6$Design)
table(pData.test.6$Design)

##############################################################################################################
#### 3. nonspecific feature prefiltering TestSet 05 ##########################################################
##############################################################################################################
fselect.6  <- genefilter(matrix.train.6, filterfun(f1))
summary(fselect.6)
matrix.train.6 <-matrix.train.6[fselect.6,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 06 ###########################################################################
##############################################################################################################

matrix.train.6 <- (t(matrix.train.6))
labels.train.6 <- as.factor(pData.train.6$Class)


#### 4.1 SetUp SVMrad for full model TestSet 06 ##############################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(678)
index.6 <- createMultiFolds(labels.train.6, k=10, times = 20)  

fullCtrl.6 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.6,
                           allowParallel = TRUE)

set.seed(721)
svmFull.6 <- train(matrix.train.6,labels.train.6,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.6)

svmFull.6  

##############################################################################################################
#### 5. SVM-RFE TestSet 05 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.6      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.6,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.6$functions         <- caretFuncs
outerctrl.6$functions$summary <- fiveStats


#### 5.6.SVM-RFE TestSet 06 ##################################################################################
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

rfe.6
write.table(rfe.6$results, file = "Results_TestSet05.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.6, type = c("g", "o"))
plot(rfe.6, type = c("g", "o"), xlim = c(0,61))

optFeatures.6 <- cbind(rfe.6$optVariables, Annotation[rfe.6$optVariables,])
write.table(optFeatures.6, file = "OptVars_rfe152_TestSet03.txt", sep="\t",col.names=NA)

#### 5.6.compare resampling performances TestSet 06 ##########################################################
##############################################################################################################

rfeResamples.6 <- resamples(list("SVM_full.6" = svmFull.6,"SVM_RFE.6" = rfe.6))
sink("Resamples_rfe152_TestSet03.txt", append = TRUE)
summary(rfeResamples.6)
sink()

modelDifferences.6 <- diff(rfeResamples.6)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet03.txt", append = TRUE)
summary(modelDifferences.6)
sink()


################################################################################################
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TestSet06 ##############################
################################################################################################

matrix.train.GA.6 <- matrix.train.6[,rfe.6$optVariables]   # subset fot the 33 optVars

#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
################################################################################################
set.seed(2144)
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

#### 6.3 run GA  ###############################################################################
################################################################################################

system.time(GA.6<- gafs(matrix.train.GA.6, labels.train.6, 
                                iters = 40,
                                popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                                gafsControl = outerctrl.GA.6,
                                ## Now we pass options to `train` via "svmGA":               
                                metric = "Accuracy",
                                method = "svmRadial",
                                tuneLength = 12,
                                trControl = trainControl(method = "repeatedcv",
                                                         repeats = 2,
                                                         allowParallel = FALSE)))

### 6.1.1 analyze results  ##################################################################
GA.6
GA.6$optIter
GA.6$optVariables
optVars.GA <- Annotation[GA.6$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.7_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)


##### 6  PCA on optimal variables TestSet 06 #################################################################
##############################################################################################################
matrix.train.opt.6     <- t(matrix.train.6[,row.names(optFeatures.6)])
matrix.test.opt.6      <- matrix.test.6[row.names(optFeatures.6),]
matrix.opt.6           <- cbind(matrix.train.opt.6, matrix.test.opt.6)
pData.opt.6            <- rbind(pData.train.6,pData.test.6)
pData.opt.6$Set        <- ifelse(row.names(pData.opt.6)%in% row.names(pData.train.6),"train","test")
pData.opt.6$Design_Color <- ifelse(pData.opt.6$Set == "train", pData.opt.6$Design_Color,"#000000")

pca.6         <- prcomp(t(matrix.opt.6))           
plot(pca.6$x, pch=16, col=pData.opt.6$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.6$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.6 and predict independent TestSet.6 ##################################
##############################################################################################################

optCtrl.6 <- trainControl(method = "repeatedcv",repeats = 5,
                          summaryFunction = fiveStats,
                          classProbs = TRUE,
                          index = index.6,
                          allowParallel = TRUE)

set.seed(721)
svmOpt.6  <- train(t(matrix.train.opt.6),labels.train.6,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = optCtrl.6)
svmOpt.6  

Prediction_optVars.6 <- predict(svmOpt.6,t(matrix.test.opt.6), type = "prob")
Prediction_optVars.6$Prediction_optVars.6 <- ifelse(Prediction_optVars.6$transforming>0.60,"transforming","untransforming")
Prediction_optVars.6 <- cbind(pData.test.6[,c(1:3)],TrueLabel=pData.test.6$Class,Prediction_optVars.6)
write.table(Prediction_optVars.6, file = paste("Predictions_TestSet 03.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 05 ##########################################
#############################################################################################################

#### 8.6 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet03.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.6$Prediction_optVars.6), as.factor(Prediction_optVars.6$TrueLabel))
sink()

#### 6.6.6 ROC on probability "transforming" TestSet 05 ########################################################

Prediction_optVars.6$Class <- as.factor(ifelse(Prediction_optVars.6$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet03.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.6$Class,                    # response vector (factor or character)
            Prediction_optVars.6$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()




#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 07 #########################################################################################################################
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
#### 3. nonspecific feature prefiltering TestSet 07 ##########################################################
##############################################################################################################
fselect.7  <- genefilter(matrix.train.7, filterfun(f1))
summary(fselect.7)
matrix.train.7 <-matrix.train.7[fselect.7,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 07 ###########################################################################
##############################################################################################################

matrix.train.7 <- (t(matrix.train.7))
labels.train.7 <- as.factor(pData.train.7$Class)


#### 4.1 SetUp SVMrad for full model TestSet 07 ##############################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(847)
index.7 <- createMultiFolds(labels.train.7, k=10, times = 20)  

fullCtrl.7 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.7,
                           allowParallel = TRUE)

set.seed(721)
svmFull.7 <- train(matrix.train.7,labels.train.7,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.7)

svmFull.7  

##############################################################################################################
#### 5. SVM-RFE TestSet 07 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.7      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.7,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.7$functions         <- caretFuncs
outerctrl.7$functions$summary <- fiveStats


#### 5.7.SVM-RFE TestSet 07 ##################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.7  <- rfe(matrix.train.7, labels.train.7, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.7,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.7
write.table(rfe.7$results, file = "Results_TestSet07.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.7, type = c("g", "o"))
plot(rfe.7, type = c("g", "o"), xlim = c(0,61))

optFeatures.7 <- cbind(rfe.7$optVariables, Annotation[rfe.7$optVariables,])
write.table(optFeatures.7, file = "OptVars_rfe152_TestSet03.txt", sep="\t",col.names=NA)

#### 5.7.compare resampling performances TestSet 07 ##########################################################
##############################################################################################################

rfeResamples.7 <- resamples(list("SVM_full.7" = svmFull.7,"SVM_RFE.7" = rfe.7))
sink("Resamples_rfe152_TestSet03.txt", append = TRUE)
summary(rfeResamples.7)
sink()

modelDifferences.7 <- diff(rfeResamples.7)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet03.txt", append = TRUE)
summary(modelDifferences.7)
sink()



################################################################################################
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TestSet07 ##############################
################################################################################################

matrix.train.GA.7 <- matrix.train.7[,rfe.7$optVariables]   # subset fot the 33 optVars

#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
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


system.time(GA.7<- gafs(matrix.train.GA.7, labels.train.7, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.7,
                        ## Now we pass options to `train` via "svmGA":               
                        metric = "Accuracy",
                        method = "svmRadial",
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 6.1.1 analyze results  ##################################################################
GA.7
GA.7$optIter
GA.7$optVariables
optVars.GA <- Annotation[GA.7$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.7_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)













##### 6  PCA on optimal variables TestSet 07 ####################################################################
##############################################################################################################
matrix.train.opt.7     <- t(matrix.train.7[,row.names(optFeatures.7)])
matrix.test.opt.7      <- matrix.test.7[row.names(optFeatures.7),]
matrix.opt.7           <- cbind(matrix.train.opt.7, matrix.test.opt.7)
pData.opt.7            <- rbind(pData.train.7,pData.test.7)
pData.opt.7$Set        <- ifelse(row.names(pData.opt.7)%in% row.names(pData.train.7),"train","test")
pData.opt.7$Design_Color <- ifelse(pData.opt.7$Set == "train", pData.opt.7$Design_Color,"#000000")

pca.7         <- prcomp(t(matrix.opt.7))           
plot(pca.7$x, pch=16, col=pData.opt.7$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.7$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.7 and predict independent TestSet.7 ##################################
##############################################################################################################

optCtrl.7 <- trainControl(method = "repeatedcv",repeats = 5,
                          summaryFunction = fiveStats,
                          classProbs = TRUE,
                          index = index.7,
                          allowParallel = TRUE)

set.seed(721)
svmOpt.7  <- train(t(matrix.train.opt.7),labels.train.7,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = optCtrl.7)
svmOpt.7  

Prediction_optVars.7 <- predict(svmOpt.7,t(matrix.test.opt.7), type = "prob")
Prediction_optVars.7$Prediction_optVars.7 <- ifelse(Prediction_optVars.7$transforming>0.70,"transforming","untransforming")
Prediction_optVars.7 <- cbind(pData.test.7[,c(1:3)],TrueLabel=pData.test.7$Class,Prediction_optVars.7)
write.table(Prediction_optVars.7, file = paste("Predictions_TestSet 03.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 07 ##########################################
#############################################################################################################

#### 8.7 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet03.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.7$Prediction_optVars.7), as.factor(Prediction_optVars.7$TrueLabel))
sink()

#### 6.7.7 ROC on probability "transforming" TestSet 07 ########################################################

Prediction_optVars.7$Class <- as.factor(ifelse(Prediction_optVars.7$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet03.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.7$Class,                    # response vector (factor or character)
            Prediction_optVars.7$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()



#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 08 #########################################################################################################################
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
#### 3. nonspecific feature prefiltering TestSet 08 ##########################################################
##############################################################################################################
fselect.8  <- genefilter(matrix.train.8, filterfun(f1))
summary(fselect.8)
matrix.train.8 <-matrix.train.8[fselect.8,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 08 ###########################################################################
##############################################################################################################

matrix.train.8 <- (t(matrix.train.8))
labels.train.8 <- as.factor(pData.train.8$Class)


#### 4.1 SetUp SVMrad for full model TestSet 08 ##############################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
index.8 <- createMultiFolds(labels.train.8, k=10, times = 20)  

fullCtrl.8 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.8,
                           allowParallel = TRUE)

set.seed(721)
svmFull.8 <- train(matrix.train.8,labels.train.8,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.8)

svmFull.8  

##############################################################################################################
#### 5. SVM-RFE TestSet 08 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.8      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.8,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.8$functions         <- caretFuncs
outerctrl.8$functions$summary <- fiveStats


#### 5.8.SVM-RFE TestSet 08 ##################################################################################
##############################################################################################################

# 200 outer resamples x 52 subsetSizes = 10,400 Predictions of HeldOuts x 30 internal resamples for tuning at each outer resampling iteration = 312,000 models
# this will take around 2 hours on an AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM

system.time(rfe.8  <- rfe(matrix.train.8, labels.train.8, 
                          sizes=FeatureNumbers,
                          rfeControl=outerctrl.8,
                          metric = "Accuracy",
                          ## Options to train()
                          method="svmRadial",
                          tuneLength = 20,
                          trControl = innerctrl))

rfe.8
write.table(rfe.8$results, file = "Results_TestSet08.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.8, type = c("g", "o"))
plot(rfe.8, type = c("g", "o"), xlim = c(0,61))

optFeatures.8 <- cbind(rfe.8$optVariables, Annotation[rfe.8$optVariables,])
write.table(optFeatures.8, file = "OptVars_rfe152_TestSet03.txt", sep="\t",col.names=NA)

#### 5.8.compare resampling performances TestSet 08 ##########################################################
##############################################################################################################

rfeResamples.8 <- resamples(list("SVM_full.8" = svmFull.8,"SVM_RFE.8" = rfe.8))
sink("Resamples_rfe152_TestSet03.txt", append = TRUE)
summary(rfeResamples.8)
sink()

modelDifferences.8 <- diff(rfeResamples.8)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet03.txt", append = TRUE)
summary(modelDifferences.8)
sink()


################################################################################################
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TestSet08 ##############################
################################################################################################

matrix.train.GA.8 <- matrix.train.8[,rfe.8$optVariables]   # subset fot the 33 optVars

#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
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


system.time(GA.8<- gafs(matrix.train.GA.8, labels.train.8, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.8,
                        ## Now we pass options to `train` via "svmGA":               
                        metric = "Accuracy",
                        method = "svmRadial",
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 6.1.1 analyze results  ##################################################################
GA.8
GA.8$optIter
GA.8$optVariables
optVars.GA <- Annotation[GA.8$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.8_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)


##### 6  PCA on optimal variables TestSet 08 ####################################################################
##############################################################################################################
matrix.train.opt.8     <- t(matrix.train.8[,row.names(optFeatures.8)])
matrix.test.opt.8      <- matrix.test.8[row.names(optFeatures.8),]
matrix.opt.8           <- cbind(matrix.train.opt.8, matrix.test.opt.8)
pData.opt.8            <- rbind(pData.train.8,pData.test.8)
pData.opt.8$Set        <- ifelse(row.names(pData.opt.8)%in% row.names(pData.train.8),"train","test")
pData.opt.8$Design_Color <- ifelse(pData.opt.8$Set == "train", pData.opt.8$Design_Color,"#000000")

pca.8         <- prcomp(t(matrix.opt.8))           
plot(pca.8$x, pch=16, col=pData.opt.8$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.8$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.8 and predict independent TestSet.8 ##################################
##############################################################################################################

optCtrl.8 <- trainControl(method = "repeatedcv",repeats = 5,
                          summaryFunction = fiveStats,
                          classProbs = TRUE,
                          index = index.8,
                          allowParallel = TRUE)

set.seed(721)
svmOpt.8  <- train(t(matrix.train.opt.8),labels.train.8,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 20,
                   trControl = optCtrl.8)
svmOpt.8  

Prediction_optVars.8 <- predict(svmOpt.8,t(matrix.test.opt.8), type = "prob")
Prediction_optVars.8$Prediction_optVars.8 <- ifelse(Prediction_optVars.8$transforming>0.80,"transforming","untransforming")
Prediction_optVars.8 <- cbind(pData.test.8[,c(1:3)],TrueLabel=pData.test.8$Class,Prediction_optVars.8)
write.table(Prediction_optVars.8, file = paste("Predictions_TestSet 03.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 08 ##########################################
#############################################################################################################

#### 8.8 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet03.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.8$Prediction_optVars.8), as.factor(Prediction_optVars.8$TrueLabel))
sink()

#### 6.8.8 ROC on probability "transforming" TestSet 08 ########################################################

Prediction_optVars.8$Class <- as.factor(ifelse(Prediction_optVars.8$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet03.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.8$Class,                    # response vector (factor or character)
            Prediction_optVars.8$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()


#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 09 #########################################################################################################################
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

intersect(row.names(pData.test.1),row.names(pData.test.9))

##############################################################################################################
#### 3. nonspecific feature prefiltering TestSet 09  #########################################################
##############################################################################################################
fselect.9  <- genefilter(matrix.train.9, filterfun(f1))
summary(fselect.9)
matrix.train.9 <-matrix.train.9[fselect.9,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 09 ###########################################################################
##############################################################################################################

matrix.train.9 <- (t(matrix.train.9))
labels.train.9 <- as.factor(pData.train.9$Class)


#### 4.1 SetUp SVMrad for full model TestSet 09 #################################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(4567)
index.9 <- createMultiFolds(labels.train.9, k=10, times = 20)  

fullCtrl.9 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.9,
                           allowParallel = TRUE)

set.seed(721)
svmFull.9 <- train(matrix.train.9,labels.train.9,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.9)

svmFull.9  

##############################################################################################################
#### 5. SVM-RFE TestSet 09 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.9      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.9,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.9$functions         <- caretFuncs
outerctrl.9$functions$summary <- fiveStats


#### 5.9.SVM-RFE TestSet 09 #####################################################################################
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

rfe.9
write.table(rfe.9$results, file = "Results_TestSet09.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.9, type = c("g", "o"))
plot(rfe.9, type = c("g", "o"), xlim = c(0,61))

optFeatures.9 <- cbind(rfe.9$optVariables, Annotation[rfe.9$optVariables,])
write.table(optFeatures.9, file = "OptVars_rfe152_TestSet09.txt", sep="\t",col.names=NA)

#### 5.4.compare resampling performances TestSet 09 #############################################################
##############################################################################################################

rfeResamples.9 <- resamples(list("SVM_full.9" = svmFull.9,"SVM_RFE.9" = rfe.9))
sink("Resamples_rfe152_TestSet09.txt", append = TRUE)
summary(rfeResamples.9)
sink()

modelDifferences.9 <- diff(rfeResamples.9)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet09.txt", append = TRUE)
summary(modelDifferences.9)
sink()


##### 6  PCA on optimal variables TestSet 09 ####################################################################
##############################################################################################################
matrix.train.opt.9     <- t(matrix.train.9[,row.names(optFeatures.9)])
matrix.test.opt.9      <- matrix.test.9[row.names(optFeatures.9),]
matrix.opt.9           <- cbind(matrix.train.opt.9, matrix.test.opt.9)
pData.opt.9            <- rbind(pData.train.9,pData.test.9)
pData.opt.9$Set        <- ifelse(row.names(pData.opt.9)%in% row.names(pData.train.9),"train","test")
pData.opt.9$Design_Color <- ifelse(pData.opt.9$Set == "train", pData.opt.9$Design_Color,"#000000")

pca.9         <- prcomp(t(matrix.opt.9))           
plot(pca.9$x, pch=16, col=pData.opt.9$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.9$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.9 and predict independent TestSet.9 ##################################
##############################################################################################################

set.seed(721)
svmOpt.9  <- train(t(matrix.train.opt.9),labels.train.9,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.9)
svmOpt.9  

Prediction_optVars.9 <- predict(svmOpt.9,t(matrix.test.opt.9), type = "prob")
Prediction_optVars.9$Prediction_optVars.9 <- ifelse(Prediction_optVars.9$transforming>0.50,"transforming","untransforming")
Prediction_optVars.9 <- cbind(pData.test.9[,c(1:3)],TrueLabel=pData.test.9$Class,Prediction_optVars.9)
write.table(Prediction_optVars.9, file = paste("Predictions_TestSet 09.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 09 #############################################
#############################################################################################################

#### 8.9 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet09.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.9$Prediction_optVars.9), as.factor(Prediction_optVars.9$TrueLabel))
sink()

#### 6.9.9 ROC on probability "transforming" TestSet 09 ########################################################

Prediction_optVars.9$Class <- as.factor(ifelse(Prediction_optVars.9$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet09.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.9$Class,                    # response vector (factor or character)
            Prediction_optVars.9$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()






#############################################################################################################################################
#############################################################################################################################################
#### II. TestSet 10 #########################################################################################################################
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

intersect(row.names(pData.test.1),row.names(pData.test.10))

##############################################################################################################
#### 3. nonspecific feature prefiltering TestSet 10  #########################################################
##############################################################################################################
fselect.10  <- genefilter(matrix.train.10, filterfun(f1))
summary(fselect.10)
matrix.train.10 <-matrix.train.10[fselect.10,]

##############################################################################################################
#### 4. SVM: FULL MODEL TestSet 10 ###########################################################################
##############################################################################################################

matrix.train.10 <- (t(matrix.train.10))
labels.train.10 <- as.factor(pData.train.10$Class)


#### 4.1 SetUp SVMrad for full model TestSet 10 #################################################################
##############################################################################################################

## create 200 resamples of the train data (10TestSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1920)
index.10 <- createMultiFolds(labels.train.10, k=10, times = 20)  

fullCtrl.10 <- trainControl(method = "repeatedcv",
                           repeats = 20,
                           summaryFunction = fiveStats,
                           classProbs = TRUE,
                           index = index.10,
                           allowParallel = TRUE)

set.seed(721)
svmFull.10 <- train(matrix.train.10,labels.train.10,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.10)

svmFull.10  

##############################################################################################################
#### 5. SVM-RFE TestSet 10 ###################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.10      <- rfeControl(method = "repeatedcv", repeats = 20, 
                               saveDetails = TRUE,
                               returnResamp="final", 
                               verbose = TRUE, 
                               index = index.10,
                               seeds = seeds.rfe,
                               allowParallel = TRUE)

outerctrl.10$functions         <- caretFuncs
outerctrl.10$functions$summary <- fiveStats


#### 5.10.SVM-RFE TestSet 10 #####################################################################################
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

rfe.10
write.table(rfe.10$results, file = "Results_TestSet10.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.10, type = c("g", "o"))
plot(rfe.10, type = c("g", "o"), xlim = c(0,61))

optFeatures.10 <- cbind(rfe.10$optVariables, Annotation[rfe.10$optVariables,])
write.table(optFeatures.10, file = "OptVars_rfe152_TestSet10.txt", sep="\t",col.names=NA)

#### 5.4.compare resampling performances TestSet 10 #############################################################
##############################################################################################################

rfeResamples.10 <- resamples(list("SVM_full.10" = svmFull.10,"SVM_RFE.10" = rfe.10))
sink("Resamples_rfe152_TestSet10.txt", append = TRUE)
summary(rfeResamples.10)
sink()

modelDifferences.10 <- diff(rfeResamples.10)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_TestSet10.txt", append = TRUE)
summary(modelDifferences.10)
sink()


################################################################################################
#### 6. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TestSet 10 #############################
################################################################################################

matrix.train.GA.10 <- matrix.train.10[,rfe.10$optVariables]   # subset fot the 33 optVars

#### 6.1 Parameters for outer resampling loop (to assess feature selection) ####################
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

#### 6.3 run GA  ###############################################################################
################################################################################################

system.time(GA.10<- gafs(matrix.train.GA.10, labels.train.10, 
                        iters = 40,
                        popSize = 40, pcrossover = 0.7, pmutation = 0.1, elite = 3,  
                        gafsControl = outerctrl.GA.10,
                        ## Now we pass options to `train` via "svmGA":               
                        metric = "Accuracy",
                        method = "svmRadial",
                        tuneLength = 12,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = 2,
                                                 allowParallel = FALSE)))

### 6.1.1 analyze results  ##################################################################
GA.10
GA.10$optIter
GA.10$optVariables
optVars.GA <- Annotation[GA.10$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.10_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)



##### 6  PCA on optimal variables TestSet 10 ####################################################################
##############################################################################################################
matrix.train.opt.10     <- t(matrix.train.10[,row.names(optFeatures.10)])
matrix.test.opt.10      <- matrix.test.10[row.names(optFeatures.10),]
matrix.opt.10           <- cbind(matrix.train.opt.10, matrix.test.opt.10)
pData.opt.10            <- rbind(pData.train.10,pData.test.10)
pData.opt.10$Set        <- ifelse(row.names(pData.opt.10)%in% row.names(pData.train.10),"train","test")
pData.opt.10$Design_Color <- ifelse(pData.opt.10$Set == "train", pData.opt.10$Design_Color,"#000000")

pca.10         <- prcomp(t(matrix.opt.10))           
plot(pca.10$x, pch=16, col=pData.opt.10$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","TestSet"), col = unique(pData.opt.10$Design_Color), pch=16, bty="n", cex=0.8)

##### 7  Train SVM on optVar FeatureSet.10 and predict independent TestSet.10 ##################################
##############################################################################################################

set.seed(721)
svmOpt.10  <- train(t(matrix.train.opt.10),labels.train.10,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneLength = 50,
                   trControl = fullCtrl.10)
svmOpt.10  

Prediction_optVars.10 <- predict(svmOpt.10,t(matrix.test.opt.10), type = "prob")
Prediction_optVars.10$Prediction_optVars.10 <- ifelse(Prediction_optVars.10$transforming>0.50,"transforming","untransforming")
Prediction_optVars.10 <- cbind(pData.test.10[,c(1:3)],TrueLabel=pData.test.10$Class,Prediction_optVars.10)
write.table(Prediction_optVars.10, file = paste("Predictions_TestSet 10.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


#### 8 Performance of opVars Classifier on test samples TestSet 10 #############################################
#############################################################################################################

#### 8.10 Confusion matrix  ##################################################################################  
sink("ConfusionMatrix_TestSet10.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_optVars.10$Prediction_optVars.10), as.factor(Prediction_optVars.10$TrueLabel))
sink()

#### 6.10.10 ROC on probability "transforming" TestSet 10 ########################################################

Prediction_optVars.10$Class <- as.factor(ifelse(Prediction_optVars.10$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_TestSet10.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Prediction_optVars.10$Class,                    # response vector (factor or character)
            Prediction_optVars.10$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()






#############################################################################################################################################
#############################################################################################################################################
#### II. FINAL MODEL WITHOUT INDEPENDENT TEST SET: Performance based only on resampling #####################################################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 3. nonspecific feature prefiltering FINAL  ##############################################################
##############################################################################################################

fselect.FINAL  <- genefilter(eset.batch, filterfun(f1))
summary(fselect.FINAL)
matrix.train.FINAL <-eset.batch[fselect.FINAL,]

##############################################################################################################
#### 4. SVM: FULL MODEL FINAL ################################################################################
##############################################################################################################

matrix.train.FINAL <- (t(matrix.train.FINAL))
labels.train.FINAL <- as.factor(pData$Class)


#### 4.1 SetUp SVMrad for full model FINAL #################################################################
##############################################################################################################

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

varImp_FINAL <- varImp(svmFull.FINAL)$importance
varImp_FINAL <- varImp_FINAL[order(varImp_FINAL$transforming,decreasing = T),]
varImp_FINAL <- cbind(varImp_FINAL,Annotation[row.names(varImp_FINAL),])
write.table(varImp_FINAL, file = "FINAL_VariableImportance_ROC.txt", sep="\t",col.names=NA)

##############################################################################################################
#### 5. SVM-RFE FINAL ########################################################################################
##############################################################################################################

#### 5.1 Parameters for outer resampling loop (to assess feature selection) ##################################
##############################################################################################################

outerctrl.FINAL      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                   saveDetails = TRUE,
                                   returnResamp="final", 
                                   verbose = TRUE, 
                                   index = index.FINAL,
                                   seeds = seeds.rfe,
                                   allowParallel = TRUE)

outerctrl.FINAL$functions         <- caretFuncs
outerctrl.FINAL$functions$summary <- fiveStats



#### 5.FINAL.SVM-RFE FINAL #####################################################################################
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

rfe.FINAL
write.table(rfe.FINAL$results, file = "Results_FINAL.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.FINAL, type = c("g", "o"))
plot(rfe.FINAL, type = c("g", "o"), xlim = c(0,61))

# Plot ROC over FeatureNumber
plot(rfe.FINAL$results$Variables,rfe.FINAL$results$ROC, type = "o",panel.first = grid())


optFeatures.FINAL <- cbind(rfe.FINAL$optVariables, Annotation[rfe.FINAL$optVariables,])
write.table(optFeatures.FINAL, file = "OptVars_rfe152_FINAL.txt", sep="\t",col.names=NA)

#### 5.FINAL.compare resampling performances FINAL #############################################################
##############################################################################################################

rfeResamples.FINAL <- resamples(list("SVM_full.FINAL" = svmFull.FINAL,"SVM_RFE.FINAL" = rfe.FINAL))
sink("Resamples_rfe152_FINAL.txt", append = TRUE)
summary(rfeResamples.FINAL)
sink()

modelDifferences.FINAL <- diff(rfeResamples.FINAL)  # paired t-test for H0: difference = 0 between the different models. 
sink("ModelDifferences_rfe152_FINAL.txt", append = TRUE)
summary(modelDifferences.FINAL)
sink()

