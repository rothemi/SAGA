#############################################################################################################################################
#################################### FeatureSelection for complete SAGA dataset of 152 samples  #############################################
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
library(PRROC)
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
#### FINAL MODEL USED IN SAGA R PACKAGE: ALL 152 SAMPLES WITHOUT INDEPENDENT TEST SET: Performance estimated only from resampling ###########
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 1. nonspecific feature prefiltering FINAL  ##############################################################
##############################################################################################################
f1       <- function(x) (IQR(x) > 0.8) 
fselect.FINAL  <- genefilter(eset.batch, filterfun(f1))
summary(fselect.FINAL)                        # 1243 features selected by nonspecific filtering: Supplementary Data 6 tab 2
matrix.train.FINAL <-eset.batch[fselect.FINAL,]

Annotation.matrix.train.FINAL <- Annotation[row.names(matrix.train.FINAL),]
write.table(Annotation.matrix.train.FINAL, file = "Annotation.matrix.train.FINAL.txt", sep="\t",col.names=NA)

##############################################################################################################
#### 2. SVM: FULL MODEL FINAL ################################################################################
##############################################################################################################
matrix.train.FINAL <- (t(matrix.train.FINAL))
labels.train.FINAL <- as.factor(pData$Class)

# calculate performance measures (accuracy, sensitivity, specificity, ROC) of external resamples
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))   


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

##  set predictors subset sizes to test: 1,2,…,40,45,50,60,….,500,1151  = 52 subsets in total
FeatureNumbers <- c(seq(1,40,by=1),45,50,60,70,80,90,100,200,300,400,500)                 

## set all seeds for reproducibility: 52 seeds for each of the 200 resamples + 1 for the complete set 
set.seed(123)
seeds.rfe <- vector(mode = "list", length = length(index.FINAL)+1)                            
for(i in 1:length(index.FINAL)) seeds.rfe[[i]] <- sample.int(10000, length(FeatureNumbers)+1) 
seeds.rfe[[length(index.FINAL)+1]] <- sample.int(10000, 1)

outerctrl.FINAL      <- rfeControl(method = "repeatedcv", repeats = 20, 
                                   index = index.FINAL,
                                   saveDetails = TRUE,
                                   returnResamp="final", 
                                   verbose = TRUE, 
                                   seeds = seeds.rfe,
                                   allowParallel = TRUE)

outerctrl.FINAL$functions         <- caretFuncs
outerctrl.FINAL$functions$summary <- fiveStats

#### 3.2 Parameters for inner (nested) cross-validation loop for hyperparameter tuning within each resample and for each subset size #######
##############################################################################################################
innerctrl <- trainControl(method = "repeatedcv",repeats = 3,           # 10CVn3 = 30 resamples
                          verboseIter = FALSE,
                          classProbs = TRUE,
                          allowParallel = FALSE)                      


#### 3.2 SVM-RFE FINAL #######################################################################################
##############################################################################################################

# This will take ca. 3-4 hrs on AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM
system.time(rfe.FINAL  <- rfe(matrix.train.FINAL, labels.train.FINAL, 
                              sizes=FeatureNumbers,
                              rfeControl=outerctrl.FINAL,
                              metric = "Accuracy",
                              method="svmRadial",
                              tuneLength = 20,
                              trControl = innerctrl))

rfe.FINAL   # 20 optVars found by rfe
write.table(rfe.FINAL$results, file = "FinalSet152_Results_rfe.txt", sep="\t",col.names=NA)

# Figure 5a
trellis.par.set(caretTheme())
plot(rfe.FINAL, type = c("g", "o"))  

# Plot ROC over FeatureNumber
plot(rfe.FINAL$results$Variables,rfe.FINAL$results$ROC, type = "o",panel.first = grid())

optFeatures.FINAL <- cbind(rfe.FINAL$optVariables, Annotation[rfe.FINAL$optVariables,]) #Supplementary Data 6 tab 3
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

# This will take ca. 4 hrs on AWS EC2 c5.18xlarge instance with 72 CPUs and 144 Gb of RAM
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
GA.FINAL  # Top11 Features selected by GA: Supplementary Data 6 tab 4 / used as final predictor set for SAGA_SVM and SAGA_GSEA
optVars.GA.FINAL <- Annotation[GA.FINAL$optVariables,]  
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

### Figure 5b plot external accuracy over the iterations ################################################################
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

### PCA on all 36,226 probes
pca.eset.batch <- prcomp(t(eset.batch))           
plot(pca.eset.batch$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(35,55, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

### Figure 5c: PCA on 11 probes from SVM-GA
matrix.train.GA.FINAL  <- eset.batch[GA.FINAL$optVariables,]
pca.FINAL              <- prcomp(t(matrix.train.GA.FINAL))           
plot(pca.FINAL$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(2,4, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

### Figure 5d: PCA on 11 random probes 
index.5d      <- sample(36226,11, replace = FALSE)
matrix.random <- eset.batch[index.5d,]
pca.random    <- prcomp(t(matrix.random))           
plot(pca.random$x, pch=16, col=pData$Design_Color, cex=1.5, asp=1)
legend(2,5, legend = c("transforming","mock","neutral"), col = unique(pData$Design_Color), pch=16, bty="n", cex=1)

################################################################################################################
### 6. ROC curve for outer resamples of SVM.FINAL / SVM.rfe.FINAL: #############################################
################################################################################################################

### 6.1 ROC curve on external SVM.rfe resamples 
a <- rfe.FINAL$pred
b <- subset(a, a$Variables == 20)  # only for the best feature subset

b$Class <- as.factor(ifelse(b$obs == "transforming","transforming","nontransforming"))
roc.b <- roc(b$Class,                    
             b$transforming,             
             percent=TRUE, levels=c("nontransforming","transforming"),
             plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
             print.auc=T,print.thres=0.5)

roc.b  # AUC = 95.3% sensitivity 88.5% specificity 93.4%


#############################################################################################################################################
#############################################################################################################################################
#### Analyze Resampling Accuracy (external resamples): Supplementary Data 6 tab 1 ###########################################################
#############################################################################################################################################
#############################################################################################################################################

library(gmodels)

### 2.1. create data.frame for all 10 TestSet/TrainingSet Splits
Results <- data.frame(Testset = "FinalSet152", FeatureNumber_full = rep(NA,1), ResamplingAccuracy_full=rep(NA,1),CI_full_lower=rep(NA,1),CI_full_upper=rep(NA,1), 
                      FeatureNumber_rfe = rep(NA,1), ResamplingAccuracy_rfe = rep(NA,1), CI_rfe_lower=rep(NA,1),CI_rfe_upper=rep(NA,1), P_full_vs_rfe=rep(NA,1),
                      FeatureNumber_GA  = rep(NA,1), ResamplingAccuracy_GA = rep(NA,1),CI_GA_lower=rep(NA,1),CI_GA_upper=rep(NA,1),
                      Accuracy_TestSet_full = rep(NA,1), CI_Test_full_lower=rep(NA,1),CI_Test_full_upper=rep(NA,1),
                      Accuracy_TestSet_rfe = rep(NA,1),  CI_Test_rfe_lower=rep(NA,1),CI_Test_rfe_upper=rep(NA,1),
                      Accuracy_TestSet_GA = rep(NA,1),   CI_Test_GA_lower=rep(NA,1),CI_Test_GA_upper=rep(NA,1),
                      Accuracy_TestSet_SAGA = rep(NA,1),AUROC_SAGA = rep(NA,1), AUPRC_SAGA = rep(NA,1),AUPRC_random = rep(NA,1))
                      
### 2.2. number of features selected by IQR filter for each training/test split
Results[1,2] <- dim(matrix.train.FINAL)[2]

### 2.3. resampling accuracy for the full model / stored in svmFull.x$results dataframe
Results[1,3] <- svmFull.FINAL$results[svmFull.FINAL$results$C == svmFull.FINAL$bestTune$C,"Accuracy"]

### 2.4. lower bound of confidence intervall for the resampling accuracy of full model / stored in rfeResamples.x$values dataframe
Results[1,4] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[2] 

### 2.5. upper bound of confidence intervall for the resampling accuracy of full model / stored in rfeResamples.x$values dataframe
Results[1,5] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[3] 

### 2.6. number of features found by SVM-rfe for each training/test split
Results[1,6] <- rfe.FINAL$optsize

### 2.7. external resampling accuracy for the svm-rfe model
Results[1,7] <- rfe.FINAL$results[rfe.FINAL$results$Variables == rfe.FINAL$optsize,5]

### 2.8. fill in lower bound confidence intervall for the resampling accuracy of rfe model / stored in rfeResamples.x$values dataframe
Results[1,8] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[2] 

### 2.9. upper bound confidence intervall for the resampling accuracy of rfe model / stored in rfeResamples.x$values dataframe
Results[1,9] <- ci(rfeResamples.FINAL$values$`SVM_full.FINAL~Accuracy`)[3] 

### 2.10. p.value for H0 = resampling accuracy of full model = resampling accuracy of full model/ stored in rfeResamples.x$values dataframe
Results[1,10] <- modelDifferences.FINAL$statistics$Accuracy[[1]][3]

### 2.11. number of features found by SVM-GA 
Results[1,11] <- length(GA.FINAL$optVariables)

### 2.12. Resampling accuracy for the GA model
Results[1,12] <- accuracy.external.FINAL

### 2.13. calculate lower bound confidence intervall for the resampling accuracy of GA model / stored in rfeResamples.x$values dataframe
Results[1,13] <- ci(accuracy.external.opt.FINAL)[2] 

### 2.14. calculate upper bound confidence intervall for the resampling accuracy of GA model / stored in rfeResamples.x$values dataframe
Results[1,14] <- ci(accuracy.external.opt.FINAL)[3] 

write.table(Results, file = "Results_resampling_rfe_GA40_FinalSet152.txt", sep="\t",col.names=NA)

