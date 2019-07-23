#############################################################################################################################################
#################################### Construction of SAGA-SVM Classifier IQR 0.8 ############################################################
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

fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...)) 

#### 1. nonspecific feature prefiltering   ###################################################################
##############################################################################################################
f1       <- function(x) (IQR(x) > 1.2)    
fselect.FINAL  <- genefilter(eset.batch, filterfun(f1))
summary(fselect.FINAL)
matrix.train.FINAL <- eset.batch[fselect.FINAL,]

#### 2. SVM: Standard model ##################################################################################
##############################################################################################################
matrix.train.FINAL <- (t(matrix.train.FINAL))
labels.train.FINAL <- as.factor(pData$Class)

## create 50 resamples of the train data (10foldCVn20) - the same index is used for SVM-rfe for comparison
set.seed(123)
index.FINAL <- createMultiFolds(labels.train.FINAL, k=10, times = 5)  

fullCtrl.FINAL <- trainControl(method = "repeatedcv",repeats = 5,
                               index = index.FINAL,
                               summaryFunction = fiveStats,
                               classProbs = TRUE,
                               allowParallel = TRUE)

set.seed(721)
svmFull.FINAL <- train(matrix.train.FINAL,labels.train.FINAL,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.FINAL)

svmFull.FINAL  # sigma = 0.003350907 and C = 64
plot(svmFull.FINAL, scales = list(x = list(log = 2)))



#### 3. SVM: fixed sigma, manual cost ########################################################################
##############################################################################################################
set.seed(123)
sigmaEstimated <- sigest(as.matrix(matrix.train.FINAL))
sigmaEstimated

svmRGridCost <- expand.grid(.sigma = sigmaEstimated[2], .C = 2^(seq(-10, 15)))

set.seed(721)
svmFull.1 <- train(matrix.train.FINAL,labels.train.FINAL,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneGrid = svmRGridCost,
                       trControl = fullCtrl.FINAL)

svmFull.1  
plot(svmFull.1, scales = list(x = list(log = 2)))


#### 4. SVM: fixed cost, manual sigma / Fig. R9 ##############################################################
##############################################################################################################

svmRGridSigma <- expand.grid(.sigma = 2^(seq(-15, 1)), .C = 64)

set.seed(721)
svmFull.2 <- train(matrix.train.FINAL,labels.train.FINAL,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneGrid = svmRGridSigma,
                   trControl = fullCtrl.FINAL)

svmFull.2  
plot(svmFull.2, scales = list(x = list(log = 2)))



#### 5. grid search / Fig. R9 ################################################################################
##############################################################################################################

svmRGrid <- expand.grid(.sigma = 2^(seq(-10, -4)), .C = 2^(seq(-2, 12)))

set.seed(721)
svmFull.3 <- train(matrix.train.FINAL,labels.train.FINAL,
                   method = "svmRadial",
                   metric = "Accuracy",
                   tuneGrid = svmRGrid,
                   trControl = fullCtrl.FINAL)

svmFull.3  
plot(svmFull.3, scales = list(x = list(log = 2)))


#### 5. centering / scaling ##################################################################################
##############################################################################################################

set.seed(721)
svmFull.FINAL.pre <- train(matrix.train.FINAL,labels.train.FINAL,
                           method = "svmRadial",
                           metric = "Accuracy",
                           tuneLength = 20,
                           preProc = c("center", "scale"),
                           trControl = fullCtrl.FINAL)
svmFull.FINAL.pre
