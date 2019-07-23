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
registerDoMC(cores = 10)    


#############################################################################################################################################
#### Read in Data / Normalized and Batch corrected log2-ExpressionSet of 152 SAGA Samples ###################################################
#############################################################################################################################################

pData      <- read.delim("SAGA_Targets_FINAL_152.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation <- read.delim("Annotation_SAGA_FINAL_KNOWN_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
eset.batch <- as.matrix(read.delim("ESET_RMA_COMBAT_KNOWN_FullSagaSet152_FINAL.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE,row.names = 1))
optVars    <- read.delim("OptVars_rfe152_FINAL.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)

#### visualize ESET
set.seed(476)  
tsne_out <- Rtsne(t(eset.batch),dims = 2, initial_dims = 20, perplexity = 11,
                  theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,
                  verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData$Design_Color, pch=16, cex=1.3)



#############################################################################################################################################
#############################################################################################################################################
#### II. Genetic Algorithm for FINAL MODEL WITHOUT INDEPENDENT TEST SET: Performance based only on resampling ###############################
#############################################################################################################################################
#############################################################################################################################################

##############################################################################################################
#### 2. Filter for optVars from SVM-RFE ######################################################################
##############################################################################################################

matrix.train.FINAL <-eset.batch[row.names(optVars),]
matrix.train.FINAL <- (t(matrix.train.FINAL))
labels.train.FINAL <- as.factor(pData$Class)

#### 2.1 SetUp SVMrad for optVars Model FINAL ################################################################
##############################################################################################################
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))

## create 50 resamples of the train data (10foldCVn5) - the same index is used for GA for comparison
set.seed(123)
index.FINAL <- createMultiFolds(labels.train.FINAL, k=10, times = 5)  

fullCtrl.FINAL <- trainControl(method = "repeatedcv",
                               repeats = 5,
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



# 152 samples, 20 predictors, 2 classes: 'transforming', 'untransforming'
# Resampling: Cross-Validated (10 fold, repeated 5 times) 
# Summary of sample sizes: 137, 137, 136, 137, 137, 136, ... 
# C          ROC        Sens       Spec       Accuracy   Kappa
# 8192.00  0.9746561  0.9238095  0.9661111  0.9479762  0.8925841

#Tuning parameter 'sigma' was held constant at a value of 0.140625
#Accuracy was used to select the optimal model using the largest value.
#The final values used for the model were sigma = 0.140625 and C = 8192.


################################################################################################
#### 6. GENETIC ALGORITHM FOR FEATURE SELECTION ################################################
################################################################################################

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

set.seed(123)
seeds.GA <- vector(mode = "integer", length = length(index.FINAL)+1)    # B+1 elements where B is the number of resamples = 101
for(i in 1:length(index.FINAL)+1) seeds.GA[[i]] <- sample.int(10000, 1)

outerctrl.GA.FINAL <- gafsControl(functions = svmGA,
                                  method = "repeatedcv", repeats = 5,
                                  index = index.FINAL,                                       
                                  seeds = seeds.GA,                                      
                                  returnResamp="all", 
                                  verbose = TRUE,
                                  maximize = c(internal = TRUE,
                                               external = TRUE),
                                 allowParallel = TRUE)                                  


#### 6.3 run GA  ###############################################################################
################################################################################################
set.seed(123)
system.time(GA.RFE_FINAL<- gafs(matrix.train.FINAL, labels.train.FINAL, 
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

### 6.1.1 analyze results  ##################################################################
GA.RFE_FINAL
GA.RFE_FINAL$optIter
GA.RFE_FINAL$optVariables
optVars.GA <- Annotation[GA.RFE_FINAL$optVariables,]
write.table(optVars.GA, file = "FINAL152_GA_optVars_CVn5S123_iter50_pop40_0.7_0.1_3_Acc.txt", sep="\t",col.names=NA)


plot(GA.RFE_FINAL) + theme_bw() + xlim(0,50)



set.seed(721)
matrix.train.GA <- matrix.train.FINAL[,GA.RFE_FINAL$optVariables]
svmFull.GA <- train(matrix.train.GA,labels.train.FINAL,
                    method = "svmRadial",
                    metric = "Accuracy",
                    tuneLength = 20,
                    trControl = fullCtrl.FINAL)

svmFull.GA


rfeResamples <- resamples(list("svmFull.FINAL" = svmFull.FINAL,"svmFull.GA" = svmFull.GA))
summary(rfeResamples)

modelDifferences <- diff(rfeResamples)  # paired t-test for H0: difference = 0 between the different models. 
summary(modelDifferences)






