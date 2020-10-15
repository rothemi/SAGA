library(limma)
library(phenoTest)
library(gridExtra)
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
rm(SAGA_RAW.full)
rm(RAW.048306,RAW.048306_Agilent,RAW.066423,RAW.066423_Agilent,RAW.084107,RAW.084107_Agilent,RAW.084956,RAW.084956_Agilent)
rm(files.048306,files.066423,files.084107,files.084956)
rm(index1,index2,index3,index4)


#############################################################################################################################################
#############################################################################################################################################
#### Figure 6a: 3 IVIMs ####################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
pData.3IVIMs  <- pData[c(1:21),]  # select batch 1-3 
RAW_3IVIMs    <- SAGA_RAW[,row.names(pData.3IVIMs)]
boxplot(log2(RAW_3IVIMs$E),col=pData.3IVIMs$IVIM_Color,names=pData.3IVIMs$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

## quantile normalization 
RMA_3IVIMs <- normalizeBetweenArrays(RAW_3IVIMs,method="quantile")      # limma: EList object with quantile normalized and log2 data
RMA_3IVIMs <- avereps(RMA_3IVIMs, ID= RMA_3IVIMs$genes$ProbeName)      # limma: EList object averaged over ProbeIDs    
boxplot(RMA_3IVIMs$E, col=pData.3IVIMs$IVIM_Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)    

## t-SNE
set.seed(40)  
tsne_out <- Rtsne(t(RMA_3IVIMs$E),dims = 2, initial_dims = 20, perplexity = 4,
                  theta = 0.5, check_duplicates = TRUE, pca = TRUE, max_iter = 1000,
                  verbose = FALSE, is_distance = FALSE)

## Figure 6a
pData.3IVIMs  <- cbind(as.data.frame(tsne_out$Y),pData.3IVIMs)
ggplot(pData.3IVIMs,aes(x=V1,y=V2, colour = as.factor(Vector), shape = Immortalized)) +
  geom_point(size=5) + 
  theme_bw(base_size = 12, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_colour_manual(values = unique(pData$Vector_Color)) +
  scale_y_reverse()


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
boxplot(RMA_train.1$E, col=pData.train.1$IVIM_Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.1$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.1 <- pData.train.1$Batch-1                       
modcombat     <- model.matrix(~1, data=pData.train.1)         
matrix.train.1.batch <- ComBat(dat=RMA_train.1$E, batch=batch.train.1, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.1.batch <- matrix.train.1.batch[row.names(Annotation.known),] # filter for 36,226 annotated probes

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.1.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.1$Design_Color, pch=16, cex=1.3) 

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

matrix.train.rfe.1 <- matrix.train.1[,rfe.1$optVariables]   # subset fot the 20 optVars from rfe

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

### 6.2 analyze results from Genetic Algorithm #################################################
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

################################################################################################
#### 7. SAGA-GSEA for TestSet 01 ###############################################################
################################################################################################
  SAGA.CORE.1 <- list(SAGA.RFE=rfe.1$optVariables, SAGA.GA=GA.1$optVariables)
  #### 8.1. Normalize, average ###################################################################
  RMA.1 <- normalizeBetweenArrays(RAW_test.1, method="quantile")       # quantil normalize TestSet 01
  RMA.1 <- avereps(RMA.1,ID= RMA.1$genes$ProbeName)                    # average replicates to one value for each probe
  matrix.gsea.1 <- RMA.1$E                                             # extract log2 expression values 
  
  #### 8.2. make ExpressionSet (Biobase) object ##################################################
  metadata.1  <- data.frame(labelDescription= rep(NA,dim(pData.test.1)[2]),row.names=colnames(pData.test.1))   # varMetadata: empty, but required 
  phenoData.1 <- new("AnnotatedDataFrame",data=pData.test.1, varMetadata=metadata.1)   # annotatedDataFrame for the annotation of the samples
  eset.gsea.1 <- ExpressionSet(assayData = matrix.gsea.1, phenoData = phenoData.1)     # this is the ExpressionSet required for phenoTest
  
  #### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
  vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
  epheno.gsea.1 <- ExpressionPhenoTest(eset.gsea.1,vars2test,p.adjust.method='BH')
  
  #### 8.4 GSEA #################################################################################
  SAGA.GSEA.1 <- gsea(x=epheno.gsea.1, gsets=SAGA.CORE.1 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                      center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
  #### 8.5 output ###############################################################################
  Output.GSEA.1 <- summary(SAGA.GSEA.1)[,c(1,2,3,5,8)]
  GSEA.RFE.1    <- subset(Output.GSEA.1,Output.GSEA.1$geneSet == "SAGA.RFE")                
  colnames(GSEA.RFE.1) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
  GSEA.GA.1     <- subset(Output.GSEA.1,Output.GSEA.1$geneSet == "SAGA.GA")
  colnames(GSEA.GA.1) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
  Output.GSEA.1 <- cbind(GSEA.RFE.1,GSEA.GA.1[,-1])
  
  Vector <- NULL    ### pull out the Vector index number from the result table                  
  for (a in 1:nrow(Output.GSEA.1)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.1$variable[a]), ".", fixed = TRUE))[2] }
  Output.GSEA.1$GSEA_Vector <- Vector
  
  pData.Test.sub.1           <- pData.test.1[pData.test.1$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
  pData.Test.sub.1$SampleID  <- row.names(pData.Test.sub.1)                   
  GSEA.result.1              <- merge(pData.Test.sub.1,Output.GSEA.1, by.x="GSEA_Vector", by.y = "GSEA_Vector") 
  
  # make pdf report: 
  pdf(file="SAGA.GSEA_Batch_1.pdf",useDingbats = F,width = 10, height = 10)  
  grid.table(summary(SAGA.GSEA.1),rows = NULL)
  plot(SAGA.GSEA.1,es.nes='nes')
  dev.off()


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

#### 2.2 quantile normalization  ############################################################################# 
RMA_train.2 <- normalizeBetweenArrays(RAW_train.2,method="quantile")      # quantile normalization of training samples
RMA_train.2 <- avereps(RMA_train.2, ID= RMA_train.2$genes$ProbeName)      # average over ProbeIDs of training samples    
boxplot(RMA_train.2$E, col=pData.train.2$IVIM_Color,boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)    

#### 2.2 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.2$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.2$IVIM_Color, pch=16, cex=1.3) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.2 <- ifelse(pData.train.2$Batch>2,pData.train.2$Batch-1,pData.train.2$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.2)         
matrix.train.2.batch <- ComBat(dat=RMA_train.2$E, batch=batch.train.2, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.2.batch <- matrix.train.2.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.2.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.2$Design_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 02 ######################################################
##############################################################################################################

fselect.2  <- genefilter(matrix.train.2.batch, filterfun(f1))
summary(fselect.2)  # 1355 features 
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

rfe.2  # 9 optimal predictors found ==> no GA necessary
write.table(rfe.2$results, file = "TrainingSet02_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.2, type = c("g", "o"))
plot(rfe.2, type = c("g", "o"), xlim = c(0,62))

optFeatures.rfe.2 <- Annotation[rfe.2$optVariables,]
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
#### 8. SAGA-GSEA for TestSet 02 ###############################################################
################################################################################################
set.2       <- data.frame(as.list(rfe.2$optVariables), row.names = "SAGA.RFE", stringsAsFactors = FALSE)   # make geneset from SVM-RFE optimal variables
SAGA.CORE.2 <- setNames(split(set.2, seq(nrow(set.2))), rownames(set.2))   

#### 8.1. Normalize, average ###################################################################
RMA.2 <- normalizeBetweenArrays(RAW_test.2, method="quantile")       # quantil normalize TestSet 01
RMA.2 <- avereps(RMA.2,ID= RMA.2$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.2 <- RMA.2$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.2  <- data.frame(labelDescription= rep(NA,dim(pData.test.2)[2]),row.names=colnames(pData.test.2))   # varMetadata: empty, but required 
phenoData.2 <- new("AnnotatedDataFrame",data=pData.test.2, varMetadata=metadata.2)   # annotatedDataFrame for the annotation of the samples
eset.gsea.2 <- ExpressionSet(assayData = matrix.gsea.2, phenoData = phenoData.2)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.2 <- ExpressionPhenoTest(eset.gsea.2,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.2 <- gsea(x=epheno.gsea.2, gsets=SAGA.CORE.2 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes=8)
Result.2    <- summary(SAGA.GSEA.2)[,c(1,2,3,5,8)]              # extract results (only NES- normalized enrichment scores)

#### 8.5 output ###############################################################################
Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Result.2)) {Vector[a] <- unlist(strsplit(as.character(Result.2$variable[a]), ".", fixed = TRUE))[2] }
Result.2$GSEA_Vector <- Vector

pData.Test.sub.2           <- pData.test.2[pData.test.2$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.2$SampleID  <- row.names(pData.Test.sub.2)                   
GSEA.result.2              <- merge(pData.Test.sub.2,Result.2, by.x="GSEA_Vector", by.y = "GSEA_Vector") 
GSEA.result.2$geneSet.GA   <- rep(NA,7)
GSEA.result.2$n.GA         <- rep(NA,7)
GSEA.result.2$nes.GA       <- rep(NA,7)
GSEA.result.2$fdr.GA       <- rep(NA,7)
colnames(GSEA.result.2)    <- colnames(GSEA.result.1)

# make pdf report:
pdf(file="SAGA.GSEA_Batch_2.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.2),rows = NULL)
plot(SAGA.GSEA.2,es.nes='nes',selGsets='SAGA.CORE')
dev.off()




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
boxplot(log2(RAW_train.3$E), col=pData.train.3$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.3 quantile normalization  ############################################################################# 
RMA_train.3 <- normalizeBetweenArrays(RAW_train.3,method="quantile")      # quantile normalization
RMA_train.3 <- avereps(RMA_train.3, ID= RMA_train.3$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.3$E, col=pData.train.3$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.3 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.3$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.3$IVIM_Color, pch=16, cex=1.4) 

#### 2.3 COMBAT batch correction ############################################################################# 
batch.train.3 <- ifelse(pData.train.3$Batch>3,pData.train.3$Batch-1,pData.train.3$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.3)         
matrix.train.3.batch <- ComBat(dat=RMA_train.3$E, batch=batch.train.3, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.3.batch <- matrix.train.3.batch[row.names(Annotation.known),]

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

rfe.3    # 22 predictors found 
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


################################################################################################
#### 8. SAGA-GSEA for TestSet 03 ###############################################################
################################################################################################
SAGA.CORE.3 <- list(SAGA.RFE=rfe.3$optVariables, SAGA.GA=GA.3$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.3 <- normalizeBetweenArrays(RAW_test.3, method="quantile")       # quantil normalize TestSet .3
RMA.3 <- avereps(RMA.3,ID= RMA.3$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.3 <- RMA.3$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.3  <- data.frame(labelDescription= rep(NA,dim(pData.test.3)[2]),row.names=colnames(pData.test.3))   # varMetadata: empty, but required 
phenoData.3 <- new("AnnotatedDataFrame",data=pData.test.3, varMetadata=metadata.3)   # annotatedDataFrame for the annotation of the samples
eset.gsea.3 <- ExpressionSet(assayData = matrix.gsea.3, phenoData = phenoData.3)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.3 <- ExpressionPhenoTest(eset.gsea.3,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.3 <- gsea(x=epheno.gsea.3, gsets=SAGA.CORE.3 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.3 <- summary(SAGA.GSEA.3)[,c(1,2,3,5,8)]
GSEA.RFE.3    <- subset(Output.GSEA.3,Output.GSEA.3$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.3) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.3     <- subset(Output.GSEA.3,Output.GSEA.3$geneSet == "SAGA.GA")
colnames(GSEA.GA.3) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.3 <- cbind(GSEA.RFE.3,GSEA.GA.3[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.3)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.3$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.3$GSEA_Vector <- Vector

pData.Test.sub.3           <- pData.test.3[pData.test.3$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.3$SampleID  <- row.names(pData.Test.sub.3)                   
GSEA.result.3              <- merge(pData.Test.sub.3,Output.GSEA.3, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_3.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.3),rows = NULL)
plot(SAGA.GSEA.3,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.4$E), col=pData.train.4$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.4 quantile normalization  ############################################################################# 
RMA_train.4 <- normalizeBetweenArrays(RAW_train.4,method="quantile")      # quantile normalization
RMA_train.4 <- avereps(RMA_train.4, ID= RMA_train.4$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.4$E, col=pData.train.4$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.4 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.4$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.4$IVIM_Color, pch=16, cex=1.4) 

#### 2.4 COMBAT batch correction ############################################################################# 
batch.train.4 <- ifelse(pData.train.4$Batch>4,pData.train.4$Batch-1,pData.train.4$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.4)         
matrix.train.4.batch <- ComBat(dat=RMA_train.4$E, batch=batch.train.4, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.4.batch <- matrix.train.4.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.4.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.4$IVIM_Color, pch=16, cex=1.4) 

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
GA.4   # yields 7 features
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 04 ###############################################################
################################################################################################
SAGA.CORE.4 <- list(SAGA.RFE=rfe.4$optVariables, SAGA.GA=GA.4$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.4 <- normalizeBetweenArrays(RAW_test.4, method="quantile")       # quantil normalize TestSet .4
RMA.4 <- avereps(RMA.4,ID= RMA.4$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.4 <- RMA.4$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.4  <- data.frame(labelDescription= rep(NA,dim(pData.test.4)[2]),row.names=colnames(pData.test.4))   # varMetadata: empty, but required 
phenoData.4 <- new("AnnotatedDataFrame",data=pData.test.4, varMetadata=metadata.4)   # annotatedDataFrame for the annotation of the samples
eset.gsea.4 <- ExpressionSet(assayData = matrix.gsea.4, phenoData = phenoData.4)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.4 <- ExpressionPhenoTest(eset.gsea.4,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.4 <- gsea(x=epheno.gsea.4, gsets=SAGA.CORE.4 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.4 <- summary(SAGA.GSEA.4)[,c(1,2,3,5,8)]
GSEA.RFE.4    <- subset(Output.GSEA.4,Output.GSEA.4$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.4) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.4     <- subset(Output.GSEA.4,Output.GSEA.4$geneSet == "SAGA.GA")
colnames(GSEA.GA.4) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.4 <- cbind(GSEA.RFE.4,GSEA.GA.4[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.4)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.4$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.4$GSEA_Vector <- Vector

pData.Test.sub.4           <- pData.test.4[pData.test.4$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.4$SampleID  <- row.names(pData.Test.sub.4)                   
GSEA.result.4              <- merge(pData.Test.sub.4,Output.GSEA.4, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_4.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.4),rows = NULL)
plot(SAGA.GSEA.4,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.5$E), col=pData.train.5$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.5 quantile normalization  ############################################################################# 
RMA_train.5 <- normalizeBetweenArrays(RAW_train.5,method="quantile")      # quantile normalization
RMA_train.5 <- avereps(RMA_train.5, ID= RMA_train.5$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.5$E, col=pData.train.5$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.5 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.5$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.5$IVIM_Color, pch=16, cex=1.4) 

#### 2.5 COMBAT batch correction ############################################################################# 
batch.train.5 <- ifelse(pData.train.5$Batch>5,pData.train.5$Batch-1,pData.train.5$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.5)         
matrix.train.5.batch <- ComBat(dat=RMA_train.5$E, batch=batch.train.5, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.5.batch <- matrix.train.5.batch[row.names(Annotation.known),]
#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.5.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.5$IVIM_Color, pch=16, cex=1.3) 

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
GA.5   # yields 23 features
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


################################################################################################
#### 8. SAGA-GSEA for TestSet 05 ###############################################################
################################################################################################
SAGA.CORE.5 <- list(SAGA.RFE=rfe.5$optVariables, SAGA.GA=GA.5$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.5 <- normalizeBetweenArrays(RAW_test.5, method="quantile")       # quantil normalize TestSet .5
RMA.5 <- avereps(RMA.5,ID= RMA.5$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.5 <- RMA.5$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.5  <- data.frame(labelDescription= rep(NA,dim(pData.test.5)[2]),row.names=colnames(pData.test.5))   # varMetadata: empty, but required 
phenoData.5 <- new("AnnotatedDataFrame",data=pData.test.5, varMetadata=metadata.5)   # annotatedDataFrame for the annotation of the samples
eset.gsea.5 <- ExpressionSet(assayData = matrix.gsea.5, phenoData = phenoData.5)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.5 <- ExpressionPhenoTest(eset.gsea.5,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.5 <- gsea(x=epheno.gsea.5, gsets=SAGA.CORE.5 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.5 <- summary(SAGA.GSEA.5)[,c(1,2,3,5,8)]
GSEA.RFE.5    <- subset(Output.GSEA.5,Output.GSEA.5$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.5) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.5     <- subset(Output.GSEA.5,Output.GSEA.5$geneSet == "SAGA.GA")
colnames(GSEA.GA.5) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.5 <- cbind(GSEA.RFE.5,GSEA.GA.5[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.5)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.5$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.5$GSEA_Vector <- Vector

pData.Test.sub.5           <- pData.test.5[pData.test.5$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.5$SampleID  <- row.names(pData.Test.sub.5)                   
GSEA.result.5              <- merge(pData.Test.sub.5,Output.GSEA.5, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_5.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.5),rows = NULL)
plot(SAGA.GSEA.5,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.6$E), col=pData.train.6$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.6 quantile normalization  ############################################################################# 
RMA_train.6 <- normalizeBetweenArrays(RAW_train.6,method="quantile")      # quantile normalization
RMA_train.6 <- avereps(RMA_train.6, ID= RMA_train.6$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.6$E, col=pData.train.6$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.6 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.6$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.6$IVIM_Color, pch=16, cex=1.4) 

#### 2.6 COMBAT batch correction ############################################################################# 
batch.train.6 <- ifelse(pData.train.6$Batch>6,pData.train.6$Batch-1,pData.train.6$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.6)         
matrix.train.6.batch <- ComBat(dat=RMA_train.6$E, batch=batch.train.6, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.6.batch <- matrix.train.6.batch[row.names(Annotation.known),]
#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.6.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.6$IVIM_Color, pch=16, cex=1.3) 

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

rfe.6  # 31 optVars found 
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
GA.6   # selectes 13 features in itation 4
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 06 ###############################################################
################################################################################################
SAGA.CORE.6 <- list(SAGA.RFE=rfe.6$optVariables, SAGA.GA=GA.6$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.6 <- normalizeBetweenArrays(RAW_test.6, method="quantile")       # quantil normalize TestSet .6
RMA.6 <- avereps(RMA.6,ID= RMA.6$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.6 <- RMA.6$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.6  <- data.frame(labelDescription= rep(NA,dim(pData.test.6)[2]),row.names=colnames(pData.test.6))   # varMetadata: empty, but required 
phenoData.6 <- new("AnnotatedDataFrame",data=pData.test.6, varMetadata=metadata.6)   # annotatedDataFrame for the annotation of the samples
eset.gsea.6 <- ExpressionSet(assayData = matrix.gsea.6, phenoData = phenoData.6)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.6 <- ExpressionPhenoTest(eset.gsea.6,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.6 <- gsea(x=epheno.gsea.6, gsets=SAGA.CORE.6 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.6 <- summary(SAGA.GSEA.6)[,c(1,2,3,5,8)]
GSEA.RFE.6    <- subset(Output.GSEA.6,Output.GSEA.6$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.6) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.6     <- subset(Output.GSEA.6,Output.GSEA.6$geneSet == "SAGA.GA")
colnames(GSEA.GA.6) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.6 <- cbind(GSEA.RFE.6,GSEA.GA.6[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.6)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.6$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.6$GSEA_Vector <- Vector

pData.Test.sub.6           <- pData.test.6[pData.test.6$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.6$SampleID  <- row.names(pData.Test.sub.6)                   
GSEA.result.6              <- merge(pData.Test.sub.6,Output.GSEA.6, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_6.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.6),rows = NULL)
plot(SAGA.GSEA.6,es.nes='nes')
dev.off()



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
boxplot(log2(RAW_train.7$E), col=pData.train.7$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.7 quantile normalization  ############################################################################# 
RMA_train.7 <- normalizeBetweenArrays(RAW_train.7,method="quantile")      # quantile normalization
RMA_train.7 <- avereps(RMA_train.7, ID= RMA_train.7$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.7$E, col=pData.train.7$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.7 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.7$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.7$IVIM_Color, pch=16, cex=1.4) 

#### 2.7 COMBAT batch correction ############################################################################# 
batch.train.7 <- ifelse(pData.train.7$Batch>7,pData.train.7$Batch-1,pData.train.7$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.7)         
matrix.train.7.batch <- ComBat(dat=RMA_train.7$E, batch=batch.train.7, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.7.batch <- matrix.train.7.batch[row.names(Annotation.known),]
#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.7.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.7$IVIM_Color, pch=16, cex=1.3) 

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

rfe.7    # 40 optVars found by rfe
write.table(rfe.7$results, file = "TrainingSet07_Results_rfe.txt", sep="\t",col.names=NA)

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
GA.7    # 8 optVars found by GA
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


################################################################################################
#### 8. SAGA-GSEA for TestSet 07 ###############################################################
################################################################################################
SAGA.CORE.7 <- list(SAGA.RFE=rfe.7$optVariables, SAGA.GA=GA.7$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.7 <- normalizeBetweenArrays(RAW_test.7, method="quantile")       # quantil normalize TestSet .7
RMA.7 <- avereps(RMA.7,ID= RMA.7$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.7 <- RMA.7$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.7  <- data.frame(labelDescription= rep(NA,dim(pData.test.7)[2]),row.names=colnames(pData.test.7))   # varMetadata: empty, but required 
phenoData.7 <- new("AnnotatedDataFrame",data=pData.test.7, varMetadata=metadata.7)   # annotatedDataFrame for the annotation of the samples
eset.gsea.7 <- ExpressionSet(assayData = matrix.gsea.7, phenoData = phenoData.7)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.7 <- ExpressionPhenoTest(eset.gsea.7,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.7 <- gsea(x=epheno.gsea.7, gsets=SAGA.CORE.7 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.7 <- summary(SAGA.GSEA.7)[,c(1,2,3,5,8)]
GSEA.RFE.7    <- subset(Output.GSEA.7,Output.GSEA.7$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.7) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.7     <- subset(Output.GSEA.7,Output.GSEA.7$geneSet == "SAGA.GA")
colnames(GSEA.GA.7) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.7 <- cbind(GSEA.RFE.7,GSEA.GA.7[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.7)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.7$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.7$GSEA_Vector <- Vector

pData.Test.sub.7           <- pData.test.7[pData.test.7$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.7$SampleID  <- row.names(pData.Test.sub.7)                   
GSEA.result.7              <- merge(pData.Test.sub.7,Output.GSEA.7, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_7.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.7),rows = NULL)
plot(SAGA.GSEA.7,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.8$E), col=pData.train.8$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.8 quantile normalization  ############################################################################# 
RMA_train.8 <- normalizeBetweenArrays(RAW_train.8,method="quantile")      # quantile normalization
RMA_train.8 <- avereps(RMA_train.8, ID= RMA_train.8$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.8$E, col=pData.train.8$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.8 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.8$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.8$IVIM_Color, pch=16, cex=1.4) 

#### 2.8 COMBAT batch correction ############################################################################# 
batch.train.8 <- ifelse(pData.train.8$Batch>8,pData.train.8$Batch-1,pData.train.8$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.8)         
matrix.train.8.batch <- ComBat(dat=RMA_train.8$E, batch=batch.train.8, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.8.batch <- matrix.train.8.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.8.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.8$IVIM_Color, pch=16, cex=1.3) 

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
GA.8   # yields 13 features
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 08 ###############################################################
################################################################################################
SAGA.CORE.8 <- list(SAGA.RFE=rfe.8$optVariables, SAGA.GA=GA.8$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.8 <- normalizeBetweenArrays(RAW_test.8, method="quantile")       # quantil normalize TestSet .8
RMA.8 <- avereps(RMA.8,ID= RMA.8$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.8 <- RMA.8$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.8  <- data.frame(labelDescription= rep(NA,dim(pData.test.8)[2]),row.names=colnames(pData.test.8))   # varMetadata: empty, but required 
phenoData.8 <- new("AnnotatedDataFrame",data=pData.test.8, varMetadata=metadata.8)   # annotatedDataFrame for the annotation of the samples
eset.gsea.8 <- ExpressionSet(assayData = matrix.gsea.8, phenoData = phenoData.8)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.8 <- ExpressionPhenoTest(eset.gsea.8,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.8 <- gsea(x=epheno.gsea.8, gsets=SAGA.CORE.8 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.8 <- summary(SAGA.GSEA.8)[,c(1,2,3,5,8)]
GSEA.RFE.8    <- subset(Output.GSEA.8,Output.GSEA.8$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.8) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.8     <- subset(Output.GSEA.8,Output.GSEA.8$geneSet == "SAGA.GA")
colnames(GSEA.GA.8) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.8 <- cbind(GSEA.RFE.8,GSEA.GA.8[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.8)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.8$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.8$GSEA_Vector <- Vector

pData.Test.sub.8           <- pData.test.8[pData.test.8$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.8$SampleID  <- row.names(pData.Test.sub.8)                   
GSEA.result.8              <- merge(pData.Test.sub.8,Output.GSEA.8, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_8.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.8),rows = NULL)
plot(SAGA.GSEA.8,es.nes='nes')
dev.off()



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
boxplot(log2(RAW_train.9$E), col=pData.train.9$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.9 quantile normalization  ############################################################################# 
RMA_train.9 <- normalizeBetweenArrays(RAW_train.9,method="quantile")      # quantile normalization
RMA_train.9 <- avereps(RMA_train.9, ID= RMA_train.9$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.9$E, col=pData.train.9$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.9 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.9$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.9$IVIM_Color, pch=16, cex=1.4) 

#### 2.9 COMBAT batch correction ############################################################################# 
batch.train.9 <- ifelse(pData.train.9$Batch>9,pData.train.9$Batch-1,pData.train.9$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.9)         
matrix.train.9.batch <- ComBat(dat=RMA_train.9$E, batch=batch.train.9, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.9.batch <- matrix.train.9.batch[row.names(Annotation.known),]
matrix.train.9.batch <- round(matrix.train.9.batch,4)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.9.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.9$IVIM_Color, pch=16, cex=1.3) 

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
GA.9   # 20 features selected at iteration 4
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 09 ###############################################################
################################################################################################
SAGA.CORE.9 <- list(SAGA.RFE=rfe.9$optVariables, SAGA.GA=GA.9$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.9 <- normalizeBetweenArrays(RAW_test.9, method="quantile")       # quantil normalize TestSet .9
RMA.9 <- avereps(RMA.9,ID= RMA.9$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.9 <- RMA.9$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.9  <- data.frame(labelDescription= rep(NA,dim(pData.test.9)[2]),row.names=colnames(pData.test.9))   # varMetadata: empty, but required 
phenoData.9 <- new("AnnotatedDataFrame",data=pData.test.9, varMetadata=metadata.9)   # annotatedDataFrame for the annotation of the samples
eset.gsea.9 <- ExpressionSet(assayData = matrix.gsea.9, phenoData = phenoData.9)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.9 <- ExpressionPhenoTest(eset.gsea.9,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.9 <- gsea(x=epheno.gsea.9, gsets=SAGA.CORE.9 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.9 <- summary(SAGA.GSEA.9)[,c(1,2,3,5,8)]
GSEA.RFE.9    <- subset(Output.GSEA.9,Output.GSEA.9$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.9) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.9     <- subset(Output.GSEA.9,Output.GSEA.9$geneSet == "SAGA.GA")
colnames(GSEA.GA.9) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.9 <- cbind(GSEA.RFE.9,GSEA.GA.9[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.9)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.9$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.9$GSEA_Vector <- Vector

pData.Test.sub.9           <- pData.test.9[pData.test.9$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.9$SampleID  <- row.names(pData.Test.sub.9)                   
GSEA.result.9              <- merge(pData.Test.sub.9,Output.GSEA.9, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_9.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.9),rows = NULL)
plot(SAGA.GSEA.9,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.10$E), col=pData.train.10$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.10 quantile normalization  ############################################################################# 
RMA_train.10 <- normalizeBetweenArrays(RAW_train.10,method="quantile")      # quantile normalization
RMA_train.10 <- avereps(RMA_train.10, ID= RMA_train.10$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.10$E, col=pData.train.10$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.10 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.10$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.10$IVIM_Color, pch=16, cex=1.4) 

#### 2.10 COMBAT batch correction ############################################################################# 
batch.train.10 <- ifelse(pData.train.10$Batch>10,pData.train.10$Batch-1,pData.train.10$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.10)         
matrix.train.10.batch <- ComBat(dat=RMA_train.10$E, batch=batch.train.10, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.10.batch <- matrix.train.10.batch[row.names(Annotation.known),]
#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.10.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.10$IVIM_Color, pch=16, cex=1.3) 

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

rfe.10   # 22 optVars selected
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
GA.10   # 7 features selected at iteration 27
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 010 ###############################################################
################################################################################################
SAGA.CORE.10 <- list(SAGA.RFE=rfe.10$optVariables, SAGA.GA=GA.10$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.10 <- normalizeBetweenArrays(RAW_test.10, method="quantile")       # quantil normalize TestSet .10
RMA.10 <- avereps(RMA.10,ID= RMA.10$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.10 <- RMA.10$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.10  <- data.frame(labelDescription= rep(NA,dim(pData.test.10)[2]),row.names=colnames(pData.test.10))   # varMetadata: empty, but required 
phenoData.10 <- new("AnnotatedDataFrame",data=pData.test.10, varMetadata=metadata.10)   # annotatedDataFrame for the annotation of the samples
eset.gsea.10 <- ExpressionSet(assayData = matrix.gsea.10, phenoData = phenoData.10)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.10 <- ExpressionPhenoTest(eset.gsea.10,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.10 <- gsea(x=epheno.gsea.10, gsets=SAGA.CORE.10 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.10 <- summary(SAGA.GSEA.10)[,c(1,2,3,5,8)]
GSEA.RFE.10    <- subset(Output.GSEA.10,Output.GSEA.10$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.10) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.10     <- subset(Output.GSEA.10,Output.GSEA.10$geneSet == "SAGA.GA")
colnames(GSEA.GA.10) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.10 <- cbind(GSEA.RFE.10,GSEA.GA.10[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.10)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.10$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.10$GSEA_Vector <- Vector

pData.Test.sub.10           <- pData.test.10[pData.test.10$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.10$SampleID  <- row.names(pData.Test.sub.10)                   
GSEA.result.10              <- merge(pData.Test.sub.10,Output.GSEA.10, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_10.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.10),rows = NULL)
plot(SAGA.GSEA.10,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.11$E), col=pData.train.11$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.11 quantile normalization  ############################################################################# 
RMA_train.11 <- normalizeBetweenArrays(RAW_train.11,method="quantile")      # quantile normalization
RMA_train.11 <- avereps(RMA_train.11, ID= RMA_train.11$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.11$E, col=pData.train.11$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.11 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.11$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.11$IVIM_Color, pch=16, cex=1.4) 

#### 2.11 COMBAT batch correction ############################################################################# 
batch.train.11 <- ifelse(pData.train.11$Batch>11,pData.train.11$Batch-1,pData.train.11$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.11)         
matrix.train.11.batch <- ComBat(dat=RMA_train.11$E, batch=batch.train.11, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.11.batch <- matrix.train.11.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.11.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.11$IVIM_Color, pch=16, cex=1.3) 

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

rfe.11   # 21 variables found by SVM-rfe
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
GA.11   # 6 features selected at iteration 6
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
##### 7 Re-Train SVM on TrainingSet reduced to optimal predictors and predict TestSet 011 #####################
##############################################################################################################

#### 7.1 Quantile normalization of training set with bapred ################################################## 
##############################################################################################################
qunorm.train.11    <- qunormtrain(t(RAW_train.11$E))                                # quantile-normalize training samples only
matrix.train.11.qn <- log2(t(qunorm.train.11$xnorm))                                # extract quantile-normalized training samples
matrix.train.11.qn <- avereps(matrix.train.11.qn, ID= RAW_train.11$genes$ProbeName)  # average quadruplicates of training samples  
colnames(matrix.train.11.qn) <- row.names(pData.train.11)

boxplot(matrix.train.11.qn,col=pData.train.11$IVIM_Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

set.seed(12)
plot(Rtsne(t(matrix.train.11.qn),dims = 2, perplexity = 16,theta = 0.5, check_duplicates = FALSE, 
           pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,col=pData.train.11$IVIM_Color, pch=16, cex=1.3) 

#### 7.2 COMBAT batch correction of training set with bapred ################################################# 
##############################################################################################################
batch.train.11         <- as.factor(ifelse(pData.train.11$Batch>11,pData.train.11$Batch-1,pData.train.11$Batch))     # create batch factor 1...18 w/o Test samples
combat.train.11        <- combatba(t(matrix.train.11.qn), batch = batch.train.11)                                  # batch training samples only
matrix.train.11.bapred <- t(combat.train.11$xadj)                                                                 # extract batch corrected training matrix 

set.seed(12)
plot(Rtsne(t(matrix.train.11.bapred),dims = 2, perplexity = 16,theta = 0.5, check_duplicates = FALSE, 
           pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,col=pData.train.11$Design_Color, pch=16, cex=1.3) 

#### 7.3 subset training set for the optimal variables #######################################################
##############################################################################################################
matrix.train.11.bapred.GA   <- t(matrix.train.11.bapred[GA.11$optVariables,])
matrix.train.11.bapred.rfe  <- t(matrix.train.11.bapred[rfe.11$optVariables,])
matrix.train.11.bapred.full <- t(matrix.train.11.bapred[colnames(matrix.train.11),])

#### 7.4 train SVM on optimal predictors found by GA, rfe and all predictors  ################################
##############################################################################################################
set.seed(721)
svmOpt.GA.11  <- train(matrix.train.11.bapred.GA,labels.train.11,
                       method = "svmRadial",
                       metric = "Accuracy",
                       tuneLength = 20,
                       trControl = fullCtrl.11)
svmOpt.GA.11  


set.seed(721)
svmOpt.rfe.11  <- train(matrix.train.11.bapred.rfe,labels.train.11,
                        method = "svmRadial",
                        metric = "Accuracy",
                        tuneLength = 20,
                        trControl = fullCtrl.11)
svmOpt.rfe.11  


set.seed(721)
svmOpt.full.11  <- train(matrix.train.11.bapred.full,labels.train.11,
                         method = "svmRadial",
                         metric = "Accuracy",
                         tuneLength = 20,
                         trControl = fullCtrl.11)
svmOpt.full.11  

##############################################################################################################
#end of operations on the training set: preprocessing, feature selection, classifier development #############
##############################################################################################################


#### 7.4 Load and Add-on quantile normalization of the test set ##############################################
##############################################################################################################
matrix.test.11.qn  <- log2(t(qunormaddon(qunorm.train.11, t(RAW_test.11$E))))        # use qunormtrain object to normalize test samples / training set is fixed
matrix.test.11.qn  <- avereps(matrix.test.11.qn, ID= RAW_test.11$genes$ProbeName)    # average quadruplicates of test samples
boxplot(cbind(matrix.train.11.qn,matrix.test.11.qn),col=c(pData.train.11$IVIM_Color,pData.test.11$IVIM_Color),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)

set.seed(12)
plot(Rtsne(t(cbind(matrix.train.11.qn,matrix.test.11.qn)),dims = 2, perplexity = 16,theta = 0.5, check_duplicates = FALSE, 
           pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.11$IVIM_Color,pData.test.11$IVIM_Color), pch=16, cex=1.3) 

#### 7.5 Add-on batch correction of the test set #############################################################
##############################################################################################################
matrix.test.11.bapred  <- t(combatbaaddon(combat.train.11, t(matrix.test.11.qn), batch = as.factor(rep(1,length(row.names(pData.test.11))))))

set.seed(2)
plot(Rtsne(t(cbind(matrix.train.11.bapred,matrix.test.11.bapred)),dims = 2, perplexity = 16,
           theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 1000,verbose = FALSE, is_distance = FALSE)$Y,
     col=c(pData.train.11$Design_Color,pData.test.11$IVIM_Color), pch=16, cex=1.3)

#### 7.6 subset test set for the optimal predictors ##########################################################
##############################################################################################################
matrix.test.11.bapred.GA   <- t(matrix.test.11.bapred[GA.11$optVariables,])
matrix.test.11.bapred.rfe  <- t(matrix.test.11.bapred[rfe.11$optVariables,])
matrix.test.11.bapred.full <- t(matrix.test.11.bapred[colnames(matrix.train.11),])

#### 7.7 project test samples into optimal predictor PCA plot spanned by training set ########################
##############################################################################################################
pca.train <- prcomp(matrix.train.11.bapred.GA, center = T, scale. = T)           
plot(pca.train$x, pch=16, col=c(pData.train.11$Design_Color), cex=1.8, asp=1)
legend(-4.2,4.2, legend=unique(pData.train.11$Design), col=unique(pData.train.11$Design_Color), pch=16, bty="n", cex=1)
summary(pca.train)

coord.pca.test <- predict(pca.train, newdata = matrix.test.11.bapred.GA)
plot(rbind(pca.train$x,coord.pca.test), pch=16, col=c(pData.train.11$Design_Color,pData.test.11$IVIM_Color), cex=1.8, asp=1)
text(coord.pca.test, labels=row.names(pData.test.11), cex= 0.4, pos=3, offset = 0.3) 

#### 7.8 predict add-on adjusted test samples  ###############################################################
##############################################################################################################

Prediction_GA.11 <- predict(svmOpt.GA.11,matrix.test.11.bapred.GA, type = "prob")
Prediction_GA.11$Prediction_GA.11 <- ifelse(Prediction_GA.11$transforming>0.50,"transforming","untransforming")
Prediction_GA.11 <- cbind(pData.test.11[,c(1:3)],TrueLabel=pData.test.11$Class,Prediction_GA.11)
write.table(Prediction_GA.11, file = paste("TestSet011_Predictions_optVars_GA.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_SVM_full.11 <- predict(svmOpt.full.11, matrix.test.11.bapred.full, type = "prob")
Prediction_SVM_full.11$Prediction_SVM_full <- ifelse(Prediction_SVM_full.11$transforming>0.50,"transforming","untransforming")
Prediction_SVM_full.11 <- cbind(pData.test.11[,c(1:3)],TrueLabel=pData.test.11$Class,Prediction_SVM_full.11)
write.table(Prediction_SVM_full.11, file = paste("TestSet011_Predictions_allVars.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

Prediction_rfe.11 <- predict(svmOpt.rfe.11,matrix.test.11.bapred.rfe, type = "prob")
Prediction_rfe.11$Prediction_rfe.11 <- ifelse(Prediction_rfe.11$transforming>0.50,"transforming","untransforming")
Prediction_rfe.11 <- cbind(pData.test.11[,c(1:3)],TrueLabel=pData.test.11$Class,Prediction_rfe.11)
write.table(Prediction_rfe.11, file = paste("TestSet011_Predictions_optVars_rfe.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

#### 7.10 Performance of optVars Classifier on  TestSet 011 #####################################################
#############################################################################################################

sink("TestSet011_ConfusionMatrix_rfe.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_rfe.11$Prediction_rfe.11), as.factor(Prediction_rfe.11$TrueLabel))
sink()

sink("TestSet011_ConfusionMatrix_GA.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_GA.11$Prediction_GA.11), as.factor(Prediction_GA.11$TrueLabel))
sink()

sink("TestSet011_ConfusionMatrix_FULL.txt", append = TRUE)
confusionMatrix(as.factor(Prediction_SVM_full.11$Prediction_SVM_full), as.factor(Prediction_SVM_full.11$TrueLabel))
sink()

#### 7.11 ROC on probability "transforming" TestSet 011 ########################################################

Prediction_rfe.11$Class <- as.factor(ifelse(Prediction_rfe.11$TrueLabel == "transforming","transforming","nontransforming"))
roc.rfe.11 <- roc(Prediction_rfe.11$Class,                    
                  Prediction_rfe.11$transforming,             
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
                  print.auc=T)

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


#### 7.12 Precision-Recall curve on probability "transforming" TestSet 011 ####################################
Prediction_rfe.11$Class_Code <- ifelse(Prediction_rfe.11$TrueLabel == "transforming",1,0)
pr.rfe.11 <- pr.curve(scores.class0 = Prediction_rfe.11$transforming , weights.class0 = Prediction_rfe.11$Class_Code, curve = TRUE,rand.compute = T)
plot(pr.rfe.11,rand.plot = TRUE, legend = F, color = 1, main = "")

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
boxplot(log2(RAW_train.12$E), col=pData.train.12$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.12 quantile normalization  ############################################################################# 
RMA_train.12 <- normalizeBetweenArrays(RAW_train.12,method="quantile")      # quantile normalization
RMA_train.12 <- avereps(RMA_train.12, ID= RMA_train.12$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.12$E, col=pData.train.12$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.12 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.12$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.12$IVIM_Color, pch=16, cex=1.4) 

#### 2.12 COMBAT batch correction ############################################################################# 
batch.train.12 <- ifelse(pData.train.12$Batch>12,pData.train.12$Batch-1,pData.train.12$Batch)                       
modcombat      <- model.matrix(~1, data=pData.train.12)         
matrix.train.12.batch <- ComBat(dat=RMA_train.12$E, batch=batch.train.12, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.12.batch <- matrix.train.12.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.12.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.12$IVIM_Color, pch=16, cex=1.3) 

##############################################################################################################
#### 3. nonspecific feature prefiltering TrainingSet 012 ######################################################
##############################################################################################################
fselect.12  <- genefilter(matrix.train.12.batch, filterfun(f1))
summary(fselect.12)   # 1249 predictors
matrix.train.12 <-matrix.train.12.batch[fselect.12,]

##############################################################################################################
#### 2. SVM: FULL MODEL TrainingSet 012 #######################################################################
##############################################################################################################
matrix.train.12 <- (t(matrix.train.12))
labels.train.12 <- as.factor(pData.train.12$Class)

## create 200 resamples of the train data (10TrainingSetCVn20) - the same index is used for SVM-rfe for comparison
set.seed(1234)
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

rfe.12   # 14 variables found by SVM-rfe
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
GA.12   # 6 features selected at iteration 11
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



################################################################################################
#### 8. SAGA-GSEA for TestSet 012 ###############################################################
################################################################################################
SAGA.CORE.12 <- list(SAGA.RFE=rfe.12$optVariables, SAGA.GA=GA.12$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.12 <- normalizeBetweenArrays(RAW_test.12, method="quantile")       # quantil normalize TestSet .12
RMA.12 <- avereps(RMA.12,ID= RMA.12$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.12 <- RMA.12$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.12  <- data.frame(labelDescription= rep(NA,dim(pData.test.12)[2]),row.names=colnames(pData.test.12))   # varMetadata: empty, but required 
phenoData.12 <- new("AnnotatedDataFrame",data=pData.test.12, varMetadata=metadata.12)   # annotatedDataFrame for the annotation of the samples
eset.gsea.12 <- ExpressionSet(assayData = matrix.gsea.12, phenoData = phenoData.12)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.12 <- ExpressionPhenoTest(eset.gsea.12,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.12 <- gsea(x=epheno.gsea.12, gsets=SAGA.CORE.12 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.12 <- summary(SAGA.GSEA.12)[,c(1,2,3,5,8)]
GSEA.RFE.12    <- subset(Output.GSEA.12,Output.GSEA.12$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.12) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.12     <- subset(Output.GSEA.12,Output.GSEA.12$geneSet == "SAGA.GA")
colnames(GSEA.GA.12) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.12 <- cbind(GSEA.RFE.12,GSEA.GA.12[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.12)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.12$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.12$GSEA_Vector <- Vector

pData.Test.sub.12           <- pData.test.12[pData.test.12$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.12$SampleID  <- row.names(pData.Test.sub.12)                   
GSEA.result.12              <- merge(pData.Test.sub.12,Output.GSEA.12, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_12.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.12),rows = NULL)
plot(SAGA.GSEA.12,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.13$E), col=pData.train.13$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.13 quantile normalization  ############################################################################# 
RMA_train.13 <- normalizeBetweenArrays(RAW_train.13,method="quantile")      # quantile normalization
RMA_train.13 <- avereps(RMA_train.13, ID= RMA_train.13$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.13$E, col=pData.train.13$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.13 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.13$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.13$IVIM_Color, pch=16, cex=1.4) 

#### 2.13 COMBAT batch correction ############################################################################# 
batch.train.13 <- ifelse(pData.train.13$Batch>13,pData.train.13$Batch-1,pData.train.13$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.13)         
matrix.train.13.batch <- ComBat(dat=RMA_train.13$E, batch=batch.train.13, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.13.batch <- matrix.train.13.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.13.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.13$IVIM_Color, pch=16, cex=1.3) 

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

rfe.13   # 27 variables found by SVM-rfe
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
GA.13   # yields 11 features
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

### plot external accuracy over the iterations ################################################################
performance.aggr.13 <- subset(performance.long.13,performance.long.13$Group =="external")

ggplot(performance.aggr.13, aes(Iter, Accuracy)) +
  geom_point() +
  geom_smooth(span = 0.7,se = T) +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        axis.text = element_text(size=14, color ="black"),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())


################################################################################################
#### 8. SAGA-GSEA for TestSet 013 ###############################################################
################################################################################################
SAGA.CORE.13 <- list(SAGA.RFE=rfe.13$optVariables, SAGA.GA=GA.13$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.13 <- normalizeBetweenArrays(RAW_test.13, method="quantile")       # quantil normalize TestSet .13
RMA.13 <- avereps(RMA.13,ID= RMA.13$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.13 <- RMA.13$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.13  <- data.frame(labelDescription= rep(NA,dim(pData.test.13)[2]),row.names=colnames(pData.test.13))   # varMetadata: empty, but required 
phenoData.13 <- new("AnnotatedDataFrame",data=pData.test.13, varMetadata=metadata.13)   # annotatedDataFrame for the annotation of the samples
eset.gsea.13 <- ExpressionSet(assayData = matrix.gsea.13, phenoData = phenoData.13)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.13 <- ExpressionPhenoTest(eset.gsea.13,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.13 <- gsea(x=epheno.gsea.13, gsets=SAGA.CORE.13 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.13 <- summary(SAGA.GSEA.13)[,c(1,2,3,5,8)]
GSEA.RFE.13    <- subset(Output.GSEA.13,Output.GSEA.13$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.13) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.13     <- subset(Output.GSEA.13,Output.GSEA.13$geneSet == "SAGA.GA")
colnames(GSEA.GA.13) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.13 <- cbind(GSEA.RFE.13,GSEA.GA.13[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.13)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.13$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.13$GSEA_Vector <- Vector

pData.Test.sub.13           <- pData.test.13[pData.test.13$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.13$SampleID  <- row.names(pData.Test.sub.13)                   
GSEA.result.13              <- merge(pData.Test.sub.13,Output.GSEA.13, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_13.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.13),rows = NULL)
plot(SAGA.GSEA.13,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.14$E), col=pData.train.14$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.14 quantile normalization  ############################################################################# 
RMA_train.14 <- normalizeBetweenArrays(RAW_train.14,method="quantile")      # quantile normalization
RMA_train.14 <- avereps(RMA_train.14, ID= RMA_train.14$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.14$E, col=pData.train.14$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.14 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.14$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.14$IVIM_Color, pch=16, cex=1.4) 

#### 2.14 COMBAT batch correction ############################################################################# 
batch.train.14 <- ifelse(pData.train.14$Batch>14,pData.train.14$Batch-1,pData.train.14$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.14)         
matrix.train.14.batch <- ComBat(dat=RMA_train.14$E, batch=batch.train.14, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.14.batch <- matrix.train.14.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.14.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.14$IVIM_Color, pch=16, cex=1.3) 

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

rfe.14   # 21 variables found by SVM-rfe
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
GA.14   # 10 features selected at iteration 33
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 014 ###############################################################
################################################################################################
SAGA.CORE.14 <- list(SAGA.RFE=rfe.14$optVariables, SAGA.GA=GA.14$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.14 <- normalizeBetweenArrays(RAW_test.14, method="quantile")       # quantil normalize TestSet .14
RMA.14 <- avereps(RMA.14,ID= RMA.14$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.14 <- RMA.14$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.14  <- data.frame(labelDescription= rep(NA,dim(pData.test.14)[2]),row.names=colnames(pData.test.14))   # varMetadata: empty, but required 
phenoData.14 <- new("AnnotatedDataFrame",data=pData.test.14, varMetadata=metadata.14)   # annotatedDataFrame for the annotation of the samples
eset.gsea.14 <- ExpressionSet(assayData = matrix.gsea.14, phenoData = phenoData.14)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.14 <- ExpressionPhenoTest(eset.gsea.14,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.14 <- gsea(x=epheno.gsea.14, gsets=SAGA.CORE.14 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.14 <- summary(SAGA.GSEA.14)[,c(1,2,3,5,8)]
GSEA.RFE.14    <- subset(Output.GSEA.14,Output.GSEA.14$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.14) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.14     <- subset(Output.GSEA.14,Output.GSEA.14$geneSet == "SAGA.GA")
colnames(GSEA.GA.14) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.14 <- cbind(GSEA.RFE.14,GSEA.GA.14[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.14)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.14$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.14$GSEA_Vector <- Vector

pData.Test.sub.14           <- pData.test.14[pData.test.14$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.14$SampleID  <- row.names(pData.Test.sub.14)                   
GSEA.result.14              <- merge(pData.Test.sub.14,Output.GSEA.14, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_14.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.14),rows = NULL)
plot(SAGA.GSEA.14,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.15$E), col=pData.train.15$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.15 quantile normalization  ############################################################################# 
RMA_train.15 <- normalizeBetweenArrays(RAW_train.15,method="quantile")      # quantile normalization
RMA_train.15 <- avereps(RMA_train.15, ID= RMA_train.15$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.15$E, col=pData.train.15$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.15 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.15$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.15$IVIM_Color, pch=16, cex=1.4) 

#### 2.15 COMBAT batch correction ############################################################################# 
batch.train.15 <- ifelse(pData.train.15$Batch>15,pData.train.15$Batch-1,pData.train.15$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.15)         
matrix.train.15.batch <- ComBat(dat=RMA_train.15$E, batch=batch.train.15, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.15.batch <- matrix.train.15.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.15.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.15$IVIM_Color, pch=16, cex=1.3) 

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

rfe.15   # 17 variables found by SVM-rfe
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
GA.15   # 8 features selected at iteration 27
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 015 ###############################################################
################################################################################################
SAGA.CORE.15 <- list(SAGA.RFE=rfe.15$optVariables, SAGA.GA=GA.15$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.15 <- normalizeBetweenArrays(RAW_test.15, method="quantile")       # quantil normalize TestSet .15
RMA.15 <- avereps(RMA.15,ID= RMA.15$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.15 <- RMA.15$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.15  <- data.frame(labelDescription= rep(NA,dim(pData.test.15)[2]),row.names=colnames(pData.test.15))   # varMetadata: empty, but required 
phenoData.15 <- new("AnnotatedDataFrame",data=pData.test.15, varMetadata=metadata.15)   # annotatedDataFrame for the annotation of the samples
eset.gsea.15 <- ExpressionSet(assayData = matrix.gsea.15, phenoData = phenoData.15)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.15 <- ExpressionPhenoTest(eset.gsea.15,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.15 <- gsea(x=epheno.gsea.15, gsets=SAGA.CORE.15 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.15 <- summary(SAGA.GSEA.15)[,c(1,2,3,5,8)]
GSEA.RFE.15    <- subset(Output.GSEA.15,Output.GSEA.15$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.15) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.15     <- subset(Output.GSEA.15,Output.GSEA.15$geneSet == "SAGA.GA")
colnames(GSEA.GA.15) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.15 <- cbind(GSEA.RFE.15,GSEA.GA.15[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.15)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.15$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.15$GSEA_Vector <- Vector

pData.Test.sub.15           <- pData.test.15[pData.test.15$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.15$SampleID  <- row.names(pData.Test.sub.15)                   
GSEA.result.15              <- merge(pData.Test.sub.15,Output.GSEA.15, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_15.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.15),rows = NULL)
plot(SAGA.GSEA.15,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.16$E), col=pData.train.16$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.16 quantile normalization  ############################################################################# 
RMA_train.16 <- normalizeBetweenArrays(RAW_train.16,method="quantile")      # quantile normalization
RMA_train.16 <- avereps(RMA_train.16, ID= RMA_train.16$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.16$E, col=pData.train.16$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.16 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.16$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.16$IVIM_Color, pch=16, cex=1.4) 

#### 2.16 COMBAT batch correction ############################################################################# 
batch.train.16 <- ifelse(pData.train.16$Batch>16,pData.train.16$Batch-1,pData.train.16$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.16)         
matrix.train.16.batch <- ComBat(dat=RMA_train.16$E, batch=batch.train.16, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.16.batch <- matrix.train.16.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.16.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.16$IVIM_Color, pch=16, cex=1.3) 

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

rfe.16   # 17 variables found by SVM-rfe
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
GA.16   # 11 features selected at iteration 3
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


################################################################################################
#### 8. SAGA-GSEA for TestSet 016 ###############################################################
################################################################################################
SAGA.CORE.16 <- list(SAGA.RFE=rfe.16$optVariables, SAGA.GA=GA.16$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.16 <- normalizeBetweenArrays(RAW_test.16, method="quantile")       # quantil normalize TestSet .16
RMA.16 <- avereps(RMA.16,ID= RMA.16$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.16 <- RMA.16$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.16  <- data.frame(labelDescription= rep(NA,dim(pData.test.16)[2]),row.names=colnames(pData.test.16))   # varMetadata: empty, but required 
phenoData.16 <- new("AnnotatedDataFrame",data=pData.test.16, varMetadata=metadata.16)   # annotatedDataFrame for the annotation of the samples
eset.gsea.16 <- ExpressionSet(assayData = matrix.gsea.16, phenoData = phenoData.16)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.16 <- ExpressionPhenoTest(eset.gsea.16,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.16 <- gsea(x=epheno.gsea.16, gsets=SAGA.CORE.16 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.16 <- summary(SAGA.GSEA.16)[,c(1,2,3,5,8)]
GSEA.RFE.16    <- subset(Output.GSEA.16,Output.GSEA.16$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.16) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.16     <- subset(Output.GSEA.16,Output.GSEA.16$geneSet == "SAGA.GA")
colnames(GSEA.GA.16) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.16 <- cbind(GSEA.RFE.16,GSEA.GA.16[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.16)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.16$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.16$GSEA_Vector <- Vector

pData.Test.sub.16           <- pData.test.16[pData.test.16$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.16$SampleID  <- row.names(pData.Test.sub.16)                   
GSEA.result.16              <- merge(pData.Test.sub.16,Output.GSEA.16, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_16.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.16),rows = NULL)
plot(SAGA.GSEA.16,es.nes='nes')
dev.off()



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
boxplot(log2(RAW_train.17$E), col=pData.train.17$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.17 quantile normalization  ############################################################################# 
RMA_train.17 <- normalizeBetweenArrays(RAW_train.17,method="quantile")      # quantile normalization
RMA_train.17 <- avereps(RMA_train.17, ID= RMA_train.17$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.17$E, col=pData.train.17$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.17 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.17$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.17$IVIM_Color, pch=16, cex=1.4) 

#### 2.17 COMBAT batch correction ############################################################################# 
batch.train.17 <- ifelse(pData.train.17$Batch>17,pData.train.17$Batch-1,pData.train.17$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.17)         
matrix.train.17.batch <- ComBat(dat=RMA_train.17$E, batch=batch.train.17, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.17.batch <- matrix.train.17.batch[row.names(Annotation.known),]
matrix.train.17.batch <- round(matrix.train.17.batch,3)

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.17.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.17$IVIM_Color, pch=16, cex=1.3) 

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

rfe.17   # all 1246 variables selected by SVM-rfe  
write.table(rfe.17$results, file = "TrainingSet017_Results_rfe.txt", sep="\t",col.names=NA)

trellis.par.set(caretTheme())
plot(rfe.17, type = c("g", "o"))
plot(rfe.17, type = c("g", "o"), xlim = c(0,61))

# manually select the second best solution (number of predictors=21)
# 1246 predictors: CV-accuracy = 0.9031
# 21 predictors :  CV-accuracy = 0.8980

varImp_17 <- varImp(svmFull.17)$importance                              # extract the variable importance (AUCROC) of all 1246 predictors 
varImp_17 <- varImp_17[order(varImp_17$transforming,decreasing = T),]   # rank by AUCROC
varImp_17 <- cbind(varImp_17,Annotation[row.names(varImp_17),])        
optFeatures.17 <- varImp_17[c(1:21),]                                   # select the 21 most important predictors
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
GA.17   # 10 features selected at iteration 3
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


################################################################################################
#### 8. SAGA-GSEA for TestSet 017 ###############################################################
################################################################################################
SAGA.CORE.17 <- list(SAGA.RFE=row.names(optFeatures.17), SAGA.GA=GA.17$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.17 <- normalizeBetweenArrays(RAW_test.17, method="quantile")       # quantil normalize TestSet .17
RMA.17 <- avereps(RMA.17,ID= RMA.17$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.17 <- RMA.17$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.17  <- data.frame(labelDescription= rep(NA,dim(pData.test.17)[2]),row.names=colnames(pData.test.17))   # varMetadata: empty, but required 
phenoData.17 <- new("AnnotatedDataFrame",data=pData.test.17, varMetadata=metadata.17)   # annotatedDataFrame for the annotation of the samples
eset.gsea.17 <- ExpressionSet(assayData = matrix.gsea.17, phenoData = phenoData.17)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.17 <- ExpressionPhenoTest(eset.gsea.17,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.17 <- gsea(x=epheno.gsea.17, gsets=SAGA.CORE.17 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.17 <- summary(SAGA.GSEA.17)[,c(1,2,3,5,8)]
GSEA.RFE.17    <- subset(Output.GSEA.17,Output.GSEA.17$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.17) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.17     <- subset(Output.GSEA.17,Output.GSEA.17$geneSet == "SAGA.GA")
colnames(GSEA.GA.17) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.17 <- cbind(GSEA.RFE.17,GSEA.GA.17[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.17)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.17$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.17$GSEA_Vector <- Vector

pData.Test.sub.17           <- pData.test.17[pData.test.17$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.17$SampleID  <- row.names(pData.Test.sub.17)                   
GSEA.result.17              <- merge(pData.Test.sub.17,Output.GSEA.17, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_17.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.17),rows = NULL)
plot(SAGA.GSEA.17,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.18$E), col=pData.train.18$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.18 quantile normalization  ############################################################################# 
RMA_train.18 <- normalizeBetweenArrays(RAW_train.18,method="quantile")      # quantile normalization
RMA_train.18 <- avereps(RMA_train.18, ID= RMA_train.18$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.18$E, col=pData.train.18$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.18 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.18$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.18$IVIM_Color, pch=16, cex=1.4) 

#### 2.18 COMBAT batch correction ############################################################################# 
batch.train.18 <- ifelse(pData.train.18$Batch>18,pData.train.18$Batch-1,pData.train.18$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.18)         
matrix.train.18.batch <- ComBat(dat=RMA_train.18$E, batch=batch.train.18, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.18.batch <- matrix.train.18.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.18.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.18$IVIM_Color, pch=16, cex=1.3) 

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
#### 4. GENETIC ALGORITHM FOR REFINED FEATURE SELECTION TrainingSet 018 ########################
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
GA.18   # 8 features selected at iteration 25
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


################################################################################################
#### 8. SAGA-GSEA for TestSet 018 ###############################################################
################################################################################################
SAGA.CORE.18 <- list(SAGA.RFE=rfe.18$optVariables, SAGA.GA=GA.18$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.18 <- normalizeBetweenArrays(RAW_test.18, method="quantile")       # quantil normalize TestSet .18
RMA.18 <- avereps(RMA.18,ID= RMA.18$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.18 <- RMA.18$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.18  <- data.frame(labelDescription= rep(NA,dim(pData.test.18)[2]),row.names=colnames(pData.test.18))   # varMetadata: empty, but required 
phenoData.18 <- new("AnnotatedDataFrame",data=pData.test.18, varMetadata=metadata.18)   # annotatedDataFrame for the annotation of the samples
eset.gsea.18 <- ExpressionSet(assayData = matrix.gsea.18, phenoData = phenoData.18)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.18 <- ExpressionPhenoTest(eset.gsea.18,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.18 <- gsea(x=epheno.gsea.18, gsets=SAGA.CORE.18 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.18 <- summary(SAGA.GSEA.18)[,c(1,2,3,5,8)]
GSEA.RFE.18    <- subset(Output.GSEA.18,Output.GSEA.18$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.18) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.18     <- subset(Output.GSEA.18,Output.GSEA.18$geneSet == "SAGA.GA")
colnames(GSEA.GA.18) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.18 <- cbind(GSEA.RFE.18,GSEA.GA.18[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.18)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.18$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.18$GSEA_Vector <- Vector

pData.Test.sub.18           <- pData.test.18[pData.test.18$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.18$SampleID  <- row.names(pData.Test.sub.18)                   
GSEA.result.18              <- merge(pData.Test.sub.18,Output.GSEA.18, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_18.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.18),rows = NULL)
plot(SAGA.GSEA.18,es.nes='nes')
dev.off()


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
boxplot(log2(RAW_train.19$E), col=pData.train.19$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.19 quantile normalization  ############################################################################# 
RMA_train.19 <- normalizeBetweenArrays(RAW_train.19,method="quantile")      # quantile normalization
RMA_train.19 <- avereps(RMA_train.19, ID= RMA_train.19$genes$ProbeName)      # average over ProbeIDs    
boxplot(RMA_train.19$E, col=pData.train.19$IVIM_Color,boxwex=0.6,cex.axis=0.4,las=2,outline=FALSE)    

#### 2.19 visualize quantile normalized data  ################################################################# 
set.seed(12)
plot(Rtsne(t(RMA_train.19$E),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.19$IVIM_Color, pch=16, cex=1.4) 

#### 2.19 COMBAT batch correction ############################################################################# 
batch.train.19 <- ifelse(pData.train.19$Batch>19,pData.train.19$Batch-1,pData.train.19$Batch)                       
modcombat     <- model.matrix(~1, data=pData.train.19)         
matrix.train.19.batch <- ComBat(dat=RMA_train.19$E, batch=batch.train.19, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
matrix.train.19.batch <- matrix.train.19.batch[row.names(Annotation.known),]

#### 2.4 t-SNE of batch corrected dataset ####################################################################
set.seed(12)
plot(Rtsne(t(matrix.train.19.batch),dims = 2, perplexity = 16,check_duplicates = FALSE,pca = TRUE, max_iter = 1000,is_distance = FALSE)$Y,
     col=pData.train.19$IVIM_Color, pch=16, cex=1.3) 

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

rfe.19   # 50 variables found by SVM-rfe
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
GA.19   # 18 features selected at iteration 6 
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

################################################################################################
#### 8. SAGA-GSEA for TestSet 019 ###############################################################
################################################################################################
SAGA.CORE.19 <- list(SAGA.RFE=rfe.19$optVariables, SAGA.GA=GA.19$optVariables)
#### 8.1. Normalize, average ###################################################################
RMA.19 <- normalizeBetweenArrays(RAW_test.19, method="quantile")       # quantil normalize TestSet .19
RMA.19 <- avereps(RMA.19,ID= RMA.19$genes$ProbeName)                    # average replicates to one value for each probe
matrix.gsea.19 <- RMA.19$E                                             # extract log2 expression values 

#### 8.2. make ExpressionSet (Biobase) object ##################################################
metadata.19  <- data.frame(labelDescription= rep(NA,dim(pData.test.19)[2]),row.names=colnames(pData.test.19))   # varMetadata: empty, but required 
phenoData.19 <- new("AnnotatedDataFrame",data=pData.test.19, varMetadata=metadata.19)   # annotatedDataFrame for the annotation of the samples
eset.gsea.19 <- ExpressionSet(assayData = matrix.gsea.19, phenoData = phenoData.19)     # this is the ExpressionSet required for phenoTest

#### 8.3. make ePheno object: contains the FCs associated with vector variable ##################
vars2test     <- list(ordinal="GSEA_Vector")                         # Variables (here: GSEA_Vectors) to test against MOCK, which is always GSEA_Vectors = 1 in the SIF 
epheno.gsea.19 <- ExpressionPhenoTest(eset.gsea.19,vars2test,p.adjust.method='BH')

#### 8.4 GSEA #################################################################################
SAGA.GSEA.19 <- gsea(x=epheno.gsea.19, gsets=SAGA.CORE.19 ,B=2000,      # calculate GSEA-scores based on the FC in the epheno object
                     center = TRUE, test = "perm", p.adjust.method='BH', minGenes = 5)
#### 8.5 output ###############################################################################
Output.GSEA.19 <- summary(SAGA.GSEA.19)[,c(1,2,3,5,8)]
GSEA.RFE.19    <- subset(Output.GSEA.19,Output.GSEA.19$geneSet == "SAGA.RFE")                
colnames(GSEA.RFE.19) <- c("variable","geneSet.RFE","n.RFE","nes.RFE","fdr.RFE")
GSEA.GA.19     <- subset(Output.GSEA.19,Output.GSEA.19$geneSet == "SAGA.GA")
colnames(GSEA.GA.19) <- c("variable","geneSet.GA","n.GA","nes.GA","fdr.GA")
Output.GSEA.19 <- cbind(GSEA.RFE.19,GSEA.GA.19[,-1])

Vector <- NULL    ### pull out the Vector index number from the result table                  
for (a in 1:nrow(Output.GSEA.19)) {Vector[a] <- unlist(strsplit(as.character(Output.GSEA.19$variable[a]), ".", fixed = TRUE))[2] }
Output.GSEA.19$GSEA_Vector <- Vector

pData.Test.sub.19           <- pData.test.19[pData.test.19$GSEA_Vector != 1, ]    # pData.Test minus the Mock samples                   
pData.Test.sub.19$SampleID  <- row.names(pData.Test.sub.19)                   
GSEA.result.19              <- merge(pData.Test.sub.19,Output.GSEA.19, by.x="GSEA_Vector", by.y = "GSEA_Vector") 

# make pdf report
pdf(file="SAGA.GSEA_Batch_19.pdf",useDingbats = F,width = 10, height = 10)  
grid.table(summary(SAGA.GSEA.19),rows = NULL)
plot(SAGA.GSEA.19,es.nes='nes')
dev.off()


#############################################################################################################################################
#############################################################################################################################################
#### Analyze GSEA leave-one-batch-out results: Figure 6 #####################################################################################
#############################################################################################################################################
#############################################################################################################################################

#### 1 aggregate all results over 18 iterations: Supplementary Data 7 tab 2 #################################################################
#############################################################################################################################################
# Note: batch 11 (IVIM #171102) excluded (no mock control available)
Results_GSEA_all <- rbind(GSEA.result.1,GSEA.result.2,GSEA.result.3,GSEA.result.4,GSEA.result.5,GSEA.result.6,GSEA.result.7,GSEA.result.8,
                          GSEA.result.9,GSEA.result.10,GSEA.result.12,GSEA.result.13,GSEA.result.14,GSEA.result.15,GSEA.result.16,
                          GSEA.result.17,GSEA.result.18,GSEA.result.19)

Results_GSEA_all$TrueLabel       <- Results_GSEA_all$Class   
Results_GSEA_all$Class           <- ifelse(Results_GSEA_all$TrueLabel == "transforming","transforming","nontransforming")

### For batch 2 no SVM-GA was performed, since RFE arrived at 9 predictors already ==> copy results from SVM-RFE batch 2 to results SVM-GA batch 2:
Results_GSEA_all[c(5:11),c(34:37)] <- Results_GSEA_all[c(5:11),c(30:33)]


#### 2 ROC / AURPRC curves ################################################################################################################## 
#############################################################################################################################################

#### 2.1 GSEA_RFE vs GSEA_GA complete data ###################################################### 
ROC.GSEA.RFE.ALL <- roc(Results_GSEA_all$Class,           
                    Results_GSEA_all$nes.RFE,     
                    percent=TRUE, levels=c("nontransforming","transforming"),
                    plot=T, auc.polygon=F, max.auc.polygon=F, col = "#4A6893",lwd=4, grid=F,
                    print.auc=T,print.thres="best" )   # best NES = 1.9 AUC = 89%

ROC.GSEA.GA.ALL <- roc(Results_GSEA_all$Class,           
                       Results_GSEA_all$nes.GA,     
                       percent=TRUE, levels=c("nontransforming","transforming"),
                       plot=T, auc.polygon=F, max.auc.polygon=F, col = "#CB4848",lwd=4, grid=F,
                       print.auc=T,print.thres="best", add=TRUE )   # best NES = 1.7 AUC = 91.4%

#### 2.2 GSEA_RFE vs GSEA_GA without LTR.RV.SFFV ################################################ 
SF91                <- "A1_LTR.RV.SF.eGFP"
Results_GSEA_woSF91 <- subset(Results_GSEA_all, !Results_GSEA_all$Vector %in% SF91 ) 

ROC.GSEA.RFE.woSF91 <- roc(Results_GSEA_woSF91$Class,           
                        Results_GSEA_woSF91$nes.RFE,     
                        percent=TRUE, levels=c("nontransforming","transforming"),
                        plot=T, auc.polygon=F, max.auc.polygon=F, col = "#4A6893",lwd=4, grid=F,
                        print.auc=T,print.thres="best" )     # best NES = 0.1 AUC = 82.4%

ROC.GSEA.GA.woSF91 <- roc(Results_GSEA_woSF91$Class,           
                       Results_GSEA_woSF91$nes.GA,     
                       percent=TRUE, levels=c("nontransforming","transforming"),
                       plot=T, auc.polygon=F, max.auc.polygon=F, col = "#CB4848",lwd=4, grid=F,
                       print.auc=T,print.thres="best", add=TRUE )   # best NES = 1.3 AUC = 85.4%


#### 2.3 Figure 6d GSEA_GA All vs LTR.RV.SFFV ################################################### 
ROC.GSEA.GA.ALL <- roc(Results_GSEA_all$Class,           
                       Results_GSEA_all$nes.GA,     
                       percent=TRUE, levels=c("nontransforming","transforming"),
                       plot=T, auc.polygon=F, max.auc.polygon=F, col = "#CB4848",lwd=4, grid=F,
                       print.auc=T,print.thres="best", add=F )   # best NES = 1.7 AUC = 91.4%

ROC.GSEA.GA.woSF91 <- roc(Results_GSEA_woSF91$Class,           
                          Results_GSEA_woSF91$nes.GA,     
                          percent=TRUE, levels=c("nontransforming","transforming"),
                          plot=T, auc.polygon=F, max.auc.polygon=F, col = "#686868",lwd=4, grid=F,
                          print.auc=T,print.thres="best", add=TRUE )   # best NES = 1.3 AUC = 85.4%


#### 2.4 Figure 6e GSEA_GA All vs IVIM ######################################################### 
ROC.GSEA.GA.ALL <- roc(Results_GSEA_all$Class,           
                       Results_GSEA_all$nes.GA,     
                       percent=TRUE, levels=c("nontransforming","transforming"),
                       plot=T, auc.polygon=F, max.auc.polygon=F, col = "#CB4848",lwd=4, grid=F,
                       print.auc=T,print.thres=1.3, add=F )    

IVIM.total <- roc(IVIM$Classification_Code,    
                  IVIM$MTT_Score,              
                  percent=TRUE, smooth = F,
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                  col = "#8285BC", grid=F, lwd = 4, cex.lab=1, 
                  print.auc=T, print.thres = 3, add = T )

roc.test(ROC.GSEA.GA.ALL, IVIM.total, alternative = "greater")   # p-value = 0.0056


### 2.5 Figure 6f GSEA_GA All vs IVIM PRROC ######################################################

Results_GSEA_all$Class_Code <- ifelse(Results_GSEA_all$TrueLabel == "transforming",1,0)
AUPRC.SAGA.GSEA <- pr.curve(scores.class0 = Results_GSEA_all$nes.GA , weights.class0 = Results_GSEA_all$Class_Code, curve = TRUE,rand.compute = T)
plot(AUPRC.SAGA.GSEA, rand.plot = TRUE, legend = F, color = "#CB4848", main = "",auc.main = FALSE)  # AUPRC = 0.944  / AUC random = 0.422


IVIM$Class_Code <- ifelse(IVIM$TrueLabel == "transforming",1,0)
AUPRC.IVIM <- pr.curve(scores.class0 = IVIM$MTT_Score, weights.class0 = as.numeric(IVIM$Class_Code), curve = TRUE,rand.compute = T)
plot(AUPRC.IVIM, rand.plot = TRUE, legend = F, color = "#8285BC", add = T) # AUPRC = 0.892  / AUC random = 0.582


#### 3 Confusion matrix ##################################################################################################################### 
#############################################################################################################################################

# best threshold for GSEA.GA = NES >1.3  
Results_GSEA_all$GSEA.GA_Prediction <- ifelse(Results_GSEA_all$nes.GA > 1.3,"transforming","untransforming")

# best threshold for GSEA.RFE = NES >0.1  
Results_GSEA_all$GSEA.RFE_Prediction <- ifelse(Results_GSEA_all$nes.RFE > 0.1,"transforming","untransforming")

#Supplementary Data 7 tab 1: 
sink("AllTestSets_ConfusionMatrix_GSEA_GA.txt", append = TRUE)
confusionMatrix(as.factor(Results_GSEA_all$GSEA.GA_Prediction), as.factor(Results_GSEA_all$TrueLabel))
sink()

sink("AllTestSets_ConfusionMatrix_GSEA_RFE.txt", append = TRUE)
confusionMatrix(as.factor(Results_GSEA_all$GSEA.RFE_Prediction), as.factor(Results_GSEA_all$TrueLabel))
sink()

write.table(Results_GSEA_all, file = "AllResults_SAGA.GSEA_leave one batch out_FINAL.txt", sep="\t",col.names=NA)


