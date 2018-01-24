###### Feature selection ######################################################################
###############################################################################################
library(e1071)
library(limma)
library(sva)
library(Rtsne)
library(RColorBrewer)
library(caret)
library(ggplot2)
library(dplyr)
library(doMC)
registerDoMC(cores = 7)


################################################################################################
#### 1. Data handling ##########################################################################
################################################################################################

pData      <- read.delim("Targets_full_170906.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation <- read.delim("SAGA_Annotation.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
Annotation <- subset(Annotation,Annotation$NAME !="Unknown")  # keep only annotated features: leaves 35492 genes
SAGA_Data  <- read.delim("SAGA_RAW_Data.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE)  # read in the RAW data 
SAGA_RAW   <- log2(as.matrix(SAGA_Data[,-1]))     
row.names(SAGA_RAW) <- SAGA_Data$PROBE_ID

### 1.1 Preprocessing ##########################################################################
################################################################################################

all(row.names(pData) == colnames(SAGA_RAW))   
boxplot(SAGA_RAW, col=pData$IVIM_Color, names=pData$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

### 1.1.1 Quantile Normalization  ##############################################################
eset.rma <- normalizeBetweenArrays(SAGA_RAW,method="quantile")      
eset.rma <- avereps(eset.rma, ID= row.names(eset.rma))             
eset.rma <- eset.rma[row.names(Annotation),]
boxplot(eset.rma,col=pData$IVIM_Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

### 1.1.2 Visualize Batch effects using tSNE ###################################################
f1       <- function(x) (IQR(x) > 1.5)
fselect  <- genefilter(eset.rma, filterfun(f1))
summary(fselect)
eset.sel <-eset.rma[fselect,]
set.seed(38)  
tsne_out <- Rtsne(t(eset.sel),dims = 2, initial_dims = 20, perplexity = 11,
                  theta = 0.5, check_duplicates = TRUE, pca = FALSE, max_iter = 1000,
                  verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData$IVIM_Color, pch=16, cex=1.3) # saved as: FullSagaSet_BH-SNE_RMA_IQR1.5_718genes_Perpl11_SEED38_over IVIMID.pdf
legend(-24,12.5, legend=unique(pData$IVIM_ID), col=unique(pData$IVIM_Color), pch=16, bty="n", cex=0.7)


### 1.1.3 COMBAT batch correction #############################################################
batch      <- pData$Batch                          
modcombat  <- model.matrix(~1, data=pData)         
eset.batch <- ComBat(dat=eset.rma, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=TRUE)   
dev.off()

### 1.1.4 Visualize COMBAT corrected training set ##############################################
f1       <- function(x) (IQR(x) > 0.9)
fselect  <- genefilter(eset.batch, filterfun(f1))
summary(fselect)
eset.sel <-eset.batch[fselect,]
set.seed(15)  
tsne_out <- Rtsne(t(eset.sel),dims = 2, initial_dims = 20, perplexity = 11,
                  theta = 0.5, check_duplicates = TRUE, pca = FALSE, max_iter = 1000,
                  verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData$IVIM_Color, pch=16, cex=1.3) # saved as: FullSagaSet_BH-SNE_RMA_COMBAT_IQR0.9_720genes_Perpl11_SEED15_over IVIMID
legend(-23,0, legend=unique(pData$IVIM_ID), col=unique(pData$IVIM_Color), pch=16, bty="n", cex=0.6)

plot(tsne_out$Y,col=pData$Design_Color, pch=16, cex=1.3) # saved as: FullSagaSet_BH-SNE_RMA_COMBAT_IQR0.9_720genes_Perpl11_SEED15_over Design
legend(-23,0, legend=unique(pData$Design), col=unique(pData$Design_Color), pch=16, bty="n", cex=0.6)


################################################################################################
#### 2. Divide into Training and Test set ######################################################
################################################################################################

matrix.training <- eset.batch[,c(8:75)] 
matrix.holdout  <- eset.batch[,-c(8:75)] 
pData.training  <- pData[c(8:75),]
pData.holdout   <- pData[-c(8:75),]

table(pData.training$Design)
table(pData.holdout$Design)

################################################################################################
#### 3. nonspecific Feature prefiltering  ######################################################
################################################################################################

f1       <- function(x) (IQR(x) > 0.78)    # 1188 genes
fselect  <- genefilter(matrix.training, filterfun(f1))
summary(fselect)
matrix.filtered <-matrix.training[fselect,]

################################################################################################
#### 4. RFE ####################################################################################
################################################################################################

################################################################################################
#### 4.1 prepare data for rfe  #################################################################
################################################################################################

matrix.train <- (t(matrix.filtered))
labels.train <- as.factor(pData.training$Label)

#### 4.2 rfe options  ##########################################################################
################################################################################################

caretFuncs$summary <- twoClassSummary
outerctrl          <- rfeControl(functions=caretFuncs,
                                 method = "cv", repeats =5, number = 10, 
                                 returnResamp="final", 
                                 verbose = TRUE, 
                                 allowParallel = TRUE)

innerctrl          <- trainControl(classProbs= TRUE, 
                                   summaryFunction = twoClassSummary)

##### 4.4. rfe SVMrad ##########################################################################
################################################################################################

## this takes around 40 min - 1 h on seven cores 
set.seed (5)
system.time(rfe.rad  <- rfe(matrix.train, labels.train, 
                sizes=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,40,50,60,70,80,90,100,150,200,300,400,500),
                rfeControl=outerctrl,
                metric = "ROC",
                method="svmRadial",
                tuneLength = 10,
                preProc = c("center", "scale"),
                trControl = innerctrl))

rfe.rad$results
rfe.rad$optVariables

# Figure 2i: saved to TrainingV2_1188features_rfe_SVMrad_tL10_SEED5_ROC over Variables_Zoom.pdf
trellis.par.set(caretTheme())
plot(rfe.rad, type = c("g", "o"),xlim = c(0,71))  

# RFE results for the different numbers of features:
write.table(rfe.rad$results, file = "TrainingV2_1188F_rfe_SVMrad_tL10_SEED5_Results.txt", sep="\t",col.names=NA)

# 17 optimal features 
optFeatures <- cbind(rfe.rad$optVariables, Annotation[rfe.rad$optVariables,1])
write.table(optFeatures, file = "TrainingV2_1188F_rfe_SVMrad_tL10_SEED5_OptVariables.txt", sep="\t",col.names=NA)

##### 4.5. Loop RFE-SVM for 20 SEEDs for 11-20 features ########################################
################################################################################################

results.rad.loop  <- NULL
features.rad.loop <- as.data.frame(matrix(data = rep(0,400), nrow = 20, ncol = 20, byrow = FALSE))

for(i in 1:20)  {set.seed (i)
  rfe.rad.loop <- rfe(matrix.train, labels.train, 
                      sizes=c(11,12,13,14,15,16,17,18,19,20),
                      rfeControl=outerctrl,
                      metric = "ROC",
                      method="svmRadial",
                      tuneLength = 10,
                      preProc = c("center", "scale"),
                      trControl = innerctrl)
  results.rad.loop <- rbind(results.rad.loop, rfe.rad.loop$results)
  features.rad.loop[c(1:length(rfe.rad.loop$optVariables)),i] <- rfe.rad.loop$optVariables
}

#### export results over the different SEEDs  ###################################################
results.rad.loop      <- subset(results.rad.loop,results.rad.loop$Variables != 1188)
results.rad.loop$SEED <- c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),
                           rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10),rep(16,10),rep(17,10),rep(18,10),rep(19,10),rep(20,10))
write.table(results.rad.loop, file = "TrainingV2_1188F_rfe_SVMrad_tL10_SEED1to20_Results.txt", sep="\t",col.names=NA)

#### most often selected features over the different seeds ########################################
BestFeatures        <- stack(features.rad.loop)[1]
BestFeaturesCounts  <- as.data.frame(table(BestFeatures$values)) # count how often each feature occurs
BestFeaturesCounts  <- arrange(BestFeaturesCounts, desc(Freq))   # arrange in descending order
BestFeaturesCounts  <- BestFeaturesCounts[-1,]                   # ditch the 0s
BestFeaturesCounts  <- cbind(BestFeaturesCounts, Annotation[as.character(BestFeaturesCounts$Var1),])
row.names(BestFeaturesCounts) <- BestFeaturesCounts[,1]
write.table(BestFeaturesCounts, file = "TrainingV2_1188F_rfe_SVMrad_tL10_SEED1to20_FeatureCounts.txt", sep="\t",col.names=NA)


################################################################################################
#### 5. Performance of RFE-selected featureSets  ###############################################
################################################################################################

#### 5.1. prepare data #########################################################################
################################################################################################

# choose all features from SVM-RFE that have been selected more than twice = Top19rfe
Top19rfe       <- subset(BestFeaturesCounts, BestFeaturesCounts$Freq>2)  

df.caret       <- as.data.frame(t(matrix.training[row.names(Top19rfe),]))
df.caret$label <- as.factor(pData.training$Label)
str(df.caret)

#### 5.2. 10CVn5 on TrainingSet Top19rfe #######################################################
################################################################################################
set.seed(2381)
intrain <- createMultiFolds(df.caret$label,k=10, times = 5)

final.frame <- NULL
svm_bestHP  <- NULL
system.time(for(i in 1:50) {
                            training   <- df.caret[intrain[[i]],]
                            testing    <- df.caret[-intrain[[i]],]
                            tunectrl   <- tune.control(sampling="cross", cross = 10, nrepeat = 3, best.model = TRUE, performances = TRUE)
                            gammalist  <- c(0.001,0.01,0.05,0.08,0.1,0.15,0.2,0.25,0.3,0.5,1,10)
                            costlist   <- c(0.5,1,2,3,4,5,10,100)
                            set.seed(3233)
                            svm_tune   <- tune(svm,label ~ ., data=training, kernel="radial", 
                                              ranges=list(gamma=gammalist,cost = costlist), tunecontrol = tunectrl)
                            svm_fit    <- svm(label ~ ., data=training, probability = TRUE, 
                                           gamma = svm_tune$best.parameters$gamma, cost = svm_tune$best.parameters$cost)
                            svm_bestHP <- rbind(svm_bestHP,svm_tune$best.parameters) 
                            # predict test set  
                            y          <- ncol(testing)-1
                            output     <- data.frame(attr(predict(svm_fit, testing[,c(1:y)], probability= T), "probabilities"))
                            output$prediction <- ifelse(output$transforming>0.5,"transforming","untransforming")   # decision is based on the calculated probs !
                            iteration  <- rep(i,nrow(testing))
                            samples    <- row.names(testing)
                            truelabels <- testing[,ncol(testing)]
                            testing.df <- data.frame(iteration,samples,output,truelabels= as.character(truelabels)) 
                            final.frame<- rbind(final.frame,testing.df)
                            })

sink("TrainingV2_ConfusionMatrix_10CVn5_Top19rfe.txt", append = TRUE)
confusionMatrix(as.factor(final.frame$prediction), as.factor(final.frame$truelabels))
sink()

#samples repeatedly wrong 
final.wrong     <- subset(final.frame,final.frame$truelabels != final.frame$prediction)
table(as.character(final.wrong$samples))                    

# best hyperparameters for the SVM
svm_bestHP$gc <- paste(svm_bestHP$gamma,svm_bestHP$cost)    
table(svm_bestHP$gc)



################################################################################################
#### 6. GENETIC ALGORITHM FOR THE BEST COMBINATION OF THE 19 PREDICTORS ########################
################################################################################################

### 6.1. perform GA  ###########################################################################

# prepare the training data 
matrix.GA.train <- as.data.frame(t(matrix.training[row.names(Top19rfe),]))
labels.GA.train <- as.factor(pData.training$Label)
svmGA           <- caretGA        # pre-made list of functions: interface to train

# conduct the GA search inside a resampling wrapper defined by index  
set.seed(104)
index <- createMultiFolds(labels.GA.train, times = 5)  

# outer control 
ga.ctrl <- gafsControl(functions = svmGA,
                       method = "repeatedcv",
                       repeats = 5,
                       index = index,
                       verbose = TRUE,
                       returnResamp = "all",
                       maximize = c(internal = TRUE,     # accuracy is default
                                    external = TRUE),
                       allowParallel = TRUE)

system.time(svm_ga.training<- gafs( x = matrix.GA.train,
                                    y = labels.GA.train,
                                    iters = 100,
                                    popSize = 50, pcrossover = 0.8, pmutation = 0.1, elite =3 ,  
                                    gafsControl = ga.ctrl,
                                    method = "svmLinear",
                                    tuneLength = 5,
                                    preProc = c("center", "scale"),
                                    trControl = trainControl(method = "repeatedcv",
                                                            repeats = 5,
                                                            allowParallel = FALSE)))

### 6.2. analyze results  ##################################################################
svm_ga.training
svm_ga.training$optIter
svm_ga.training$optVariables

#### Figure 2j:
plot(svm_ga.training) + theme_bw() # saved as: TrainingV2_GA_SVMlin_inputTop19rfe_100Gen_maxAccuracy.pdf

sink("TrainingV2_Results_GA_Top19rfe_100_50_0.8_0.1_3_maxAcc.txt", append = TRUE)
svm_ga.training
sink()

### 6.3 output of GA #######################################################################

performance.external <- svm_ga.training$external                  # external fitness in each iteration
performance.external <- arrange(performance.external, Iter) 

performance          <- svm_ga.training$ga$internal               # internal fitness in each iteration           
performance$external <- aggregate(performance.external$Accuracy, by=list(Iter=performance.external$Iter),mean)[,2]
bestGAsubsets <- subset(performance, performance$external > 0.97) # the three best iterations are 43,47,56: nearly identical performance

ga.subsets <- svm_ga.training$ga$subsets
Top9GA     <- BestFeaturesCounts[ga.subsets[[56]],]               
write.table(BestSubset, file = "Top9_I56_GA_Top19rfe_100_50_0.8_0.1_3_maxAcc.txt", sep="\t",col.names=NA)

#### 6.4. 10CVn5 on TrainingSet Top9GA #########################################################
################################################################################################
df.caret.GA       <- as.data.frame(t(matrix.training[row.names(Top9GA),]))
df.caret.GA$label <- as.factor(pData.training$Label)
str(df.caret.GA)

final.frame.GA <- NULL
svm_bestHP.GA  <- NULL
system.time(for(i in 1:50) {
        training   <- df.caret.GA[intrain[[i]],]
        testing    <- df.caret.GA[-intrain[[i]],]
        tunectrl   <- tune.control(sampling="cross", cross = 10, nrepeat = 3, best.model = TRUE, performances = TRUE)
        gammalist  <- c(0.01,0.05,0.08,0.1,0.15,0.2,0.25,0.3,0.5,1)
        costlist   <- c(0.5,1,2,3,4,5,10)
        set.seed(3233)
        svm_tune   <- tune(svm,label ~ ., data=training, kernel="radial", 
                           ranges=list(gamma=gammalist,cost = costlist), tunecontrol = tunectrl)
        svm_fit    <- svm(label ~ ., data=training, probability = TRUE, 
                           gamma = svm_tune$best.parameters$gamma, cost = svm_tune$best.parameters$cost)
        svm_bestHP.GA <- rbind(svm_bestHP.GA,svm_tune$best.parameters) 
        # predict test set  
        y          <- ncol(testing)-1
        output     <- data.frame(attr(predict(svm_fit, testing[,c(1:y)], probability= T), "probabilities"))
        output$prediction <- ifelse(output$transforming>0.5,"transforming","untransforming")   # decision is based on the calculated probs !
        iteration  <- rep(i,nrow(testing))
        samples    <- row.names(testing)
        truelabels <- testing[,ncol(testing)]
        testing.df <- data.frame(iteration,samples,output,truelabels= as.character(truelabels)) 
        final.frame.GA <- rbind(final.frame.GA,testing.df)
})

sink("TrainingV2_ConfusionMatrix_10CVn5_Top9GA.txt", append = TRUE)
confusionMatrix(as.factor(final.frame.GA$prediction), as.factor(final.frame.GA$truelabels))
sink()

#samples repeatedly wrong 
final.wrong.GA <- subset(final.frame.GA,final.frame.GA$truelabels != final.frame.GA$prediction)
table(as.character(final.wrong.GA$samples))                    

# best hyperparameters for the SVM
svm_bestHP.GA$gc <- paste(svm_bestHP.GA$gamma,svm_bestHP.GA$cost)    
table(svm_bestHP.GA$gc)


