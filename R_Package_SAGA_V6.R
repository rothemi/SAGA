library(RColorBrewer)
library(limma)
library(sva)
library(e1071)
library(bapred)   

################################################################################################
#### 1. Data handling ##########################################################################
################################################################################################

### the following files should be build into the R package as data files:
pData       <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_75Samples.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation  <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_Annotation.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
Top12       <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_Top9_GA.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1) # new in V6: Top9_GA
SAGA_Data   <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_75Data_AVE.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE)
SAGA_RAW    <- as.matrix(SAGA_Data[,-1])     
row.names(SAGA_RAW) <- SAGA_Data$PROBE_ID

### Read in files from user: assume that they are in the working directory #####################
################################################################################################
SIF                     <- read.delim("SampleInformation.txt",row.names=1,header=TRUE,sep="\t", stringsAsFactors = F) 
pData.user              <- read.delim("SampleInformation.txt",row.names=1,header=TRUE,sep="\t", stringsAsFactors = F) 
pData.user$Batch        <- pData.user$Batch + 8                  
pData.user$IVIM_Color   <- rep("#000000", nrow(pData.user))      
pData.user$Design_Color <- rep("#000000", nrow(pData.user))      

eset.user <- read.maimages(files=pData.user$Filename, path=".", source="agilent.median", green.only=T,
                       columns=list(G="gMedianSignal"), annotation=c("ProbeName", "GeneName"))

colnames(eset.user)    <- row.names(pData.user)   
matrix.user            <- eset.user$E                    
row.names(matrix.user) <- eset.user$genes$ProbeName     
matrix.user            <- avereps(matrix.user, ID= row.names(matrix.user)) 

stopifnot(all(row.names(SAGA_RAW) %in% row.names(matrix.user))) 
matrix.user <- matrix.user[row.names(SAGA_RAW),] 


### Make joint sample information file #########################################################
################################################################################################
stopifnot(all(colnames(pData.user) == colnames(pData)))   
pData.joint <- rbind(pData,pData.user)
boxplot(log2(cbind(SAGA_RAW,matrix.user)),col=pData.joint$IVIM_Color, names=pData.joint$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

################################################################################################
#### 2.  Addon Quantile Normalization   ########################################################
###################################################################################ä############
qunorm.SAGA  <- qunormtrain(t(SAGA_RAW))                          
matrix.SAGA  <- log2(t(qunorm.SAGA$xnorm))                        
matrix.user  <- log2(t(qunormaddon(qunorm.SAGA, t(matrix.user)))) 
boxplot(cbind(matrix.SAGA,matrix.user),col=pData.joint$IVIM_Color,names=pData.joint$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

###############################################################################################
#### 3. Addon COMBAT  #########################################################################
###############################################################################################
batch.SAGA   <- as.factor(pData$Batch)                                            
combat.SAGA  <- combatba(t(matrix.SAGA), batch = batch.SAGA)   
matrix.SAGA  <- t(combat.SAGA$xadj)                            
colnames(matrix.SAGA) <- row.names(pData)
matrix.user  <- t(combatbaaddon(combat.SAGA, t(matrix.user), batch = as.factor(SIF$Batch)))  

################################################################################################
#### 4. PCA on CORE Genes  #####################################################################
################################################################################################
matrix.Top12 <- cbind(matrix.SAGA, matrix.user)[row.names(Top12),]  
index        <- nrow(pData)+nrow(pData.user)   # new in V6: nrow(pData)               

pdf(file="PCA_SAGA.pdf",useDingbats = F,width = 6, height = 5)  
pca     <- prcomp(t(matrix.Top12))           
plot(pca$x, pch=16, col=pData.joint$Design_Color, cex=1, asp=1)
legend(1,-2, legend = c("transforming","mock","neutral","new samples"), col = unique(pData.joint$Design_Color), pch=16, bty="n", cex=0.8)
text(pca$x[c((nrow(pData)+1):index),c(1:2)], labels=pData.user$Filename, cex= 0.3, pos=3, offset = 0.3) # new in V6: (nrow(pData)+1)
dev.off()

################################################################################################
#### 5. split into Prediction and Known Sets ###################################################
################################################################################################
matrix.train   <- t(matrix.SAGA[row.names(Top12),]) 
labels.train   <- as.factor(pData$Class)   
matrix.unknown <- t(matrix.user[row.names(Top12),])            

################################################################################################
#### 7. e1071 SVM training, tuning & prediction ################################################
################################################################################################

# new in V6: tune model and take best hyperparameters instead of using fixed gamma and cost:
tunectrl_1  <- tune.control(sampling="cross", cross = 10, nrepeat = 5, best.model = TRUE, performances = TRUE) # new in V6
gammalist_1 <- c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5)    # new in V6
costlist_1  <- c(0.1,0.5,1,2,3,4,5,10)                                                # new in V6      

set.seed (1234)                                                                    
svm_tune_1  <- tune(svm, matrix.train, labels.train, kernel="radial",                      # new in V6            
                  ranges=list(gamma=gammalist_1,cost = costlist_1), tunecontrol = tunectrl_1)  

cat("Tuning successful, Gamma=",svm_tune_1$best.parameters$gamma,", Cost=", svm_tune_1$best.parameters$cost,   # new in V6
    ", Accuracy 10fold CV n5 =", 1-svm_tune_1$best.performance)

model  <- svm(matrix.train,labels.train, probability = TRUE,type="C-classification",         # new in V6
                gamma = svm_tune_1$best.parameters$gamma, cost = svm_tune_1$best.parameters$cost)

output <- data.frame(attr(predict(model, matrix.unknown, probability= T), "probabilities"),
                     predict(model, matrix.unknown))
write.table(output, file = paste("Predictions_SAGA.SVMrad.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)

################################################################################################
#### 8. e1071 train new SVM  with grid search ##################################################
################################################################################################

#### 8.1 new split into training and unknown data ##############################################
################################################################################################

pData.known    <- subset(pData.joint, !is.na(pData.joint$Class))  
pData.unknown  <- subset(pData.joint, is.na(pData.joint$Class))   

matrix.known   <- t(matrix.Top12[,row.names(pData.known)])    
labels.known   <- as.factor(pData.known$Class)                
matrix.unknown <- t(matrix.Top12[,row.names(pData.unknown)])   

#### 8.2 tune and train SVM ####################################################################
################################################################################################
tunectrl  <- tune.control(sampling="cross", cross = 10, nrepeat = 5, best.model = TRUE, performances = TRUE)
gammalist <- c(0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5)   # changed in V4
costlist  <- c(0.5,1,2,3,4,5,10)                                    # changed in V4

set.seed (1)
svm_tune  <- tune(svm,matrix.known, labels.known, kernel="radial", 
                  ranges=list(gamma=gammalist,cost = costlist), tunecontrol = tunectrl)

cat("Tuning successful, Gamma=",svm_tune$best.parameters$gamma,", Cost=", svm_tune$best.parameters$cost,
    ", Accuracy 10fold CV n5 =", 1-svm_tune$best.performance)

svm_fit <- svm(matrix.known,labels.known, probability = TRUE,type="C-classification",  
               gamma = svm_tune$best.parameters$gamma, 
               cost = svm_tune$best.parameters$cost)

#### 8.3 e1071 SVM prediction    ################################################################
################################################################################################

prediction <- data.frame(attr(predict(svm_fit, matrix.unknown, probability= T), "probabilities"),
                         predict(svm_fit, matrix.unknown))

write.table(prediction, file = paste("Predictions_SAGA.SVMrad_tuned.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


################################################################################################
############################### 9. SAGA-GSEA  ##################################################
################################################################################################
library(phenoTest)
library(gridExtra)

################################################################################################
#### 1. Data handling ##########################################################################
################################################################################################

# all the following is new in V6: 
sets     <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_CORE_18RFE_ProbeID.txt",header=FALSE,sep="\t",stringsAsFactors =FALSE,row.names = 1)
SAGA.CORE<- setNames(split(sets, seq(nrow(sets))), rownames(sets))   # GeneSets have to be stored in a list object

### 1.1. Read in files from user and loop over all batches separately from here on ##############
#################################################################################################
maxBatch   <- max(as.integer(SIF$Batch))   # how many assays / batches 

for(i in 1:maxBatch) {   
  SIF.i <- SIF[SIF$Batch==i,]
  RAW.i <- read.maimages(files=SIF.i$Filename, path=".", source="agilent.median", green.only=T,
                         columns=list(G="gMedianSignal"), annotation=c("ProbeName", "GeneName"))
  colnames(RAW.i) <- row.names(SIF.i)
  #### 2.1. Normalize, average ###################################################################
  RMA.i <- normalizeBetweenArrays(RAW.i, method="quantile")                 # quantil normalize
  RMA.i <- avereps(RMA.i,ID= RMA.i$genes$ProbeName)                                # average replicates to one value for each probe
  matrix.gsea <- RMA.i$E                                               # extract log2 expression values 
  
  #### 2.3. make ExpressionSet (Biobase) object ##################################################
  metadata  <- data.frame(labelDescription= rep(NA,dim(SIF.i)[2]),row.names=colnames(SIF.i))   # varMetadata: empty, but required 
  phenoData <- new("AnnotatedDataFrame",data=SIF.i, varMetadata=metadata)     # annotatedDataFrame for the annotation of the samples
  eset.gsea <- ExpressionSet(assayData = matrix.gsea, phenoData = phenoData)  # this is the ExpressionSet required for phenoTest
  
  #### 2.4. make ePheno object: contains the FCs associated with Group variable ##################
  vars2test   <- list(ordinal="Group")    # Variables (here: Groups) to test against MOCK, which are always Group = 1 in the SIF 
  epheno.gsea <- ExpressionPhenoTest(eset.gsea,vars2test,p.adjust.method='BH')
  
  #### 2.5 GSEA #################################################################################
  SAGA.GSEA <- gsea(x=epheno.gsea, gsets=SAGA.CORE ,B=2000,                  # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH')
  
  result            <- summary(SAGA.GSEA)[,c(1,2,3,5,8)]                     # extract results (only NES- normalized enrichment scores)
  result$pred.class <- ifelse(result$nes>0,"transforming","nontransforming") # prediction based on NES
  
  #### 2.6 output ###############################################################################
  Group <- NULL    ### pull out the Group index number from the result table
  for (a in 1:nrow(result)) {Group[a] <- unlist(strsplit(as.character(result$variable[a]), ".", fixed = TRUE))[2] }
  result$Group     <- Group
  
  SIF.sub           <- SIF.i[SIF.i$Group != 1, c(3,4,1) ]                     # pull out info of tested Groups
  SIF.sub$SampleID  <- row.names(SIF.sub)
  result.m          <- merge(SIF.sub,result, by.x="Group", by.y = "Group") # merge result with SIF for SampleIDs and FileNames
  write.table(result.m, file = paste("Results_SAGA.GSEA_Batch_",i,".txt",sep = ""), sep="\t",row.names = FALSE)
  # make pdf report
  pdf(file=paste("SAGA.GSEA_Batch_",i,".pdf",sep = ""),useDingbats = F,width = 10, height = 10)
  grid.table(result.m,rows = NULL)
  plot(SAGA.GSEA,es.nes='nes',selGsets='SAGA.CORE')
  dev.off()
}







