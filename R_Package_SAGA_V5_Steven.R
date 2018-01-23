library(RColorBrewer)
library(limma)
library(sva)
library(e1071)
library(bapred)   # new in V5: library for addon normalization and batch correction

################################################################################################
#### 1. Data handling ##########################################################################
################################################################################################
setwd("C:\Users\Zool\R_data\MyPackages\SAGA\SAGA_INBUILD")
### the following files should be build into the R package as data files:
pData       <- read.delim("SAGA_INBUILD_Samples.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation  <- read.delim("SAGA_INBUILD_Annotation.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
Top12       <- read.delim("SAGA_INBUILD_Top12_Manual.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
SAGA_Data   <- read.delim("SAGA_INBUILD_Data_AVE.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE)
SAGA_RAW    <- as.matrix(SAGA_Data[,-1])
row.names(SAGA_RAW) <- SAGA_Data$PROBE_ID

### Read in files from user: assume that they are in the working directory #####################
################################################################################################
SIF                     <- read.delim("SampleInformation.txt",row.names=1,header=TRUE,sep="\t", stringsAsFactors = F) # moved up in V5: read in again and leave the batches as is
pData.user              <- read.delim("SampleInformation.txt",row.names=1,header=TRUE,sep="\t", stringsAsFactors = F)
pData.user$Batch        <- pData.user$Batch + 9
pData.user$IVIM_Color   <- rep("#000000", nrow(pData.user))
pData.user$Design_Color <- rep("#000000", nrow(pData.user))

eset.user <- read.maimages(files=pData.user$Filename, path=".", source="agilent.median", green.only=T,
                       columns=list(G="gMedianSignal"), annotation=c("ProbeName", "GeneName"))

colnames(eset.user)    <- row.names(pData.user)
matrix.user            <- eset.user$E
row.names(matrix.user) <- eset.user$genes$ProbeName
matrix.user            <- avereps(matrix.user, ID= row.names(matrix.user))

stopifnot(all(row.names(SAGA_RAW) %in% row.names(matrix.user)))
matrix.user <- matrix.user[row.names(SAGA_RAW),] # new in V5: leave test set as matrix.user and normalize separately


### Make joint sample information file #########################################################
################################################################################################
stopifnot(all(colnames(pData.user) == colnames(pData)))
pData.joint <- rbind(pData,pData.user)

# new in V5: no joint data set!
boxplot(log2(cbind(SAGA_RAW,matrix.user)),col=pData.joint$IVIM_Color, names=pData.joint$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)




################################################################################################
#### 2.  Addon Quantile Normalization   ########################################################
###################################################################################ä############
# new in V5:
qunorm.SAGA  <- qunormtrain(t(SAGA_RAW))    # make qunorm object containing normalized data and parameters for addon normalization
matrix.SAGA  <- log2(t(qunorm.SAGA$xnorm))  # extract normalized data
matrix.user  <- log2(t(qunormaddon(qunorm.SAGA, t(matrix.user)))) # addon normalization of test data

#### new in V5: QC after normalization - separate sets
boxplot(cbind(matrix.SAGA,matrix.user),col=pData.joint$IVIM_Color,names=pData.joint$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)


###############################################################################################
#### 3. Addon COMBAT  #########################################################################
###############################################################################################
# new in V5:
batch.SAGA   <- as.factor(pData$Batch)                         # batches of trainings data only
combat.SAGA  <- combatba(t(matrix.SAGA), batch = batch.SAGA)   # create COMBAT object containing batch corrected data and parameters for addon correction
matrix.user  <- t(combatbaaddon(combat.SAGA, t(matrix.user), batch = as.factor(SIF$Batch)))  # addon COMBAT of test data
matrix.SAGA  <- t(combat.SAGA$xadj)
colnames(matrix.SAGA) <- row.names(pData)

################################################################################################
#### 4. PCA on CORE Genes  #####################################################################
################################################################################################

## new in V5: make joint expression matrix here!
matrix.Top12 <- cbind(matrix.SAGA, matrix.user)[row.names(Top12),]
index        <- 83+nrow(pData.user)


# new in V5: export PCA to pdf
pdf(file="PCA_SAGA.pdf",useDingbats = F,width = 6, height = 5)
pca     <- prcomp(t(matrix.Top12))
plot(pca$x, pch=16, col=pData.joint$Design_Color, cex=1, asp=1)
legend(5,-2, legend = c("transforming","mock","neutral","new samples"), col = unique(pData.joint$Design_Color), pch=16, bty="n", cex=0.8)
text(pca$x[c(84:index),c(1:2)], labels=pData.user$Filename, cex= 0.3, pos=3, offset = 0.3) # changed into 84
dev.off()



################################################################################################
#### 5. split into Prediction and Known Sets ###################################################
################################################################################################

## new in V5: train and test data directly from the separate datasets- no indexing necessary
matrix.train   <- t(matrix.SAGA[row.names(Top12),]) # Changed in V5 to matrix.SAGA
labels.train   <- as.factor(pData$Class)
matrix.unknown <- t(matrix.user[row.names(Top12),]) # Changed in V5 to matrix.user

################################################################################################
#### 7. e1071 SVM prediction    ################################################################
################################################################################################

model <- svm(matrix.train, labels.train, kernel="radial", cost=2, gamma=0.08,cross = 10,
             type="C-classification", probability=TRUE)

output <- data.frame(attr(predict(model, matrix.unknown, probability= T), "probabilities"),
                     predict(model, matrix.unknown))

# new in V5: col.names = NA option ==> otherwise the colnames in the txt file are shifted to the left
write.table(output, file = paste("Predictions_SAGA.SVMrad_fixed.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


################################################################################################
#### 8. e1071 train new SVM  with grid search ##################################################
################################################################################################

#### 8.1 new split into training and unknown data ##############################################
################################################################################################

pData.known    <- subset(pData.joint, !is.na(pData.joint$Class))
pData.unknown  <- subset(pData.joint, is.na(pData.joint$Class))

matrix.known   <- t(matrix.Top12[,row.names(pData.known)])    # changed in V4 to "matrix.known" since "matrix.train has been used above already
labels.known   <- as.factor(pData.known$Class)                # changed in V4 to "labels.known"
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

# new in V5: col.names = NA option ==> otherwise the colnames in the txt file are shifted to the left
write.table(prediction, file = paste("Predictions_SAGA.SVMrad_tuned.txt",sep = ""), sep="\t",row.names = TRUE,col.names=NA)


################################################################################################
############################### 9. GSEA approach  ##############################################
################################################################################################
library(phenoTest)
library(gridExtra)

################################################################################################
#### 1. Data handling ##########################################################################
################################################################################################

### the following  should be build into the R package: genes (char vector); SAGA.CORE (list)

genes    <- unique(read.delim("./SAGA_INBUILD/Allgenes in GSEA-Sets for SAGA final.txt",header=FALSE,sep="\t",stringsAsFactors =FALSE)$V1)
sets     <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_CORESET.txt",header=FALSE,sep="\t",stringsAsFactors =FALSE,row.names = 1)
SAGA.CORE<- setNames(split(sets, seq(nrow(sets))), rownames(sets))   # GeneSets have to be stored in a list object

### 1.1. Read in files from user and loop over all batches separately from here on ##############
#################################################################################################
maxBatch   <- max(as.integer(SIF$Batch))   # how many assays / batches

for(i in 1:maxBatch) {
        SIF.i <- SIF[SIF$Batch==i,]                                # SampleInfo batchwise
        RAW.i <- eset.user[,row.names(SIF.i)]                      # split RAW data batchwise

        #### 2.1. Normalize, average ###################################################################
        RMA.i <- normalizeBetweenArrays(RAW.i, method="quantile")  # quantile normalize batchwise: log2 is done automatically since this is an RAWELIST object!
        RMA.i <- avereps(RMA.i,ID= RMA.i$genes$ProbeName)          # average replicates to one value for each probe
        matrix.rma          <- RMA.i$E                             # extract log2 expression values
        row.names(matrix.rma) <- toupper(RMA.i$genes$GeneName)                         # name rows according to Genes

        #### 2.2. GSEA-matrix ##########################################################################
        matrix.gsea     <- subset(matrix.rma, row.names(matrix.rma) %in% genes)    # take all genes with entries in GSEA-Sets
        matrix.gsea.SD  <- apply(matrix.gsea, 1, sd)                               # calculate the absolute SD for each gene
        matrix.gsea     <- matrix.gsea[order(matrix.gsea.SD,decreasing=TRUE),]     # order eset.gsea according to SD in descending order
        matrix.gsea     <- subset(matrix.gsea,!duplicated(row.names(matrix.gsea))) # throw out duplicates

        #### 2.3. make ExpressionSet (Biobase) S4 object for input into phenoTest ######################
        metadata  <- data.frame(labelDescription= rep(NA,dim(SIF.i)[2]),row.names=colnames(SIF.i))   # varMetadata: empty, but required
        phenoData <- new("AnnotatedDataFrame",data=SIF.i, varMetadata=metadata)     # annotatedDataFrame for the annotation of the samples
        eset.gsea <- ExpressionSet(assayData = matrix.gsea, phenoData = phenoData)  # this is the ExpressionSet required for phenoTest

        #### 2.4. make ePheno object: contains the FCs associated with Group variable ##################
        vars2test   <- list(ordinal="Group")                                        # Variables (here: Groups) to test against MOCK, which is always = 1 in the SIF
        epheno.gsea <- ExpressionPhenoTest(eset.gsea,vars2test,p.adjust.method='BH')# ePheno object for input into gsea function

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







