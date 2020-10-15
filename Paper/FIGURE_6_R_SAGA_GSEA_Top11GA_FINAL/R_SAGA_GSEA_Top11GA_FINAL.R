################################################################################################
############################### GSEA approach  #################################################
################################################################################################
library(phenoTest)
library(limma)
library(gridExtra)
library(pROC)
library(caret)

################################################################################################
#### 1. Data handling ##########################################################################
################################################################################################

### load the Top 11 SAGA features from SVM-GA for FullSet152: 
sets     <- read.delim("./SAGA_INBUILD/SAGA_INBUILD_Top11GA_GSEA.txt",header=FALSE,sep="\t",stringsAsFactors =FALSE,row.names = 1)
SAGA.CORE<- setNames(split(sets, seq(nrow(sets))), rownames(sets))   # GeneSets have to be stored in a list object

### 1.1. Read in files from user and loop over all batches separately ###########################
#################################################################################################

# Note: IVIM #170906 was excluded due to no MOCK control available for this assay
# IVIM ID 180523 used as is (no split and no mock duplicates) ==> 159 samples alltogether / 17 batches
# all original .txt files have to be in the directory and annotated by the SampleInformation File 

SIF        <- read.delim("SampleInformation.txt",row.names=1,header=TRUE,sep="\t", stringsAsFactors = F)
maxBatch   <- max(as.integer(SIF$Batch))   # how many assays / batches 
result.all <- NULL                         # dummy data.frame

################################################################################################
#### 2. SAGA-GSEA Loop over all batches ########################################################
################################################################################################

for(i in 1:maxBatch) {   
  #### 2.0 read in .txt of each batch  ###########################################################
  SIF.i <- SIF[SIF$Batch==i,]
  RAW.i <- read.maimages(files=SIF.i$Filename, path=".", source="agilent.median", green.only=T,
                         columns=list(G="gMedianSignal"), annotation=c("ProbeName", "GeneName"))
  colnames(RAW.i) <- row.names(SIF.i)
  #### 2.1. Normalize, average ###################################################################
  RMA.i <- normalizeBetweenArrays(RAW.i, method="quantile")            # quantil normalize
  RMA.i <- avereps(RMA.i,ID= RMA.i$genes$ProbeName)                    # average replicates to one value for each probe
  matrix.gsea <- RMA.i$E                                               # extract log2 expression values 
  
  #### 2.3. make ExpressionSet (Biobase) object ##################################################
  metadata  <- data.frame(labelDescription= rep(NA,dim(SIF.i)[2]),row.names=colnames(SIF.i))   # varMetadata: empty, but required 
  phenoData <- new("AnnotatedDataFrame",data=SIF.i, varMetadata=metadata)     # annotatedDataFrame for the annotation of the samples
  eset.gsea <- ExpressionSet(assayData = matrix.gsea, phenoData = phenoData)  # this is the ExpressionSet required for phenoTest
  
  #### 2.4. make ePheno object: contains the FCs associated with vector variable ##################
  vars2test   <- list(ordinal="Vector")    # Variables (here: Vectors) to test against MOCK, which is always = 1 in the SIF 
  epheno.gsea <- ExpressionPhenoTest(eset.gsea,vars2test,p.adjust.method='BH')
  
  #### 2.5 GSEA #################################################################################
  SAGA.GSEA <- gsea(x=epheno.gsea, gsets=SAGA.CORE ,B=2000,                  # calculate GSEA-scores based on the FC in the epheno object
                    center = TRUE, test = "perm", p.adjust.method='BH')
  
  result            <- summary(SAGA.GSEA)[,c(1,2,3,5,8)]                     # extract results (only NES- normalized enrichment scores)
  
  #### 2.6 output ###############################################################################
  Vector <- NULL    ### pull out the Vector index number from the result table                  
  for (a in 1:nrow(result)) {Vector[a] <- unlist(strsplit(as.character(result$variable[a]), ".", fixed = TRUE))[2] }
  result$Vector     <- Vector
  
  SIF.sub           <- SIF.i[SIF.i$Vector != 1, c(1,6) ]                     # pull out info of tested vectors                    
  SIF.sub$SampleID  <- row.names(SIF.sub)                   
  result.m          <- merge(SIF.sub,result, by.x="Vector", by.y = "Vector") # merge result with SIF for SampleIDs and FileNames
  result.all        <- rbind(result.all, result.m)                           # merge all results to one file
  write.table(result.m, file = paste("Results_SAGA.GSEA_Batch_",i,".txt",sep = ""), sep="\t",row.names = FALSE)
  # make pdf report / Figure 6b and 6c are from batch #2 / IVIM 150128
  pdf(file=paste("SAGA.GSEA_Batch_",i,".pdf",sep = ""),useDingbats = F,width = 10, height = 10)  
  grid.table(result.m,rows = NULL)
  plot(SAGA.GSEA,es.nes='nes',selGsets='SAGA.CORE')
  dev.off()
}


################################################################################################
#### 3. merge data  ############################################################################
################################################################################################

SIF$SampleID <- row.names(SIF)
result.red   <- result.all[, c(3:8)]
result.red$Prediction_GSEA <- ifelse(result.red$nes>1.0,"transforming","untransforming")    # best threshold w/o LTR.SF vector (derived from ROC below)

Results      <- merge(SIF, result.red, by.x = "SampleID",by.y = "SampleID", all.x = TRUE )
row.names(Results) <- Results$SampleID
Results      <- Results[,-1]
Results      <- Results[row.names(SIF),]   # same order as in SampleInformation File

# exclude Mock controls = 127 samples 
Results.woMock <- subset(Results, Results$Name != "A2_MOCK")
Results.woMock$Class <- ifelse(Results.woMock$TrueLabel == "transforming", "transforming","nontransforming")
write.table(Results.woMock, file = "Results_SAGA.GSEA_AllBatches.txt", sep="\t",row.names = FALSE)

Results.Rothe <- Results.woMock[Rothe_Vectors$V1,]
write.table(Results.Rothe, file = "Results_SAGA.GSEA_Rothe_Vectors for Figure 6.txt", sep="\t",row.names = FALSE)

################################################################################################
#### 4. results without LTR.RV.SF (positive control) - to determine best cut-off ###############
################################################################################################

Results.woMock <- subset(Results.woMock, !is.na(Results.woMock$TrueLabel)) # exclude samples w/o ground truth available ==> 112 samples

#### 4.1 filter for LTR.RV.SFFV ("pRSF91") ##################################################### 
SF91                   <- "A1_LTR.RV.SF.eGFP"
Results.SF91           <- subset(Results.woMock, Results.woMock$Name %in% SF91 ) 
min(Results.SF91$nes)  
max(Results.SF91$fdr) ## all vectors have nes > 1.77 & fdr <= 0.005 ==> LTR.SF is always correctly classified

#### 4.2 filter without LTR.RV.SFFV ("pRSF91") ################################################ 
Results.woSF91 <- subset(Results.woMock, !Results.woMock$Name %in% SF91 ) 

#### 4.3 ROC curve without LTR.RV.SFFV to determine best nes cutoff ############################ 
roc.woSF91 <- roc(Results.woSF91$Class,                    
                  Results.woSF91$nes,             
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#686868", lwd=4, grid=F, add=F,
                  print.auc=T,print.thres="best")   # best threshold = 1.0  ==> 88% sensitivity, 74% specificity

#### 4.4 ConfusionMatrix ######################################################################
sink("ConfusionMatrix_SAGA_GSEA_Top11GA_1.0Threshold_woSF91.txt", append = TRUE)
confusionMatrix(as.factor(Results.woSF91$Prediction_GSEA), as.factor(Results.woSF91$TrueLabel))
sink()

################################################################################################
#### 5. results for complete dataset using nes = 1 as cutoff ###################################
################################################################################################

### ROC using cutoff from above for non-LTR vectors
ROC.GSEA.ALL <- roc(Results.woMock$Class,   # TrueLabel over               
                    Results.woMock$nes,     # numeric enrichment score       
                    percent=TRUE, levels=c("nontransforming","transforming"),
                    plot=T, auc.polygon=F, max.auc.polygon=F, col = "#CB4848", lwd=4, grid=F,
                    print.auc=T,print.thres=1.0 )   # Specificity 74 % Sensitivity 95.2 % for NES > 1.0 AUC 94.9%

roc.woSF91 <- roc(Results.woSF91$Class,                    
                  Results.woSF91$nes,             
                  percent=TRUE, levels=c("nontransforming","transforming"),
                  plot=T, auc.polygon=F, max.auc.polygon=F, col = "#686868", lwd=4, grid=F, add=T,
                  print.auc=T,print.thres=1)   # best threshold = 1.0  ==> 88% sensitivity, 74% specificity


#### confusion matrix for complete dataset using NES > 1.0 ######################################### 
sink("ConfusionMatrix_SAGA_GSEA_Top11GA_1.0Threshold_all.txt", append = TRUE)
confusionMatrix(as.factor(Results.woMock$Prediction_GSEA), as.factor(Results.woMock$TrueLabel))
sink()

################################################################################################
#### 5. results for complete dataset vs IVIM ###################################################
################################################################################################
IVIM  <- read.delim("IVIM_MTT_ROC.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE)
IVIM$Classification_Code <- as.factor(IVIM$Classification_Code)
IVIM$Prediction <- ifelse(IVIM$MTT_Score>=3,"transforming","untransforming")  
IVIM$TrueLabel  <- ifelse(IVIM$Classification_Code==1,"transforming","untransforming")
confusionMatrix(as.factor(IVIM$Prediction), as.factor(IVIM$TrueLabel))

### 2.6.2. IVIM vs SAGA AUROC (for all vectors) Figure 6g ###################################################################################
#############################################################################################################################################
ROC.GSEA.ALL <- roc(Results.woMock$Class,   # TrueLabel over               
                    Results.woMock$nes,     # numeric enrichment score       
                    percent=TRUE, levels=c("nontransforming","transforming"),
                    plot=T, auc.polygon=F, max.auc.polygon=F, col = "#CB4848", lwd=4, grid=F,
                    print.auc=T,print.thres=1.0 )   # Specificity 74 % Sensitivity 95.2 % for NES > 1.0 AUC 94.9%

IVIM.total <- roc(IVIM$Classification_Code,    # response vector (factor or character)
                  IVIM$MTT_Score,              # predictor vector (numeric)
                  percent=TRUE, smooth = F,
                  plot=TRUE, auc.polygon=F, max.auc.polygon=F,
                  col = "#8285BC", grid=F, lwd = 4, cex.lab=1, 
                  print.auc=T, print.thres = 3, add=T )

roc.test(ROC.GSEA.ALL, IVIM.total, alternative = "greater")   # p-value = 4.516e-07


#### 6.  Rank plot of IVIM confirmed samples ###################################################
################################################################################################

Results.woMock.o <- Results.woMock[order(Results.woMock$nes,decreasing = FALSE),]
Results.woMock.o <- subset(Results.woMock.o,!is.na(Results.woMock.o$Class))
plot(c(1:112),Results.woMock.o$nes, col = Results.woMock.o$Design_Color, pch=16, cex=1.6)
abline(1.0,0, lty = 2, col = "grey", cex =2)
legend(80,-1, legend=unique(Results.woMock.o$Design), col=unique(Results.woMock.o$Design_Color), pch=16, bty="n", cex=0.8)













