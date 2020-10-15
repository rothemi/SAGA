library(limma)
library(UsingR)
library(sva)
library(genefilter)
library(Rtsne)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(rgl)
library(BiocParallel)
bpparam <- MulticoreParam(workers = 10)
bpparam

# Note: The probe A_55_P2337033 was deleted from all datasets/annotations 
# due to cross-hybridization with EGFP
# Code for Figure 2a-d at the end of the script

################################################################################################
#### 1. Read in data into EListRaw objects #####################################################
################################################################################################

pData.full       <- read.delim("SAGA_Targets_FINAL_169.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
Annotation       <- read.delim("Annotation_SAGA_FINAL_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)       # all 39428 probes 
Annotation.known <- read.delim("Annotation_SAGA_FINAL_KNOWN_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1) # 36226 probes with Annotation available


#### 1.1 Read in SAGA Assays on Agilent-048306  ################################################
################################################################################################
files.048306 <- pData.full$Filename[1:64]
system.time(RAW.048306   <- read.maimages(files=files.048306, path=".", source="agilent.median", green.only=T,
                              columns=list(G="gMedianSignal"), annotation=c("ProbeName")))

colnames(RAW.048306) <- row.names(pData.full)[1:64]  

## subset for the 39,428 probes from Agilent ID026655  
index1             <- RAW.048306$genes$ProbeName %in% row.names(Annotation)
RAW.048306_Agilent <- RAW.048306[index1,]                                            # 157,712 probes / 4 = 39,428 quadruplicate probes
RAW.048306_Agilent <- RAW.048306_Agilent[order(RAW.048306_Agilent$genes$ProbeName),] # order prior cbind
str(RAW.048306_Agilent)

#### 1.2 Read in SAGA Assays on Agilent-066423  ################################################
################################################################################################
files.066423 <- pData.full$Filename[65:67]
RAW.066423   <- read.maimages(files=files.066423, path=".", source="agilent.median", green.only=T,
                              columns=list(G="gMedianSignal"), annotation=c("ProbeName"))

colnames(RAW.066423) <- row.names(pData.full)[65:67]  

## subset for the 39428 probes from Agilent ID026655  
index2             <- RAW.066423$genes$ProbeName %in% row.names(Annotation)
RAW.066423_Agilent <- RAW.066423[index2,]  
RAW.066423_Agilent <- RAW.066423_Agilent[order(RAW.066423_Agilent$genes$ProbeName),] # order prior cbind
str(RAW.066423_Agilent)

#### 1.3 Read in SAGA Assays on Agilent-084107  ################################################
################################################################################################
files.084107 <- pData.full$Filename[68:91]
RAW.084107   <- read.maimages(files=files.084107, path=".", source="agilent.median", green.only=T,
                              columns=list(G="gMedianSignal"), annotation=c("ProbeName"))

colnames(RAW.084107) <- row.names(pData.full)[68:91]  

## subset for the 39,428 probes from Agilent ID026655  
index3             <- RAW.084107$genes$ProbeName %in% row.names(Annotation)
RAW.084107_Agilent <- RAW.084107[index3,]   
RAW.084107_Agilent <- RAW.084107_Agilent[order(RAW.084107_Agilent$genes$ProbeName),] # order prior cbind
str(RAW.084107_Agilent)

#### 1.4 Read in SAGA Assays on Agilent-084956  ################################################
################################################################################################
files.084956 <- pData.full$Filename[92:169]
RAW.084956   <- read.maimages(files=files.084956, path=".", source="agilent.median", green.only=T,
                              columns=list(G="gMedianSignal"), annotation=c("ProbeName"))

colnames(RAW.084956) <- row.names(pData.full)[92:169]  

## subset for the 39,428 probes from Agilent ID026655  
index4             <- RAW.084956$genes$ProbeName %in% row.names(Annotation)
RAW.084956_Agilent <- RAW.084956[index4,]  
RAW.084956_Agilent <- RAW.084956_Agilent[order(RAW.084956_Agilent$genes$ProbeName),] # order prior cbind
str(RAW.084956_Agilent)


################################################################################################
#### 2. combine datasets based on probes present on all arrays using cbind.EList ###############
################################################################################################
table(RAW.048306_Agilent$genes$ProbeName == RAW.066423_Agilent$genes$ProbeName) # check whether all probeNames match
table(RAW.048306_Agilent$genes$ProbeName == RAW.084107_Agilent$genes$ProbeName)
table(RAW.048306_Agilent$genes$ProbeName == RAW.084956_Agilent$genes$ProbeName)

SAGA_RAW    <- cbind.EList(RAW.048306_Agilent,RAW.066423_Agilent,RAW.084107_Agilent,RAW.084956_Agilent)
SAGA_matrix <- data.frame(ProbeID = RAW.048306_Agilent$genes$ProbeName, SAGA_RAW$E) # export
write.table(SAGA_matrix, file = "RAW_FullSagaSet169_quadruplicates_FINAL.txt", sep="\t",row.names = TRUE, col.names = NA)

eset.RAW     <- SAGA_RAW$E
eset.RAW.AVE <- avereps(eset.RAW, ID= SAGA_RAW$genes$ProbeName)  
write.table(eset.RAW.AVE, file = "RAW_FullSagaSet169_averaged_FINAL.txt", sep="\t",row.names = TRUE, col.names = NA)

eset.RAW.AVE.152 <- eset.RAW.AVE[,row.names(pData)]
write.table(eset.RAW.AVE.152, file = "SAGA_INBUILD_Data_AVE_152.txt", sep="\t",row.names = TRUE, col.names = NA)

## Supplementary Figure 2a: RAW data intensities ############################################### 
boxplot(log2(SAGA_RAW$E), col=pData.full$IVIM_Color, names=row.names(pData.full),boxwex=0.6,cex.axis=0.35,las=2,outline=FALSE)   #saved as: Boxplot_RAW_FullSagaSet115.pdf    

RAW.DF <- as.data.frame(log2(SAGA_RAW$E))
simple.densityplot(RAW.DF, do.legend=FALSE)  

################################################################################################
#### 3. Quantile Normalization and Averaging  ##################################################
################################################################################################
all(row.names(pData.full) == colnames(SAGA_RAW))   

SAGA_RMA <- normalizeBetweenArrays(SAGA_RAW,method="quantile")      # quantile normalization
SAGA_RMA <- avereps(SAGA_RMA, ID= SAGA_RMA$genes$ProbeName)         # average over ProbeIDs    
eset.rma <- SAGA_RMA$E                                              # export from EList to matrix
write.table(eset.rma, file = "ESET_RMA_FullSagaSet169_FINAL.txt", sep="\t",row.names = TRUE, col.names = NA)

## Supplementary Figure 2b: normalized and averaged intensities   ##############################
boxplot(eset.rma,col=pData.full$IVIM_Color,boxwex=0.6,cex.axis=0.35,names=row.names(pData.full), las=2,outline=FALSE) # saved as Boxplot_RMA_FullSagaSet115.pdf         

RMA.DF <- as.data.frame(eset.rma)
simple.densityplot(RMA.DF,do.legend=FALSE)                    

#### 3.2 t-SNE: Reveal batch effects in quantile normalized DataSet ############################
################################################################################################

#### 3.2.1 function for BH-SNE on complete eset.rma without prior PCA:##########################
Kullback <- function(i) {set.seed(i)  
                        tsne_out <- Rtsne(t(eset.rma), dims = 2, perplexity = 16,theta = 0.5, pca = FALSE, 
                                          max_iter = 1000, verbose = FALSE, is_distance = FALSE, check_duplicates = FALSE)
                        min(tsne_out$itercosts)}            # return final value of Kullback-Leibler divergence

#### 3.2.2 run Kullback 1000 times (in parallel on 10 cores) and record Kullback-Leibler divergence
KL <- bplapply(seq_len(1000), Kullback, BPPARAM = bpparam)  
KL <- unlist(KL, use.names=FALSE)    
seed.index <- which(KL==min(KL))  # get the iteration with the minimum Kullback-Leibler divergence: 520

#### 3.2.3 Figure 3a:  ######################################################################### 
set.seed(seed.index)  
tsne_out <- Rtsne(t(eset.rma),dims = 2, perplexity = 16,theta = 0.5, pca = FALSE, max_iter = 1000,
                  verbose = FALSE, is_distance = FALSE, check_duplicates = FALSE)
plot(tsne_out$Y,col=pData.full$IVIM_Color, pch=16, cex=1.4)    
legend(17.5,16, legend=unique(pData.full$IVIM_ID), col=unique(pData.full$IVIM_Color), pch=16, bty="n", cex=0.7)
text(tsne_out$Y, pData.full$Filename, cex=0.4, offset = 0.3, pos=3)


################################################################################################
#### 4. Batch correction   #####################################################################
################################################################################################

#### 4.1 Combat ################################################################################
################################################################################################
batch           <- pData.full$Batch                          
modcombat       <- model.matrix(~1, data=pData.full)         
eset.batch.full <- ComBat(dat=eset.rma, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=TRUE)
dev.off()
write.table(eset.batch.full, file = "ESET_RMA_COMBAT_FullSagaSet169_FINAL.txt", sep="\t",row.names = TRUE, col.names = NA)

#### 4.2 tSNE of COMBAT corrected ESET #######################################################
################################################################################################

#### 4.2.1 function for BH-SNE on complete eset.batch without prior PCA:########################
Kullback.batch <- function(i) {set.seed(i)  
                               tsne_out <- Rtsne(t(eset.batch.full), dims = 2, perplexity = 16,theta = 0.5, pca = FALSE, 
                                                 max_iter = 1000, verbose = FALSE, is_distance = FALSE, check_duplicates = FALSE)
                               min(tsne_out$itercosts)}            # return final value of Kullback-Leibler divergence

KLb <- bplapply(seq_len(1000), Kullback.batch, BPPARAM = bpparam) # run 1000 times (in parallel on 10 cores) and record Kullback-Leibler divergence
KLb <- unlist(KLb, use.names=FALSE)    
seed.index.batch <- which(KLb==min(KLb))    # get the iteration with the minimum Kullback-Leibler divergence:270

#### 4.2.2 Figure 3b:  #########################################################################
set.seed(seed.index.batch)  
tsne_out <- Rtsne(t(eset.batch.full),dims = 2, perplexity = 16, theta = 0.5,pca = FALSE, 
                  max_iter = 1000,verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData.full$IVIM_Color, pch=16, cex=1.4)      
legend(-13,0, legend=unique(pData.full$IVIM_ID), col=unique(pData.full$IVIM_Color), pch=16, bty="n", cex=0.5)

#### 4.2.3 visualize over Vector Design: Figure 3c  ###########################################
plot(tsne_out$Y,col=pData.full$Design_Color, pch=16, cex=1.4)  
legend(-13,-5, legend=unique(pData.full$Design), col=unique(pData.full$Design_Color), pch=16, bty="n", cex=1)


#### 4.3 Output of FINAL DataSet ###############################################################
################################################################################################

#### 4.3.1 remove two duplicate mock assays of IVIM ID 180523B (splitted preprocessing due to severe class imbalance)
pData.167      <- subset(pData.full, !pData.full$Filename %in% c("6374.1.txt","6379.1.txt"))  # remove 2 duplicate mock samples introduced for split-batch correction of #180523 due to severe class imbalance 
eset.batch.167 <- eset.batch.full[,row.names(pData.167)]  
write.table(eset.batch.167,file="ESET_RMA_COMBAT_FullSagaSet167_FINAL.txt",   sep="\t",col.names=NA)

#### 4.3.2 subset for the 36,226 annotated probes ###############################################
eset.batch.known.167 <- eset.batch.167[row.names(Annotation.known),]
write.table(eset.batch.known.167,file="ESET_RMA_COMBAT_KNOWN_FullSagaSet167_FINAL.txt",   sep="\t",col.names=NA)

#### 4.3.3 tSNE of FINAL DataSet / Figure 2e ###################################################
Kullback.FINAL <- function(i) {set.seed(i)  
                               tsne_out <- Rtsne(t(eset.batch.known.167), dims = 2, perplexity = 16,theta = 0.5, pca = FALSE, 
                                           max_iter = 1000, verbose = FALSE, is_distance = FALSE, check_duplicates = FALSE)
                                           min(tsne_out$itercosts)}            # return final value of Kullback-Leibler divergence

KLf <- bplapply(seq_len(1000), Kullback.FINAL, BPPARAM = bpparam) 
KLf <- unlist(KLf, use.names=FALSE)    
seed.index.final <- which(KLf==min(KLf)) # get the iteration with the minimum Kullback-Leibler divergence:476

set.seed(seed.index.final)  
tsne_out <- Rtsne(t(eset.batch.known.167),dims = 2, perplexity = 16, theta = 0.5,pca = FALSE, 
                  max_iter = 1000,verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData.167$Design_Color, pch=16, cex=2)  
legend(-13,-8, legend=unique(pData.167$Design), col=unique(pData.167$Design_Color), pch=16, bty="n", cex=1)


###############################################################################################
#### 5. limma - differentially expressed genes between vectors with known IVIM behaviour ######
###############################################################################################

#### 5.0. remove vectors with unknown ground truth (to few or inconclusive IVIM assays) ######
pData            <- subset(pData.167, !is.na(pData.167$Class)) 
eset.batch       <- eset.batch.167[,row.names(pData)]
eset.batch.known <- eset.batch.known.167[,row.names(pData)] 
table(pData$Design)  # 152 samples / 65 mutagenic 32 mock and 55 safe vectors

write.table(pData, file = "SAGA_Targets_FINAL_152.txt", sep="\t",row.names = TRUE, col.names = NA)

#### 5.1. design matrix #######################################################################
groups           <- factor(pData$Design)
design           <- model.matrix(~0+groups)
rownames(design) <- row.names(pData)
colnames(design) <- c("transforming","mock", "safe")
design

#### 5.2. contrast matrix ####################################################################
cont.matrix <- makeContrasts(c1 = transforming - mock,
                             c2 = safe - mock, 
                             c3 = transforming - (mock+safe)/2,
                             c4 = transforming - safe,
                             levels=design)
cont.matrix

#### 5.3. limma ##############################################################################
fit   <- lmFit(eset.batch.known,design)
fit   <- contrasts.fit(fit,contrasts=cont.matrix)
fit2  <- eBayes(fit)

#### 5.4. output: Supplementary Data 2 ######################################################

# Supplementary Data 2 tab 1
top.3 <- topTable(fit2, coef=3, n=Inf, adjust.method="BH", sort.by="logFC")
top.3 <- cbind(top.3,Annotation.known[row.names(top.3),])
write.table(top.3,file="Toplist_FullSagaSet152_Transforming vs Mock and Safe_FINAL.txt",   sep="\t",col.names=NA)

# Supplementary Data 2 tab 3: transforming - mock
top.1 <- topTable(fit2, coef=1, n=Inf, adjust.method="BH", sort.by="logFC")
top.1 <- cbind(top.1,Annotation.known[row.names(top.1),])
write.table(top.1,file="Toplist_FullSagaSet152_Transforming vs Mock_FINAL.txt",   sep="\t",col.names=NA)

# Supplementary Data 2 tab 4: safe vs mock
top.2 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", sort.by="logFC")
top.2 <- cbind(top.2,Annotation.known[row.names(top.2),])
write.table(top.2,file="Toplist_FullSagaSet152_Safe vs Mock_FINAL.txt",   sep="\t",col.names=NA)

# Supplementary Data 2 tab 5
top.4 <- topTable(fit2, coef=4, n=Inf, adjust.method="BH", sort.by="logFC")
top.4 <- cbind(top.4,Annotation.known[row.names(top.4),])
write.table(top.4,file="Toplist_FullSagaSet152_Transforming vs Safe_FINAL.txt",   sep="\t",col.names=NA)


###############################################################################################
#### 7. GSEA ##################################################################################
###############################################################################################

# Note: GSEA needs a collapsed matrix with only one probe for each gene
# since most genes have only one probe, this probe is selected
# in genes with more than one probe, one of the probes is often not functional/optimal ==> averaging quenches signal from functional probe
# ==> in case a gene has more than one probe, probes are not averaged, but the probe with the higher variance (SD) is selected

#### 7.1 import all gene symbols that occur in MSigDB GeneSets ################################
###############################################################################################
genes <- read.delim("Allgenes in C23567AS_20181205.txt",header=FALSE,sep="\t",stringsAsFactors =FALSE)
genes <- genes$V1   # 24,664 gene symbols that occur in at least one of the tested gene sets (MSigDB + own)

#### 7.2 make DataSet with one probe for each gene (select probe with highest variation) ######
###############################################################################################
eset           <- data.frame(eset.batch.known, SYMBOL=Annotation.known$SYMBOL)
eset.gsea      <- subset(eset, eset$SYMBOL %in% genes)                  # take all probes with entries in GSEA-Sets
matrix.gsea    <- as.matrix(eset.gsea[,c(1:152)])                       # export expression matrix
matrix.gsea.SD <- apply(matrix.gsea, 1, sd)                             # calculate the absolute SD for each probe
eset.gsea      <- eset.gsea[order(matrix.gsea.SD,decreasing=TRUE),]     # order eset.gsea according to SD in descending order
eset.gsea      <- subset(eset.gsea,!duplicated(eset.gsea$SYMBOL))       # throw out duplicates with lower SD
matrix.gsea    <- as.matrix(eset.gsea[,c(1:152)])                       # 15,376 unique Gencode annotated probes 
Annotation.GSEA<- Annotation.known[row.names(matrix.gsea),]

#### 7.2.1 export CHIP annotation file for GSEA ##############################################
GSEA.chip <- data.frame(ProbeSetID=row.names(matrix.gsea),GeneSymbol=Annotation.GSEA$SYMBOL,GeneTitle = Annotation.GSEA$GeneName_FINAL)
write.table(GSEA.chip,file="Annotation_GSEA.chip",sep="\t",row.names = FALSE) 

#### 7.2.2 export .gct files for GSEA Broad Tool #############################################
mock         <- subset(pData, pData$Design == "A2_Mock")      # select 32 mock samples
transforming <- subset(pData, pData$Design == "A1_Mutagenic") # select 65 transforming samples
safe         <- subset(pData, pData$Design == "A3_Safe")      # select 55 safe samples

matrix.gsea.transforming <- matrix.gsea[,row.names(transforming)]   
matrix.gsea.mock         <- matrix.gsea[,row.names(mock)]           
matrix.gsea.safe         <- matrix.gsea[,row.names(safe)]           

# data used as input for BROAD GSEA tool - output in Supplementary Data 3 tab 2:
df.TvM  <- data.frame(NAME=row.names(matrix.gsea),Description=Annotation.GSEA$SYMBOL,matrix.gsea.transforming,matrix.gsea.mock)
# data used as input for BROAD GSEA tool - output in Supplementary Data 3 tab 5:
df.TvS  <- data.frame(NAME=row.names(matrix.gsea),Description=Annotation.GSEA$SYMBOL,matrix.gsea.transforming,matrix.gsea.safe)
# data used as input for BROAD GSEA tool - output in Supplementary Data 3 tab 8:
df.SvM  <- data.frame(NAME=row.names(matrix.gsea),Description=Annotation.GSEA$SYMBOL,matrix.gsea.safe,matrix.gsea.mock)
# data used as input for BROAD GSEA tool - output in Supplementary Data 3 tab 11
df.TvSM <- data.frame(NAME=row.names(matrix.gsea),Description=Annotation.GSEA$SYMBOL,matrix.gsea.transforming,matrix.gsea.mock,matrix.gsea.safe)

numberofgenes   <- length(row.names(matrix.gsea))
numberofsamples <- dim(df.TvS)[2]-2

f <- file("GSEA_MATRIX_TvS.gct", "w")
writeLines(paste(c("#1.2", numberofgenes), c("", numberofsamples),sep="\t"),f) 
write.table(df.TvS, f,sep="\t",row.names = FALSE)
close(f)

###### 7.3 GSEA Plots Figure 2 f,g,h ##########################################################
###############################################################################################

# uses output from BROAD GSEA Tool on the .gct files created above 

# transforming vs mock (Figure 2f)
GSEA.TvM <- read.delim("GSEA_Results_TvM_STEM.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)

plot(log10(GSEA.TvM$FDR),GSEA.TvM$NES, ylab="Normalized enrichment score (NES)",
     xlab="FDR (log10)", col=GSEA.TvM$Color, pch=16, cex=3)        
segments(log10(0.10),-5,log10(0.10),5, lty = 2)

# safe vs mock (Figure 2g)
GSEA.SvM <- read.delim("GSEA_Results_SvM_STEM.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)

plot(log10(GSEA.SvM$FDR),GSEA.SvM$NES, ylab="Normalized enrichment score (NES)",
     xlab="FDR (log10)", col=GSEA.SvM$Color, pch=16, cex=3)        
segments(log10(0.10),-5,log10(0.10),5, lty = 2)

# transforming vs safe (Figure 2h)
GSEA.TvS <- read.delim("GSEA_Results_TvS_STEM.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
plot(log10(GSEA.TvS$FDR),GSEA.TvS$NES, ylab="Normalized enrichment score (NES)",
     xlab="FDR (log10)", col=GSEA.TvS$Color, pch=16, cex=3)        
segments(log10(0.10),-5,log10(0.10),5, lty = 2)


###############################################################################################
#### 8. self-contained GSEA with ROAST ########################################################
###############################################################################################

#### 8.1. read in genesets and convert genes to row indices of expression matrix ##############
###############################################################################################

sets        <- scan("c2.C5.hm.AML.STEM_v6.2.symbols.txt", what="", sep="\n")  # 10845 GeneSets
sets        <- scan("StemCellGeneSets_Fig2g_20181205.txt", what="", sep="\n") # 106 GeneSets

sets        <- strsplit(sets, "[[:space:]]+")                                 # convert to list of character vectors 
names(sets) <- sapply(sets, `[[`, 1)
sets        <- lapply(sets, `[`, -1)
sets        <- sets[lengths(sets)<301]  # no larger sets then 300
sets        <- sets[lengths(sets)>14]   # no smaller sets then 15  

id.vector        <- as.character(Annotation.GSEA$SYMBOL)  # named vector (ProbeID-SYMBOL)
names(id.vector) <- row.names(Annotation.GSEA)
idx              <- ids2indices(sets, id=id.vector)       # converts gene symbols in the geneset list into indices of the rows of matrix.GSEA

#### 8.2. ROAST: self contained geneset testing ###############################################
###############################################################################################

# for the >8000 sets from c2c5... set nrot = 5000

roast.TvM  <- mroast(matrix.gsea,idx,design,contrast=cont.matrix[,1], nrot = 50000,set.statistic="mean")    
roast.SvM  <- mroast(matrix.gsea,idx,design,contrast=cont.matrix[,2], nrot = 50000,set.statistic="mean")    
roast.TvSM <- mroast(matrix.gsea,idx,design,contrast=cont.matrix[,3], nrot = 50000,set.statistic="mean")   
roast.TvS  <- mroast(matrix.gsea,idx,design,contrast=cont.matrix[,4], nrot = 50000,set.statistic="mean")   

#Supplementary Data 3 tab 3
write.table(roast.TvM, file="ROAST_TvM_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)
#Supplementary Data 3 tab 6
write.table(roast.TvS, file="ROAST_TvS_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)
#Supplementary Data 3 tab 9
write.table(roast.SvM, file="ROAST_SvM_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)
#Supplementary Data 3 tab 12
write.table(roast.TvSM, file="ROAST_TvSM_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)

#### 8.3. CAMERA: competitive geneset testing #################################################
###############################################################################################
cam.TvM  <- camera(matrix.gsea,idx,design,contrast=cont.matrix[,1],inter.gene.cor=0.01)   
cam.SvM  <- camera(matrix.gsea,idx,design,contrast=cont.matrix[,2],inter.gene.cor=0.01) 
cam.TvSM <- camera(matrix.gsea,idx,design,contrast=cont.matrix[,3],inter.gene.cor=0.01)   
cam.TvS  <- camera(matrix.gsea,idx,design,contrast=cont.matrix[,4],inter.gene.cor=0.01)

#Supplementary Data 3 tab 4
write.table(cam.TvM, file="CAMERA_TvM_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)
#Supplementary Data 3 tab 7
write.table(cam.TvS, file="CAMERA_TvS_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)
#Supplementary Data 3 tab 10
write.table(cam.SvM, file="CAMERA_SvM_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)
#Supplementary Data 3 tab 13
write.table(cam.TvSM, file="CAMERA_TvSM_C2C5HMAMLSTEM.txt",sep="\t",col.names=NA)


################################################################################################
#### Figure 2 a & 2b: IVIM 120411 ##############################################################
################################################################################################
pData.120411  <- pData[pData$Batch==1,]               # select batch 1 / IVIM #120411 

#### 1.1 read in .txt files ####################################################################
################################################################################################
eset.raw.120411 <- read.maimages(files=pData.120411$Filename, path=".", source="agilent.median", green.only=T,
                                 columns=list(G="gMedianSignal"), annotation=c("ProbeName"))
colnames(eset.raw.120411)  <- row.names(pData.120411)
boxplot(log2(eset.raw.120411$E),col=pData.120411$IVIM_Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

#### 1.2. Preprocessing & Subsetting  ##########################################################
################################################################################################

RMA.120411      <- normalizeBetweenArrays(eset.raw.120411,method="quantile")    # quantile normalization
RMA.120411      <- avereps(RMA.120411, ID= RMA.120411$genes$ProbeName)                 # average quadruplicates
eset.rma.120411 <- RMA.120411$E                                                 # extract log2 expression values 

boxplot(eset.rma.120411,col=pData.120411$IVIM_Color,names=pData.120411$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

## subset for 36,226 annotated genes ############################################################ 
eset.rma.120411 <- eset.rma.120411[row.names(Annotation.known),]
eset.120411     <- cbind(eset.rma.120411,Annotation.known)
write.table(eset.120411,file="ESET_RMA_120411_FINAL.txt", sep="\t",col.names=NA)

################################################################################################
#### 2. unsupervised analysis ##################################################################
################################################################################################

#### 2.1 BH-SNE: Figure 2a  ####################################################################
################################################################################################

### define function for BH-SNE on eset without prior PCA
Kullback.120411 <- function(i) {set.seed(i)  
                                tsne_out <- Rtsne(t(eset.rma.120411), dims = 2, perplexity = 2, theta = 0.5, pca = FALSE, 
                                           max_iter = 1000, verbose = FALSE, is_distance = FALSE)
                                min(tsne_out$itercosts) # return final value of Kullback-Leibler divergence
                               }   

## run 100 times and record Kullback Leibler divergence for each run
KL.120411 <- bplapply(seq_len(1000), Kullback.120411, BPPARAM = bpparam)  # run 1000 times
KL.120411 <- unlist(KL.120411, use.names=FALSE)

## get the iteration with the minimum Kullback-Leibler divergence
seed.index.120411 <- which(KL.120411==min(KL.120411))  

## perform t-SNE
set.seed(seed.index.120411)                                                     
tsne_out.120411 <- Rtsne(t(eset.rma.120411), dims = 2, perplexity = 2, theta = 0.5, pca = FALSE, 
                         max_iter = 1000, verbose = FALSE, is_distance = FALSE)
### plot Figure 2a: 
pData.120411  <- cbind(as.data.frame(tsne_out.120411$Y),pData.120411)

ggplot(pData.120411,aes(x=V1,y=V2, colour = as.factor(Vector), shape = Immortalized)) +
  geom_point(size=7) + 
  theme_bw(base_size = 12, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_colour_manual(values = unique(pData$Vector_Color)) + 
  scale_x_reverse() 

#### 2.2 Heatmap Figure 2b #####################################################################
################################################################################################
hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(256)
hmcol <- rev(hmcol)

IQRs  <- apply(eset.rma.120411, 1, function(x) IQR(x))   # apply filter function IQR to eset.rma  
q     <- quantile(IQRs, probs = 0.99)                    # select Top 1% of genes based on variance
index.IQR <- which(IQRs > q)   

eset.sel.120411 <-eset.rma.120411[index.IQR,]
## Figure 2b: 
heatmap.2(eset.sel.120411, col = hmcol, cexRow=0.7, cexCol=0.6,
          scale="row", trace="none", ColSideColors=pData.120411$Design_Color, labRow = NULL, labCol=pData.120411$Name )

################################################################################################
#### Figure 2 c & d: IVIM 150128 ##############################################################
################################################################################################
pData.150128  <- pData[pData$Batch==2,]               # select batch 1 / IVIM #150128 

#### 1.1 read in .txt files ####################################################################
################################################################################################
eset.raw.150128 <- read.maimages(files=pData.150128$Filename, path=".", source="agilent.median", green.only=T,
                                 columns=list(G="gMedianSignal"), annotation=c("ProbeName"))
colnames(eset.raw.150128)  <- row.names(pData.150128)
boxplot(log2(eset.raw.150128$E),col=pData.150128$IVIM_Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

#### 1.2. Preprocessing & Subsetting  ##########################################################
################################################################################################
RMA.150128      <- normalizeBetweenArrays(eset.raw.150128,method="quantile")    # quantile normalization
RMA.150128      <- avereps(RMA.150128, ID= RMA.150128$genes$ProbeName)          # average quadruplicates
eset.rma.150128 <- RMA.150128$E                                                 # extract log2 expression values 
boxplot(eset.rma.150128,col=pData.150128$IVIM_Color,names=pData.150128$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       

## subset for 36,226 annotated genes ############################################################ 
eset.rma.150128 <- eset.rma.150128[row.names(Annotation.known),]
eset.150128     <- cbind(eset.rma.150128,Annotation.known)
write.table(eset.150128,file="ESET_RMA_150128_FINAL.txt", sep="\t",col.names=NA)

################################################################################################
#### 2. unsupervised analysis ##################################################################
################################################################################################

#### 2.1 BH-SNE: Figure 2c  ####################################################################
################################################################################################

### define function for BH-SNE on eset without prior PCA
Kullback.150128 <- function(i) {set.seed(i)  
  tsne_out <- Rtsne(t(eset.rma.150128), dims = 2, perplexity = 2, theta = 0.5, pca = FALSE, 
                    max_iter = 1000, verbose = FALSE, is_distance = FALSE)
  min(tsne_out$itercosts) # return final value of Kullback-Leibler divergence
}   

## run 1000 times and record Kullback Leibler divergence for each run
KL.150128 <- bplapply(seq_len(1000), Kullback.150128, BPPARAM = bpparam)  # run 1000 times
KL.150128 <- unlist(KL.150128, use.names=FALSE)

## get the iteration with the minimum Kullback-Leibler divergence
seed.index.150128 <- which(KL.150128==min(KL.150128))  

## perform t-SNE
set.seed(seed.index.150128)                                                     
tsne_out.150128 <- Rtsne(t(eset.rma.150128), dims = 2, perplexity = 2, theta = 0.5, pca = FALSE, 
                         max_iter = 1000, verbose = FALSE, is_distance = FALSE)
### plot Figure 2c: 
pData.150128  <- cbind(as.data.frame(tsne_out.150128$Y),pData.150128)

ggplot(pData.150128,aes(x=V1,y=V2, colour = as.factor(Vector), shape = Immortalized)) +
  geom_point(size=7) + 
  theme_bw(base_size = 12, base_family = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_colour_manual(values = unique(pData$Vector_Color))+ 
  scale_y_reverse() 

#### 2.2 Heatmap Figure 2d #####################################################################
################################################################################################
IQRs.150128      <- apply(eset.rma.150128, 1, function(x) IQR(x))   # apply filter function IQR to eset.rma  
q.150128         <- quantile(IQRs.150128, probs = 0.99)             # select Top 1% of genes based on variance
index.IQR.150128 <- which(IQRs.150128 > q.150128)   
eset.sel.150128  <- eset.rma.150128[index.IQR.150128,]

## Figure 2d: 
heatmap.2(eset.sel.150128, col = hmcol, cexRow=0.7, cexCol=0.6,
          scale="row", trace="none", ColSideColors=pData.150128$Design_Color, labRow = NULL, labCol=pData.150128$Name )





