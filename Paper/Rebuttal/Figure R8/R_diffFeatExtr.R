library(limma)
library(UsingR)
library(sva)
library(genefilter)
library(Rtsne)
library(BiocParallel)

bpparam <- MulticoreParam(workers = 10)
bpparam

# Note: The probe A_55_P2337033 was deleted from all datasets/annotations 
# due to cross-hybridization with EGFP

################################################################################################
#### 1. Read in data into EListRaw objects #####################################################
################################################################################################
pData      <- read.delim("Targets.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)
RAW.084956 <- read.maimages(files=pData$FileName, path=".", source="agilent.median", green.only=T,
                            columns=list(G="gMedianSignal"), annotation=c("ProbeName"))

eset.RAW            <- RAW.084956$E
row.names(eset.RAW) <- RAW.084956$genes$ProbeName  
colnames(eset.RAW)  <- row.names(pData)
 
## Figure R8: RAW data intensities ############################################################## 
boxplot(log2(eset.RAW), col=pData$Color, names=row.names(pData),boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)   #saved as: Boxplot_RAW_FullSagaSet115.pdf    

RAW.DF <- as.data.frame(log2(eset.RAW))
simple.densityplot(RAW.DF[,21:22], do.legend=TRUE, color = c("red","black"))  

ggplot(RAW.DF, aes(X6251, X6251_V12,)) +
  geom_point(size = 0.1)
  

 

