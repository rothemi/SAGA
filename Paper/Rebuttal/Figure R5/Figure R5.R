###### FIGURE 2A: Gene Expression changes induced by LTR.SF vectors ###########################
###############################################################################################
library(gplots)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(UsingR)
library(limma)
library(Rtsne)
library(BiocParallel)

bpparam <- MulticoreParam(workers = 8)
bpparam

################################################################################################
#### 1. Data handling & Preprocessing using EList objects ######################################
################################################################################################

pData      <- read.delim("TARGETS_Fig2a_FINAL.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
Annotation <- read.delim("Annotation_SAGA_FINAL_KNOWN_20181128.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)

#### 1.1 read in .txt files ####################################################################
################################################################################################

eset.raw  <- read.maimages(files=pData$Filename, path=".", source="agilent.median", green.only=T,
                                          columns=list(G="gMedianSignal"), annotation=c("ProbeName"))
colnames(eset.raw)  <- row.names(pData)

boxplot(log2(eset.raw$E),col=pData$Color,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       
RAW <- as.data.frame(eset.raw$E)
simple.densityplot(log2(RAW[,c(1:7)]))

#### 1.2. Preprocessing & Subsetting  ##########################################################
################################################################################################

RMA      <- normalizeBetweenArrays(eset.raw,method="quantile")    # quantile normalization
RMA      <- avereps(RMA, ID= RMA$genes$ProbeName)                 # average quadruplicates
eset.rma <- RMA$E                                                 # extract log2 expression values 

boxplot(eset.rma,col=pData$Color,names=pData$Name,boxwex=0.6,cex.axis=0.5,las=2,outline=FALSE)       
RMA.df <- as.data.frame(eset.rma)
simple.densityplot(RMA.df)

## subset for 36,226 annotated genes ############################################################ 
eset.rma <- eset.rma[row.names(Annotation),]
eset     <- cbind(eset.rma,Annotation)
write.table(eset,file="ESET_RMA_120411_FINAL.txt", sep="\t",col.names=NA)

################################################################################################
#### 3. unsupervised analysis ##################################################################
################################################################################################

#### 3.1 BH-SNE: Figure 2a  ####################################################################
################################################################################################

### define function for BH-SNE on eset without prior PCA:  
Kullback <- function(i) {set.seed(i)  
                        tsne_out <- Rtsne(t(eset.rma), dims = 2, perplexity = 2, theta = 0.5, pca = FALSE, 
                                    max_iter = 1000, verbose = FALSE, is_distance = FALSE)
                        min(tsne_out$itercosts)}   # return final value of Kullback-Leibler divergence

system.time(KL <- bplapply(seq_len(1000), Kullback, BPPARAM = bpparam))  # run 1000 times
KL <- unlist(KL, use.names=FALSE)
seed.index <- which(KL==min(KL))                                         # get the iteration with the minimum Kullback-Leibler divergence

set.seed(seed.index)                                                     # Set seed with minimal KL
tsne_out <- Rtsne(t(eset.rma), dims = 2, perplexity = 2, theta = 0.5, pca = FALSE, 
                  max_iter = 1000, verbose = FALSE, is_distance = FALSE)
plot(tsne_out$Y,col=pData$Color, pch=16, cex=2)                          
legend(-40,70, legend=unique(pData$Vector), col=unique(pData$Color), pch=16, bty="n",cex=1)


### output Figure 2a: 
pData  <- cbind(as.data.frame(tsne_out$Y),pData)

ggplot(pData,aes(x=V1,y=V2, colour = as.factor(Vector), shape = Immortalized)) +
     geom_point(size=7) + 
     theme_bw(base_size = 12, base_family = "") +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
     scale_colour_manual(values = unique(pData$Vector_Color)) + 
     scale_x_reverse() # rotate plot to make the depiction similar to Figure 2c


#### 3.2 Heatmap Figure 2b / R5 ################################################################
################################################################################################
hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(256)
hmcol <- rev(hmcol)

IQRs  <- apply(eset.rma, 1, function(x) IQR(x))   # define filter function: IQR  
q     <- quantile(IQRs, probs = 0.99)             # select Top 1% of genes based on variance / For Figure R5 select 1%-10%
index <- which(IQRs > q)   

eset.sel <-eset.rma[index,]
heatmap.2(eset.sel, col = hmcol, cexRow=0.7, cexCol=0.6,
          scale="row", trace="none", ColSideColors=pData$Color, labRow = NULL, labCol=pData$Name )

