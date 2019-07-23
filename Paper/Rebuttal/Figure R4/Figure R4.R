library(limma)
library(UsingR)
library(sva)
library(genefilter)
library(Rtsne)
library(BiocParallel)

bpparam <- MulticoreParam(workers = 10)
bpparam


################################################################################################
#### 1. Data  ##################################################################################
################################################################################################

# eset.rma:  39,428 x 169 quantile normalized and averaged intensitiy values (R_SAGA_ESET_169_FINAL)
# eset.batch.full: 39,428 x 169 quantile normalized, averaged and batch corrected intensitiy values (R_SAGA_ESET_169_FINAL)  
# pData.full: Annotation for the 169 SAGA samples

boxplot(eset.rma,col=pData.full$IVIM_Color,boxwex=0.6,cex.axis=0.35,names=row.names(pData.full), las=2,outline=FALSE) # saved as Boxplot_RMA_FullSagaSet115.pdf         
boxplot(eset.batch.full,col=pData.full$IVIM_Color,boxwex=0.6,cex.axis=0.35,names=row.names(pData.full), las=2,outline=FALSE) # saved as Boxplot_RMA_FullSagaSet115.pdf         


####### 2 t-SNE with different parameters for IQR and perplexitiy (for rebuttal) ###############
################################################################################################

#### nonspecific filtering using inter-quartile range ##########################################
f1       <- function(x) (IQR(x) > 0)          # range from 0...2
fselect  <- genefilter(eset.rma, filterfun(f1))
summary(fselect)
eset.sel <-eset.rma[fselect,]

set.seed(520)  
tsne_out <- Rtsne(t(eset.sel),dims = 2, perplexity = 3, # range from 3...30
                  theta = 0.5, pca = FALSE, max_iter = 1000,verbose = FALSE, is_distance = FALSE, check_duplicates = FALSE)
plot(tsne_out$Y,col=pData.full$IVIM_Color, pch=16, cex=1.4)    

####### 3 t-SNE with different parameters for IQR and perplexitiy (for rebuttal) on batch corrected set
################################################################################################

#### nonspecific filtering using inter-quartile range ##########################################
f1       <- function(x) (IQR(x) > 0)          # range from 0...2
fselect  <- genefilter(eset.batch.full, filterfun(f1))
summary(fselect)
eset.sel.batch <-eset.batch.full[fselect,]

set.seed(270)  
tsne_out <- Rtsne(t(eset.sel.batch),dims = 2, perplexity = 3, # range from 3...30
                  theta = 0.5, pca = FALSE, max_iter = 1000,verbose = FALSE, is_distance = FALSE, check_duplicates = FALSE)
plot(tsne_out$Y,col=pData.full$Design_Color, pch=16, cex=1.4)    

