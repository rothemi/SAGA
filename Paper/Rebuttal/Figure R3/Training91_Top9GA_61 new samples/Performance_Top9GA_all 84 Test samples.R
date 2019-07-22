#### 1. Performance of 9GA Classifier on all 84 test samples ##################################
###############################################################################################

# 9 predictors were derived by SVM-rfe and GA from training set comprised of 68 samples (IVIM ID: 150128,150304,150318,151014,160210,160413,160525,160706,161102)
# results were aggregated from the hold-out set (IVIM ID: 120411,170125,171102) = 23 samples predicted by SVM trained on the 68 training samples and 9 predictors
# and from newly generated SAGA Assays (IVIM ID: 171115,180110,180131,180620,180523,180801,180822 = 61 samples, predicted with SVMrad trained on the complete (old) training set of 91 samples and 9 predictors

Results_84 <- read.delim("Results_SAGA_all84TestSamples.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)  # annotation of training samples

#### 1.1 Confusion matrix for SVM #############################################################  
###############################################################################################
sink("ConfusionMatrix_Top9GA_all84TestSamples.txt", append = TRUE)
confusionMatrix(as.factor(Results_84$Prediction.SVM.Caret), as.factor(Results_84$TrueLabel))
sink()

#### 1.2 ROC on probability "transforming" for SVM ############################################
###############################################################################################
Results_84$Class <- as.factor(ifelse(Results_84$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_SAGA_Top9GA_all84TestSamples.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Results_84$Class,                    # response vector (factor or character)
            Results_84$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=F)
dev.off()

#### 1.3 Confusion matrix for GSEA ############################################################  
###############################################################################################
Results_GSEA <- subset(Results_84, !is.na(Results_84$pred.class.GSEA))   # all mock samples are NA

sink("ConfusionMatrix_GSEA.CORE18_65TestSamples.txt", append = TRUE)
confusionMatrix(as.factor(Results_GSEA$pred.class.GSEA), as.factor(Results_GSEA$TrueLabel))
sink()

#### 1.2 ROC on GSEA ##########################################################################
###############################################################################################

pdf(file="ROC_GSEA.CORE18_65TestSamples.pdf",useDingbats = F,width = 6, height = 5)  
roc1 <- roc(Results_GSEA$Class,                    # response vector (factor or character)
            Results_GSEA$nes,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres="best")
dev.off()

