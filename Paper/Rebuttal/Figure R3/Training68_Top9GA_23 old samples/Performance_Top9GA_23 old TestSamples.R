#### 1. Performance of 9GA Classifier on th 23 original test samples ##########################
###############################################################################################

# 9 predictors were derived by SVM-rfe and GA from training set comprised of 68 samples (IVIM ID: 150128,150304,150318,151014,160210,160413,160525,160706,161102)
# results were aggregated from the hold-out set (IVIM ID: 120411,170125,171102) = 23 samples predicted by SVM trained on the 68 training samples and 9 predictors

library(caret)
library(pROC)

Results_23 <- read.delim("Predictions_SAGA91_Top9GA_120411_170125_171102.txt",row.names=1,header=TRUE,sep="\t",stringsAsFactors =FALSE)  # annotation of training samples

#### 1.1 Confusion matrix for SVM #############################################################  
###############################################################################################
sink("ConfusionMatrix_Top9GA_23TestSamples.txt", append = TRUE)
confusionMatrix(as.factor(Results_23$Prediction.SVM.Caret), as.factor(Results_23$TrueLabel))
sink()

#### 1.2 ROC on probability "transforming" for SVM ############################################
###############################################################################################
Results_23$Class <- as.factor(ifelse(Results_23$TrueLabel == "transforming","transforming","nontransforming"))

pdf(file="ROC_SAGA_Top9GA_23TestSamples.pdf",useDingbats = F,width = 5, height = 5)  
roc1 <- roc(Results_23$Class,                    # response vector (factor or character)
            Results_23$transforming,             # predictor vector (numeric)
            percent=TRUE, levels=c("nontransforming","transforming"),
            plot=T, auc.polygon=F, max.auc.polygon=F, col = "#000000", grid=F,
            print.auc=T,print.thres=F)
dev.off()


