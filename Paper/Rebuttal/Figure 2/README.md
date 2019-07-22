Figure R2 is the same as Figure 3 d-f.

The main code is described in [Feature Selection](Paper/Feature Selection/20190306_rfe_GA40_152_10TestSets_IQR0.8_FINAL_median.R)

Below is the code to generate Figure R2:
```
#######################################################
### 3. Plot resampling results vs TestSet results: Figure 3 d-e #############################################################################

df.results <- Results[c(1:10),]   # w/o FinalSet

### for full model / Figure 3d: #############################################################################
Accuracy_TestSet_full.median   <- median(df.results$Accuracy_TestSet_full)
ResamplingAccuracy_full.median <- median(df.results$ResamplingAccuracy_full)

ggplot(df.results, aes(Accuracy_TestSet_full,ResamplingAccuracy_full, label = Testset)) +
  geom_point(size=7) +
  geom_errorbarh(aes(xmin = CI_Test_full_lower,  xmax= CI_Test_full_upper, height = 0.003),colour = "black",size = 0.2) +
  geom_errorbar (aes(ymin = CI_full_lower, ymax = CI_full_upper, width = 0.003),colour = "black",size = 0.3) +
  geom_point(aes(x=Accuracy_TestSet_full.median, y=ResamplingAccuracy_full.median), colour="red", size = 7)+
  coord_fixed(ratio = 1, xlim=c(0.65,1),ylim=c(0.65,1)) +
  geom_abline(intercept = 0, slope = 1, color = "red",linetype = "dashed") + 
  geom_text(aes(label=Testset),hjust=0.5, vjust=0.5, colour = "white", size = 5,fontface="bold") +
  scale_y_continuous(breaks=seq(0.65,1, 0.05)) + 
  scale_x_continuous(breaks=seq(0.65,1, 0.05)) + 
  xlab("Test set accuracy" ) +
  ylab("Cross-validation accuracy") +
  theme_bw() +
  theme(axis.title.x = element_text(size=17),axis.title.y = element_text(size=17),
        axis.text = element_text(size=14, color ="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())


### for SVM-rfe / Figure 3e:  #############################################################################
Accuracy_TestSet_rfe.median   <- median(df.results$Accuracy_TestSet_rfe)
ResamplingAccuracy_rfe.median <- median(df.results$ResamplingAccuracy_rfe)

ggplot(df.results, aes(Accuracy_TestSet_rfe,ResamplingAccuracy_rfe, label = Testset)) +
  geom_point(size=7) +
  geom_errorbarh(aes(xmin = CI_Test_rfe_lower,  xmax= CI_Test_rfe_upper, height = 0.003),colour = "black",size = 0.2) +
  geom_errorbar (aes(ymin = CI_rfe_lower, ymax = CI_rfe_upper, width = 0.003),colour = "black",size = 0.3) +
  geom_point(aes(x=Accuracy_TestSet_rfe.median, y=ResamplingAccuracy_rfe.median), colour="red", size = 7)+
  coord_fixed(ratio = 1, xlim=c(0.65,1),ylim=c(0.65,1)) +
  geom_abline(intercept = 0, slope = 1, color = "red",linetype = "dashed") + 
  geom_text(aes(label=Testset),hjust=0.5, vjust=0.5, colour = "white", size = 5,fontface="bold") +
  scale_y_continuous(breaks=seq(0.65,1, 0.05)) + 
  scale_x_continuous(breaks=seq(0.65,1, 0.05)) + 
  xlab("Test set accuracy" ) +
  ylab("Cross-validation accuracy") +
  theme_bw() +
  theme(axis.title.x = element_text(size=17),axis.title.y = element_text(size=17),
        axis.text = element_text(size=14, color ="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())


### for SVM-GA / Figure 3f  #############################################################################
df.results.GA              <- df.results[c(2,3,6,7,8,10),]    # subset for the Training / TestSet splits in which GA was performed 
Accuracy_TestSet_GA.median   <- median(df.results.GA$Accuracy_TestSet_GA)
ResamplingAccuracy_GA.median <- median(df.results.GA$ResamplingAccuracy_GA)

ggplot(df.results.GA, aes(Accuracy_TestSet_GA,ResamplingAccuracy_GA, label = Testset)) +
  geom_point(size=7) +
  geom_errorbarh(aes(xmin = CI_Test_GA_lower,  xmax= CI_Test_GA_upper, height = 0.003),colour = "black",size = 0.2) +
  geom_errorbar (aes(ymin = CI_GA_lower, ymax = CI_GA_upper, width = 0.003),colour = "black",size = 0.3) +
  geom_point(aes(x=Accuracy_TestSet_GA.median, y=ResamplingAccuracy_GA.median), colour="red", size = 7)+
  coord_fixed(ratio = 1, xlim=c(0.65,1),ylim=c(0.65,1)) +
  geom_abline(intercept = 0, slope = 1, color = "red",linetype = "dashed") + 
  geom_text(aes(label=Testset),hjust=0.5, vjust=0.5, colour = "white", size = 5,fontface="bold") +
  scale_y_continuous(breaks=seq(0.65,1, 0.05)) + 
  scale_x_continuous(breaks=seq(0.65,1, 0.05)) + 
  xlab("Test set accuracy" ) +
  ylab("Cross-validation accuracy") +
  theme_bw() +
  theme(axis.title.x = element_text(size=17),axis.title.y = element_text(size=17),
        axis.text = element_text(size=14, color ="black"),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())
```
