# Feature Selection for the SAGA-SVM classifier

The Github folder [SAGA_Input_Feature_Selection](./SAGA_Input_Feature_Selection)  contains the R-script [R_Classification_FeatureSelection_Paper.R](./SAGA_Input_Feature_Selection/R_Classification_FeatureSelection_Paper.R)
and all necessary files to run the script used in the paper to perform feature
selection for the SAGA classifier. It starts with the RAW data of 83 SAGA
samples which are annotated in the file “Targets_full_170906.txt” followed by
normalization and batch correction. The SVM-RFE and genetic algorithm search
procedures are implemented with the package “caret”. All output from the script
is stored in the folder [SAGA_Output_Feature_Selection](./SAGA_Output_Feature_Selection). Note that the SVM-RFE
loop will take around 20 hrs and the genetic algorithm run for 100 iterations
takes around 80 hrs to finish on an 6-8 core desktop computer. 

