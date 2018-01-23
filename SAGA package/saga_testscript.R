
# SAGA Script

library(saga)

samplepath     <- "somewhere useful"

# saga_gentargets; automatically generates an empty (!) sample information file
# Modify manually according to your requirements
#targets        <- saga_gentargets(samplepath)

################################################################################
### Wrapper function - all in one
### to use: uncomment function & execute
################################################################################
#mySAGAres      <- saga_wrapper(samplepath, showPCA=1, doGESEA=1)


################################################################################
### Or use these single functions
################################################################################
# saga_import
rawdata        <- saga_import(samplepath, showjoint=1)
eset.user      <- rawdata$eset.user
pData.user     <- rawdata$pData.user
pData.joint    <- rawdata$pData.joint
SIF            <- rawdata$SIF

# normalize saga data
normalized     <- saga_norm(rawdata$SAGA_RAW, rawdata$pData.joint, rawdata$matrix.user, normplot=1)
matrix.SAGA    <- normalized$matrix.SAGA
matrix.user    <- normalized$matrix.user

# remove batch effects
batchnorm      <- saga_batch(matrix.SAGA, matrix.user, SIF)
matrix.SAGA    <- batchnorm$matrix.SAGA
matrix.user    <- batchnorm$matrix.user

# filtering & model building
NN             <- saga_sampling(matrix.SAGA, matrix.user, pData.joint, pData.user, showPCA=1)
matrix.train   <- NN$matrix.train
labels.train   <- NN$labels.train
matrix.unknown <- NN$matrix.unknown
matrix.Top12   <- NN$matrix.Top12

# Array predictions with optimized SVM parameters (default settings)
output         <- saga_predict(matrix.train, labels.train, matrix.unknown, writeFile=1)
output

# Predictions with own neg/pos controls and grid optimization; may be used without controls
optPred        <- saga_optmodel(pData.joint, matrix.Top12, showbest=0, writeFile=1)
optPred

# GESEA
gesea_results  <- saga_gesea(samplepath, SIF)



############################
######## HELP Files ########
############################
# to see the help overview
help(package = "saga", help_type = "html")

# to see the vignette
vignette("saga_vignette")





