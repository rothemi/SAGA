# SAGA
Surrogate Assay for Genotoxicity Assessment

# Using SAGA in R - first steps
In order to use the SAGA package some Bioconductor dependencies need to
be installed. First, the connection to the Bioconductor database must be
established. Then the required packages can be loaded independently. Please
note that these packages are only available over Bioconductor and cannot be
loaded using CRAN.

source("https://bioconductor/biocLite.R")
biocLite("Biobase")
biocLite("BioGenerics")
biocLite("affy")
biocLite("affyPLM")

The following packages need to be loaded as well. They can be installed via CRAN via the R console.
install.packages("limma")
install.packages("sva")
install.packages("e1071")
install.packages("phenoTest")
install.packages("gridExtra")
install.packages("bapred")

For more help on using the SAGA package see the Vignette:
vignette("saga_vignette")

or the individual help files for each function.
help(package = "saga", help_type = "html")

There is also a full example in the Vignette. 






