# SAGA
**S**urrogate **A**ssay for **G**enotoxicity **A**ssessment

# How to use SAGA with provided test files
The user needs to download the following files into one single folder:

* SAGA script [“R_Package_SAGA_V6.R” ](/R_Package_SAGA_V6.R)
* [“SampleInformation.txt”](/SampleInformation.txt) file
* The folder [“SAGA_INBUILD”](/SAGA_INBUILD)
* The sample files (to test SAGA, download test files [here](https://www.dropbox.com/sh/v0xjkgibxq8btgr/AABUx5l0e0qcVHtGqegCA79ca?dl=0) (~800 MB))

Before using SAGA, the user needs to specify this new folder in R by the following command:

setwd("c:/name of your folder")

Now, the script can be sourced or executed row by row. By following this instruction, the test files can be analysed. Alternatively, the user can provide his own sample files.

(Important: the Sampleinformation.txt must be adapted in that case!).

A detailed description of all functions used in the R-script, instructions for creating a SampleInformation file and how to interpret the SAGA results can be found in the [SAGA_Vignette](./SAGA_Vignette.pdf).

Here, the [SAGA results](./SAGA_test_files_Results) for the provided test files can be found.
