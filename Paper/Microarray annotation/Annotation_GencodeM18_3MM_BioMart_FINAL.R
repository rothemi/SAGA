################################################################################################
##### 1. Annotation of Agilent Arrays ID 026655  ###############################################
################################################################################################
library(Biostrings)
library(BiocParallel)
bpparam <- MulticoreParam(workers = 11)
bpparam

# Note: The probe A_55_P2337033  was deleted from all datasets/annotations 
# due to cross-hybridization with EGFP

#### 1.1 original annotation for ID026655 supplied by Agilent ##################################
################################################################################################
Annotation.INBUILD  <- read.delim("SAGA_INBUILD_Annotation.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
table(Annotation.INBUILD$Description =="Unknown")                       # 35492 of 39428 probes annotated

#### 1.2 Annotation for ID026655 downloaded from Agilent eArray 2018-07-12 ############## 
################################################################################################
Annotation.ID026655 <- read.delim("Agilent_Annotation_026655_D_AA_20171030_wo.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1, na.strings = "")
table(row.names(Annotation.ID026655) == row.names(Annotation.INBUILD))  # all Probes are ident
table(is.na(Annotation.ID026655$GeneSymbol))                            # 33706 of 39428 probes annotated with GeneSymbols
table(is.na(Annotation.ID026655$Description))                           # 35126 of 39428 probes annotated with Gene Description

################################################################################################
#### 2. De novo mapping of probe sequences to Gencode M18 (GRCm38.p6) Transcriptome ############
################################################################################################

####  2.1. Data ################################################################################
################################################################################################

# fasta files downloaded from https://www.gencodegenes.org/mouse/release_M18.html
Gencode   <- readDNAStringSet("gencode.vM18.transcripts.fa", format="fasta",use.names=TRUE)      # all 136535 Gencode M17 transcripts
Coding    <- readDNAStringSet("gencode.vM18.pc_transcripts.fa", format="fasta",use.names=TRUE)   # all 64732 coding transcripts

# split set to prevent memory exhaust
Annotation.1  <- Annotation.ID026655[1:20000,]
Annotation.2  <- Annotation.ID026655[20001:39428,]

probes.1      <- DNAStringSet(Annotation.1$Sequence)            # put probeseqs into DNAStringSet container
pdict.1       <- PDict(probes.1, tb.start = 25, tb.end = 35)    # convert into PDict library; define trusted band as 10 bp in the middle of the probe where no mismatches are allowed

probes.2      <- DNAStringSet(Annotation.2$Sequence)            
pdict.2       <- PDict(probes.2, tb.start = 25, tb.end = 35)    

################################################################################################
####  2.2. 1st priority: Mapping to the coding transcriptome ###################################
################################################################################################

# note: to prevent memory exhaust first do everything for Annotation.1, then remove pc.inc.matrix.1 
# and start again with Annotation.2 / change entries in function etc. accordingly and remove pc.inc.matrix.2

##### 2.2.1 mapping: incidence matrix with a row for each probe and col for each transcript ####
################################################################################################
pc.inc.matrix.1    <- vcountPDict(pdict.1, Coding, max.mismatch = 3)   
pc.inc.matrix.2    <- vcountPDict(pdict.2, Coding, max.mismatch = 3)    

##### 2.2.2 extract names of matched coding transcripts from fasta file ########################
################################################################################################

#                                                        1                   2                    3                    4                  5         6   7      8          9            10                    
# Namestring in Coding fasta file has format: ENSMUST00000070533.4|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|AC157543.1-001|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634|"
# that is:                                     ensembl_transcript     ensembl_gene          havanna_gene       havanna_transcript  transcript_name Symbol       5'UTR       CDS       3'UTR              

getGeneID.pc            <- function(i) {ifelse(max(pc.inc.matrix.2[i,])==0,NA,
                                               unlist(strsplit(names(Coding[which(pc.inc.matrix.2[i,]>0)[1]]), "|", fixed = TRUE))[2])} # extract 2nd element from Namestring = Ensembl Gene ID  
pc.ensembl_gene_ids     <- bplapply(seq_len(length(pdict.2)), getGeneID.pc, BPPARAM = bpparam)
pc.ensembl_gene_ids     <- unlist(pc.ensembl_gene_ids, use.names=FALSE)    

Annotation_pc_Gencode_1 <- data.frame(row.names = row.names(Annotation.1),pc.ensembl_gene_ids = pc.ensembl_gene_ids)
rm(pc.inc.matrix.1) 
rm(ensembl_gene_ids) ### repeat with pc.inc.matrix.2

Annotation_pc_Gencode_2 <- data.frame(row.names = row.names(Annotation.2),pc.ensembl_gene_ids = pc.ensembl_gene_ids)
rm(pc.inc.matrix.2)

Annotation_pc_Gencode   <- rbind(Annotation_pc_Gencode_1,Annotation_pc_Gencode_2)
write.table(Annotation_pc_Gencode, file = "Annotation_pc_Gencode.vM18.txt", sep="\t",col.names=NA)

table(is.na(Annotation_pc_Gencode$pc.ensembl_gene_ids))  # 28481 coding genes annotated de novo with 3 mismatches
rm(Annotation_pc_Gencode_2)
rm(Annotation_pc_Gencode_1)

################################################################################################
####  2.3. 2nd priority: Mapping to the whole transcriptome ####################################
################################################################################################

# note: to prevent memory exhaust first do everything for Annotation.1, then remove pc.inc.matrix.1 
# and start again with Annotation.2 / change entries in function etc. accordingly and remove pc.inc.matrix.2

##### 2.3.1 mapping: incidence matrix ##########################################################
################################################################################################
inc.matrix.1    <- vcountPDict(pdict.1, Gencode, max.mismatch = 3)  # cave: >10 Gb large!  
inc.matrix.2    <- vcountPDict(pdict.2, Gencode, max.mismatch = 3)    

##### 2.3.2 extract names of matched coding genes from fasta file ##############################
################################################################################################

#          1                    2                    3                  4                    5              6           7   8  
# ENSMUST00000193812.1|ENSMUSG00000102693.1|OTTMUSG00000049935.1|OTTMUST00000127109.1|RP23-271O17.1-001|RP23-271O17.1|1070|TEC|

getGeneID <- function(i) {ifelse(max(inc.matrix.2[i,])==0,NA,
                                 unlist(strsplit(names(Gencode[which(inc.matrix.2[i,]>0)[1]]), "|", fixed = TRUE))[2])}

ensembl_gene_ids     <- bplapply(seq_len(length(pdict.2)), getGeneID, BPPARAM = bpparam)
ensembl_gene_ids     <- unlist(ensembl_gene_ids, use.names=FALSE)    

Annotation_Gencode_1 <- data.frame(row.names = row.names(Annotation.1),ensembl_gene_ids = ensembl_gene_ids)
rm(inc.matrix.1) 
rm(ensembl_gene_ids) ### repeat with inc.matrix.2

Annotation_Gencode_2 <- data.frame(row.names = row.names(Annotation.2),ensembl_gene_ids = ensembl_gene_ids)
rm(inc.matrix.2)

Annotation_Gencode   <- rbind(Annotation_Gencode_1,Annotation_Gencode_2)
write.table(Annotation_Gencode, file = "Annotation_Gencode.vM18.txt", sep="\t",col.names=NA)

table(is.na(Annotation_Gencode$ensembl_gene_ids))  # 33361 genes annotated de novo using all Gencode.M18 transcripts
rm(Annotation_Gencode_1)
rm(Annotation_Gencode_2)
rm(Annotation.1)
rm(Annotation.2)


################################################################################################
#### 3. Final Annotation  ######################################################################
################################################################################################
Annotation_Gencode    <- read.delim("Annotation_Gencode.vM18.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1, na.strings = "NA")
Annotation_pc_Gencode <- read.delim("Annotation_pc_Gencode.vM18.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1, na.strings = "NA")

table(row.names(Annotation_Gencode)    == row.names(Annotation.ID026655))
table(row.names(Annotation_pc_Gencode) == row.names(Annotation.ID026655))

#### 3.1 1st priority: Gencode.M18 coding transcripts ##########################################  
################################################################################################
Annotation_FINAL              <- Annotation_pc_Gencode   # take the 28481 probes de novo mapped to coding genes
colnames(Annotation_FINAL)    <- c("GeneID_FINAL")
Annotation_FINAL$Source_FINAL <- ifelse(is.na(Annotation_FINAL$GeneID_FINAL) ,NA,"Gencode.M18_pc") 
table(is.na(Annotation_FINAL$GeneID_FINAL))

#### 3.2 2nd priority: Gencode.M18 all transcripts #############################################
################################################################################################
index.1 <- which(is.na(Annotation_FINAL$GeneID_FINAL) & !is.na(Annotation_Gencode$ensembl_gene_ids))  # 4880 probes are not in coding transcripts, but have annotation in all Gencode M18 transcripts
Annotation_FINAL$GeneID_FINAL[index.1]       <- Annotation_Gencode$ensembl_gene_ids[index.1]    
Annotation_FINAL$Source_FINAL[index.1]       <- "Gencode M18"
table(is.na(Annotation_FINAL$GeneID_FINAL))  # 33361 probes annotated by de novo mapping against Gencode 

### change ensembl_gene_id.version to ensembl_gene_id only ######################################
GENE.ID<- NULL                   
for (i in 1:39428) {GENE.ID[i] <- unlist(strsplit(Annotation_FINAL$GeneID_FINAL[i], ".", fixed = TRUE))[1] }        
Annotation_FINAL$GeneID_FINAL <- GENE.ID


#### 3.3. Annotate ensembl_gene_ids with BiomaRt ###############################################
################################################################################################
library(biomaRt)

ensembl=useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
BM <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","gene_biotype"),
               filters = 'ensembl_gene_id', values = Annotation_FINAL$GeneID_FINAL, mart = ensembl) # cave unordered output in the order the results are returned from the Biomart server

BM         <- BM[!duplicated(BM$ensembl_gene_id),]   # retain 24047 unique ensembl_gene_ids 
BM[BM==""] <- NA

### merge Annotations
Annotation_FINAL$ProbeID <- row.names(Annotation_FINAL)
df.merge                 <- merge(Annotation_FINAL,BM,by.x = "GeneID_FINAL",by.y = "ensembl_gene_id", all.x=TRUE)   
row.names(df.merge)      <- df.merge$ProbeID
df.merge                 <- df.merge[row.names(Annotation_FINAL),]

table(df.merge$GeneID_FINAL == Annotation_FINAL$GeneID_FINAL)
table(row.names(df.merge)   == row.names(Annotation_FINAL)) 

Annotation_FINAL <- df.merge
rm(df.merge)

table(is.na(Annotation_FINAL$GeneID_FINAL))
table(is.na(Annotation_FINAL$external_gene_name))
table(row.names(Annotation_FINAL) == row.names(Annotation.ID026655))

Annotation_FINAL <- Annotation_FINAL[,-3]
colnames(Annotation_FINAL) <- c("GeneID_FINAL","Source_FINAL", "GeneSymbol_FINAL", "GeneName_FINAL", "GeneType_FINAL")


#### 3.3 3rd priority: Agilent eArray Annotation when Gencode.M18 is NA   ######################
################################################################################################
index.2 <- which(is.na(Annotation_FINAL$GeneSymbol_FINAL) & !is.na(Annotation.ID026655$GeneSymbol))  # 2872 probes
Annotation_FINAL$GeneSymbol_FINAL[index.2] <- Annotation.ID026655$GeneSymbol[index.2]              
Annotation_FINAL$GeneName_FINAL[index.2]   <- Annotation.ID026655$GeneName[index.2]    
Annotation_FINAL$Source_FINAL[index.2]     <- "Agilent eArray"

table(is.na(Annotation_FINAL$GeneSymbol_FINAL))  # 36226 probes annotated until here

#### 3.4 split Annotation file for step 5 ######################################################
################################################################################################
Annotation_FINAL.eArray  <- Annotation_FINAL[index.2,]
Annotation_FINAL.Gencode <- Annotation_FINAL[-index.2,]

#### 3.4 get GeneType_FINAL and GeneID from GencodeM18 for probes with eArray Annotation #######
################################################################################################
library(rtracklayer)
library(dplyr)

# make Annotation Matrix from genomic Gencode.M18 gtf file with Symbols as row.names:
GenCodeM18  <- import("gencode.vM18.annotation.gtf")                        # gtf file downloaded from https://www.gencodegenes.org/mouse_releases/current.html
GenCodeM18  <- GenCodeM18[mcols(GenCodeM18)$type=="gene"]                   # subset for genes only :54146  
GenCodeM18  <- as.data.frame(GenCodeM18)
GenCodeM18  <- select(GenCodeM18,seqnames:gene_name,-score,-phase)
GenCodeM18  <- subset(GenCodeM18,!duplicated(GenCodeM18$gene_name))         # 54087 unique GeneSymbols
row.names(GenCodeM18) <- GenCodeM18$gene_name

GenCodeM18.sel <- GenCodeM18[Annotation_FINAL.eArray$GeneSymbol_FINAL,]     # subset for the 2872 probes with eArray Annotation
Annotation_FINAL.eArray$GeneType_FINAL <- GenCodeM18.sel$gene_type          
Annotation_FINAL.eArray$GeneID_FINAL   <- GenCodeM18.sel$gene_id

#### get rid of version number in Ensemble GeneIDs
GENE.ID<- NULL                   
for (i in 1:2872) {GENE.ID[i] <- unlist(strsplit(Annotation_FINAL.eArray$GeneID_FINAL[i], ".", fixed = TRUE))[1] }        
Annotation_FINAL.eArray$GeneID_FINAL <- GENE.ID

#### union of both sets again
Annotation_FINAL <- rbind(Annotation_FINAL.Gencode,Annotation_FINAL.eArray) 
Annotation_FINAL <- Annotation_FINAL[order(row.names(Annotation_FINAL)),]

table(row.names(Annotation_FINAL) == row.names(Annotation.ID026655))
table(is.na(Annotation_FINAL$GeneID_FINAL))
table(is.na(Annotation_FINAL$GeneSymbol_FINAL))

#### 3.6 Attach Agilent eArray metadata  #######################################################
################################################################################################
Annotation_FINAL$Sequence                          <- Annotation.ID026655$Sequence
Annotation_FINAL$Agilent_eArray_GeneSymbol         <- Annotation.ID026655$GeneSymbol
Annotation_FINAL$Agilent_eArray_GeneName           <- Annotation.ID026655$GeneName
Annotation_FINAL$Agilent_eArray_GenomicCoordinates <- Annotation.ID026655$GenomicCoordinates
Annotation_FINAL$Agilent_eArray_GenbankAccession   <- Annotation.ID026655$GenbankAccession
Annotation_FINAL$Agilent_eArray_EntrezGeneID       <- Annotation.ID026655$EntrezGeneID

#### 3.7 export ################################################################################
################################################################################################
write.table(Annotation_FINAL, file = "Annotation_SAGA_FINAL_20181128.txt", sep="\t",col.names=NA)

Annotation.known <- subset(Annotation_FINAL,!is.na(Annotation_FINAL$GeneSymbol_FINAL))
write.table(Annotation.known, file = "Annotation_SAGA_FINAL_KNOWN_20181128.txt", sep="\t",col.names=NA)

#### 3.10 characterize annotation ##############################################################
################################################################################################
GeneTypes         <- table(Annotation_FINAL$GeneType_FINAL)                   
write.table(GeneTypes, file = "GeneTypes_Annotation_SAGA_FINAL_39428.txt", sep="\t",col.names=NA)

SAGA.Top19.old    <- read.delim("SAGA_INBUILD_CORE19.txt",header=TRUE,sep="\t",stringsAsFactors =FALSE, row.names = 1)
SAGA.Top19.FINAL  <- Annotation_FINAL[row.names(SAGA.Top19.old),]
table(SAGA.Top19.old$GeneSymbol == SAGA.Top19.FINAL$GeneSymbol_FINAL)  

table(duplicated(Annotation_FINAL$GeneSymbol_FINAL))
table(duplicated(Annotation.ID026655$GeneSymbol))

mismatches.1 <- subset(Annotation_FINAL, Annotation_FINAL$GeneSymbol_FINAL != Annotation_FINAL$Agilent_eArray_GeneSymbol)  # 1638 differences between latest eArray and novel annotation
mismatches.2 <- subset(Annotation_FINAL, Annotation_FINAL$GeneSymbol_FINAL != Annotation.INBUILD$GeneSymbol)               # 6653 differences between old annotation and novel annotation
mismatches.3 <- Annotation.INBUILD[row.names(mismatches.2),]

## how many probes per gene 
A <- subset(Annotation_FINAL, !is.na(Annotation_FINAL$GeneSymbol_FINAL))

table(duplicated(A$GeneSymbol_FINAL))  # 25976 genes with only one probe
a <- table(A$GeneSymbol_FINAL)
max(a)
hist(a)
a.2 <- subset(a, a==2)
a.3 <- subset(a, a>3)
a.4 <- subset(a, a==4)
a.x <- subset(a, a>10)



