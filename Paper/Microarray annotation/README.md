# Microarray annotation

Since the annotation of microarray probes may change as the annotation of the genome advances, we re-annotated the 39,428 probes on the Agilent Whole Mouse Genome Oligo Microarray 4x44K v2 (Design ID 026655) by mapping the 60mer sequences to a recent release of the murine transcriptome (Gencode version M18, GRCm38.p6, release 07/2018). The transcript databases were downloaded as FASTA files for the 64,732 protein coding transcripts (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.pc_transcripts.fa.gz) and for all 136,535 coding and noncoding transcripts of the reference transcriptome (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.transcripts.fa.gz). The annotation was performed using R 3.5.1, Bioconductor 3.735 and the R package “Biostrings” 36.  The 60mers were first aligned to the protein coding transcripts of Gencode M18, second to all transcripts of Gencode M18 allowing a maximum of 3 mismatches. Using these parameters 33,361 out of 39428 probes could be successfully mapped to the Gencode M18 transcriptome. The mapping process retrieved an Ensembl-GeneID (e.g. “ENSMUSG00000020743”) for each probe with a hit in the Gencode transcriptome. The GeneID was further annotated using the “BiomaRt” R package to retrieve gene symbols, description and gene type from the Ensembl 94 database.  For probes that could not be annotated by Gencode, the annotation was taken from latest annotation file for the Whole Mouse Genome Oligo Microarray 4x44K v2 downloaded from Agilent eArray webservice (https://earray.chem.agilent.com/earray/, ID 026655, released October 2017) resulting in annotation of 2872 additional probes, resulting in 36,226 annotated probes in total.

## Availability of raw data and R workspace file

The following files can be downloaded [here](https://owncloud.gwdg.de/index.php/s/83axjKOt1JrviwS):
*	Agilent_Annotation_026655_D_AA_20171030_wo.txt
*	Microarray_annotation.RData
*	Annotation_Gencode.vM18.txt
*	Annotation_GencodeM18_3MM_BioMart_FINAL.R
*	Annotation_pc_Gencode.vM18.txt
*	Annotation_SAGA_FINAL_20181128.txt
*	Annotation_SAGA_FINAL_KNOWN_20181128.txt
*	gencode.vM18.annotation.gtf
*	gencode.vM18.pc_transcripts.fa
*	gencode.vM18.transcripts.fa
*	GeneTypes_Annotation_SAGA_FINAL_39428.txt
*	SAGA_INBUILD_Annotation.txt
