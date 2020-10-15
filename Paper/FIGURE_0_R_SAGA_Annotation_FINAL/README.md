# FIGURE 0 - R SAGA Annotation FINAL

All 39,428 probes on the Agilent Whole Mouse Genome Oligo Microarray 4x44K v2 (Design ID 026655) were re-annotated by mapping the 60mer sequences to Gencode version M18 murine reference transcriptome (GRCm38.p6, release 07/2018). The alignment of the 60mer probes to the reference transcriptome was performed using the R package “Biostrings”. Protein coding sequences were prioritized by first aligning all 60mers to protein coding transcripts of Gencode M18, and second to all transcripts of Gencode M18 allowing a maximum of 3 mismatches per 60mer. Using these parameters 33,361 out of 39,428 probes were successfully mapped to the Gencode M18 transcriptome. The mapping process retrieved an Ensembl-GeneID (e.g. “ENSMUSG00000020743”) for each probe with a hit in the Gencode transcriptome. The Ensembl-GeneID was further annotated using the “BiomaRt” R package to retrieve gene symbols, description and gene type from the Ensembl 94 database. For probes that could not be annotated by Gencode, annotation was taken from the latest annotation file for the Whole Mouse Genome Oligo Microarray 4x44K v2 downloaded from Agilent eArray web service (https://earray.chem.agilent.com/earray/, ID 026655, released October 2017) resulting in annotation of 2,872 additional probes and 36,226 annotated probes in total.

## Availability of R script and workspace file

*	[FIGURE_0_R_SAGA_Annotation_FINAL.RData](https://www.dropbox.com/s/93mupp9o2x7wvjj/.RData?dl=0)
*	[Annotation_GencodeM18_3MM_BioMart_FINAL.R](./Annotation_GencodeM18_3MM_BioMart_FINAL.R)
