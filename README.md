# _In vivo_ Chemical genetics of Mycobacterium tuberculosis

Codes and commands used in the Oluoch et al In vivo CGI paper

# Requirements
## Bash code:
* subread (version 1.6.2)
* samtools (version 1.3)
## R code:
* R (Version 4.3.1)
* dplyr (version 1.1.2)
* ggplot2 (3.5.1)
* pheatmap (1.0.12)
* tidyr (version 1.3.0)
* tidyverse (version 2.0.0)
* DESeq2 (version 1.40.2)
* DEBRA (version 1.01)

## Structure:
### ./
Scripts for barcode count extractions from fastq files and DEBRA analysis for differential barcode representation.

### Data
* Counts - Raw barcode counts of individual hypomorphs across the four depletion levels and conditions from the main spleen-based screen.
* GSE - Folder with a modified KEGG gene function database; used with the gene-set enrichment analysis script.
* barcodes.fa - Unique barcode sequences assigned to each strain; used to extract strain counts from the raw fastq files.
