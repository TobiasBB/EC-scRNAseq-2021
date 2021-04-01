# EC-scRNAseq-2021
ScRNAseq analysis for the paper: **Decoding the molecular landscape of the developing spatial processing system and production of entorhinal stellate cell-like cells by a direct programming approach**.

## Code for analysis
The R code for the bioinformatics analysis is broken into 4 part, corresponding to the first 4 main figures:
1. Figure 1: Merging, clustering and analysis of entire dataset: Fig1_Analysis of parent dataset.R
2. Figure 2: Subclustering and analysis of inhibitory interneurons: Fig2_IN_subanalysis.R
3. Figure 3 Subclustering and analysis of oligodendrocytes: Fig3_Oligodendrocyte_subanalysis.R
4. Figure 4 and 5A: Subclustering and analysis of intermediate progenitors and excitatory neurons: Fig4-5C_IP-neuron_subanalysis.R

Required version of all dependent packages are stated as the library is loaded. See methods section of the paper for further desciption of the analysis.

## Availability of data
The datasets generated and used for the analysis are available at the NCBI repository with GEO accession number [GSE134482](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134482). The submitted data includes the raw sequencing data as fastq files together with the processed count matrix used in this study.

## Citation
If you use the code, or data from the publication, then please cite:

