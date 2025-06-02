# Spatial-GT-seq

This repository contains the code for the manuscript "Spatial genome and transcriptome co-mapping reveals clonal heterogeneity in tumor tissue." In this work, we introduce Spatial-GT-seq, a spatial multi-omic technology that enables genome-wide co-profiling of the genome and transcriptome at near single-cell resolution. Simultaneous mapping of the genome and transcriptome within the same tumor tissue allows for direct linkage of clonal architecture to transcriptional states, offering deeper insights into both intrinsic and extrinsic sources of tumor heterogeneity within the spatial context of the tissue.



![illustration for github](https://github.com/user-attachments/assets/0b5e3163-86e3-4958-bbe9-35fbfc6ed779)



It includes: 1) codes for processing and analyzing spatial-DNA-seq data, 2) codes for processing and analyzing spatial-RNA-seq data, and 3) codes for processing and analyzing single-cell DNA-seq data. Specifically, this repository is organized as following folders:


| Folder                | Purpose                                                                                                                                             | Visibility |
|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| `install`             | Installing the required command line tools and packages of R and Python languages.                                                          | public     |
| `Spatial-GT-seq DNA`     | Scripts for alignment, deduplication, barcode parsing, and fragment quantification to generate BAM files compatible with CopyKit. Additionally, this folder contains code for downstream analyses used in Spatial-GT-seq DNA.        | public     |
| `Spatial-GT-seq RNA`     | Scripts to preprocess Spatial-GT-seq RNA data, as well as code for downstream analysis and visualization.                                  | public     |
| `Single-cell DNA-seq` | Scripts to preprocess single-cell DNA-seq data, as well as code for downstream analysis and visualization.                                  | public   |
| `testData`            | A sample dataset for testing the spatial-DNA-seq pipeline.                                                                                  | public    |


---

