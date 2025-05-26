# Spatial-GT-seq

This repository includes the codes for the manuscript:'Spatial genome and transcriptome co-mapping reveals clonal heterogeneity in tumor tissue'.

It includes: 1) codes for processing and analyzing spatial-DNA-seq data, 2) codes for processing and analyzing spatial-RNA-seq data, and 3) codes for processing and analyzing single-cell DNA-seq data. Specifically, this repository is organized as following folders:

| Folder        | Purpose                                                                                                                                     | Visibility |
|---------------|---------------------------------------------------------------------------------------------------------------------------------------------|------------|
| `install`             | Installing the required command line tools and packages of R and Python languages.                                                          | public     |
| `Spatial-DNA-seq`     | This tutorial uses a real sample in the manuscript to demonstrate how to preprocess and analyze wellDA-seq data.                            | public     |
| `Spatial-RNA-seq`     | Scripts to take the BCL files of wellDA-seq to create the single-cell CNA and ATAC matrices files.                                        | public     |
| `Single-cell DNA-seq` | Codes of computational analysis used in the manuscript. This folder is organized in the order of applications.                             | private    |
| `testData`            | Codes of data visualization shown in the manuscript. This folder is organized in the order of figures.                                     | private    |

---
# Spatial and Single-Cell DNA Sequencing Data Processing Workflow

## Requirements

### Programming Language

* **Python 3.11** is required for compatibility with all included scripts.

### Python Packages

Install the required Python dependencies using `pip`:

```bash
pip install pysam natsort numpy
```

### External Software Tools

The following tools are required and should be installed separately. These tools are used via absolute paths defined in the scripts.

```bash
conda install -c bioconda bwa samtools sambamba
```

You must modify the absolute paths in the scripts (`preprocess.py`, `preprocess1.py`) to reflect the location of these tools on your system:

```python
BWA_PATH = "/path/to/bwa"
SAMTOOLS_PATH = "/path/to/samtools"
SAMBAMBA_PATH = "/path/to/sambamba"
```

---

## Spatial DNA Sequencing Data Processing

This part of the pipeline processes spatially barcoded bulk DNA sequencing data and outputs barcode-specific BAM files and fragment count matrices.

### Step 1: Preprocessing (`preprocess.py`)

This script performs:

* Read alignment using BWA MEM
* SAM to BAM conversion, sorting, and indexing using Samtools
* Duplicate marking and removal using Sambamba
* Splitting deduplicated SAM files by barcode (extracted from FASTQ read 2)
* Per-barcode SAM to BAM conversion and autosomal read count filtering

#### Usage

```bash
python preprocess.py \
  --fq1 path/to/sample_R1.fq.gz \
  --fq2 path/to/sample_R2.fq.gz \
  --ref path/to/reference.fa \
  --output path/to/output_directory \
  --threshold 5
```

The final filtered BAM files will be located in:

```
<output_directory>/bamFile/
```

### Step 2: Fragment Counting (`getFragmentsInfo.py`)

This script counts unique UMI-tagged fragments for each BAM file in 5 million base pair bins across all chromosomes.

#### Usage

```bash
python getFragmentsInfo.py \
  --directory path/to/bamFile \
  --fq2path path/to/sample_R2.fq.gz \
  --output_tsv fragments_matrix.tsv \
  --output_sum_csv genome_bin_summary.csv \
  --genome mouse
```

Output:

* `fragments_matrix.tsv`: Each row represents a barcode spot; each column is a genomic bin.
* `genome_bin_summary.csv`: Total fragment count per genomic bin.

---

## Single-Cell DNA Sequencing Data Processing

This section handles high-resolution DNA sequencing data from individual cells. The barcode-aware preprocessing allows reconstruction of single-cell data from pooled reads.

### Step 1: Alignment and Preprocessing (`preprocess1.py`)

Aligns single-end reads, performs BAM conversion and deduplication, and outputs a cleaned SAM file.

#### Usage

```bash
python preprocess1.py \
  --fq1 path/to/sample.fq.gz \
  --ref path/to/reference.fa \
  --output path/to/output_directory
```

Output:

```
<output_directory>/preProcess/sort_<sample>_clean.sam
```

### Step 2: Barcode Extraction and SAM Splitting (`preprocess2.py`)

Extracts three barcodes per read, performs error correction (1 mismatch tolerated), filters by read count, and splits SAM file into per-cell files.

#### Required Input

* **Barcode whitelist file**: Found in this repository under:

```
testData/barcode/barcode_list.txt
```

* **FASTQ file**: R2 read in `.fq.gz` format with barcodes encoded at positions 0–8, 38–46, and 76–84

#### Usage

```bash
python preprocess2.py \
  path/to/sample_R2.fq.gz \
  path/to/sort_<sample>_clean.sam \
  testData/barcode/barcode_list.txt \
  path/to/output_split_sam
```

Each resulting SAM file corresponds to one cell barcode and is named `<barcode>.sam` in the output directory.

---

## Output Explanation

* **BAM files (spatial)**: Ready for downstream CNV calling or genome-wide coverage analysis.
* **Per-cell SAM files (single-cell)**: Enable independent cell-level variant calling or fragment reconstruction.
* **Fragment matrix**: Useful for clustering, dimensionality reduction, and signal heatmaps.

---

## Performance Optimization

* Scripts utilize Python's `multiprocessing.Pool` to maximize CPU utilization.
* Ensure adequate system memory and thread capacity when processing large datasets.

---

## Contact

For questions or contributions, please file an issue or submit a pull request. We encourage collaborative development and constructive feedback to further improve this project.
