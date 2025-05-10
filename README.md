# CoprofillingPipeline
# Spatial and Single-Cell DNA Sequencing Data Processing Workflow

This repository provides comprehensive and automated pipelines for processing spatial and single-cell DNA sequencing data. The workflows are optimized for performance and ease of use, ensuring accurate alignment, deduplication, and quantification of genomic fragments.

---

## Requirements

### Programming Languages

* **Python 3.11**

### Python Dependencies

The following Python libraries are required:

```bash
pip install pysam natsort
```

### External Software Dependencies

You must install and configure the following external bioinformatics tools:

```bash
conda install -c bioconda bwa samtools sambamba
```

---

## Setting Up the Environment

Modify the paths to external software clearly in each script (`preprocess.py` and `preprocess1.py`):

```python
# Modify these paths to match your environment
BWA_PATH = "/path/to/your/bwa"
SAMTOOLS_PATH = "/path/to/your/samtools"
SAMBAMBA_PATH = "/path/to/your/sambamba"
```

---

## Spatial DNA Data Processing

### Step 1: Preprocessing (`preprocess.py`)

Performs alignment, deduplication, and BAM file generation.

#### Usage

```bash
python preprocess.py --fq1 <path_to_fastq1> --fq2 <path_to_fastq2> --ref <path_to_reference_genome> --output <output_directory> --threshold 5
```

### Step 2: Fragment Counting (`getFragmentsInfo.py`)

Counts unique UMI-tagged DNA fragments within 5M base genomic windows.

#### Usage

```bash
python getFragmentsInfo.py --directory <bam_files_directory> --fq2path <path_to_fastq2> --output_tsv <output_tsv_path> --output_sum_csv <output_summary_csv_path> --genome human
```

---

## Single-Cell DNA Data Processing

The single-cell DNA data processing pipeline consists of two main preprocessing steps:

### Step 1: Alignment and Initial Preprocessing (`preprocess1.py`)

Align reads using BWA MEM, convert SAM to sorted, deduplicated BAM, and back to SAM format.

#### Usage

```bash
python preprocess1.py --fq1 <path_to_fastq1> --ref <path_to_reference_genome> --output <output_directory>
```

### Step 2: Barcode Extraction and SAM Splitting (`preprocess2.py`)

Extracts barcodes from FASTQ files, filters barcodes based on read count, and splits SAM files according to barcodes.

#### Required Barcode File

Barcode file can be found in the repository at:

```
testData/barcode/barcode_list.txt
```

#### Usage

```bash
python preprocess2.py <path_to_fastq2.gz> <path_to_cleaned_sam_file> <path_to_barcode_file> <output_directory>
```

* `<path_to_fastq2.gz>`: Path to FASTQ file containing barcode information.
* `<path_to_cleaned_sam_file>`: Path to the cleaned SAM file from step 1.
* `<path_to_barcode_file>`: Path to the candidate barcode file.
* `<output_directory>`: Directory to store split SAM files.

---

## Output Explanation

* **Detailed TSV/CSV files**:

  * Provides comprehensive counts and metrics suitable for downstream analyses.
* **Split SAM files**:

  * Individual SAM files per barcode, facilitating single-cell resolution analysis.

---

## Performance Optimization

* All workflows utilize parallelization via Python's `multiprocessing` module.
* Ensure your computational environment effectively utilizes available CPU cores.

---

## Contact

For issues or contributions, please create an issue or pull request on GitHub. We welcome community engagement and feedback to improve these pipelines.
