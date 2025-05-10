# CoprofillingPipeline
# Spatial DNA Sequencing Data Processing Workflow

This repository provides a comprehensive and automated pipeline for processing spatial DNA sequencing data, covering preprocessing steps from raw FASTQ files through fragment counting per genomic bins. The workflow is optimized for performance and ease of use, ensuring accurate alignment, deduplication, and quantification of genomic fragments.

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

Before running the scripts, please set the paths to the external software in the `preprocess.py` script clearly:

```python
# Modify these paths to match your environment
BWA_PATH = "/path/to/your/bwa"
SAMTOOLS_PATH = "/path/to/your/samtools"
SAMBAMBA_PATH = "/path/to/your/sambamba"
```

---

## Pipeline Overview

The spatial DNA sequencing data processing pipeline comprises two main steps:

### Step 1: Preprocessing (`preprocess.py`)

This step performs the following operations:

* Alignment of FASTQ reads using BWA MEM.
* Conversion, sorting, and indexing of SAM/BAM files using Samtools.
* Deduplication using Sambamba.
* Splitting of deduplicated SAM files by barcode.
* Conversion of split SAM files into sorted BAM files.
* Filtering BAM files based on autosomal read count threshold.

#### Usage

```bash
python preprocess.py --fq1 <path_to_fastq1> --fq2 <path_to_fastq2> --ref <path_to_reference_genome> --output <output_directory> --threshold 5
```

* `--fq1`: Path to FASTQ file (read 1).
* `--fq2`: Path to FASTQ file (read 2, contains barcodes).
* `--ref`: Path to reference genome file.
* `--output`: Directory to save processed files.
* `--threshold`: Minimum number of autosomal reads required to retain BAM files.

The output will be organized into:

```
output_directory/
├── preProcess/
├── separatedSam/
├── separatedBam/
└── bamFile/  # Final filtered BAM files
```

---

### Step 2: Fragment Counting (`getFragmentsInfo.py`)

This step counts unique UMI-tagged DNA fragments within 5M base genomic windows.

#### Usage

```bash
python getFragmentsInfo.py --directory <bam_files_directory> --fq2path <path_to_fastq2> --output_tsv <output_tsv_path> --output_sum_csv <output_summary_csv_path> --genome human
```

* `--directory`: Path to directory containing sorted BAM files.
* `--fq2path`: Path to FASTQ file (read 2, contains UMIs).
* `--output_tsv`: Output path for detailed fragment counts (TSV format).
* `--output_sum_csv`: Output path for summary fragment counts per genomic bin (CSV format).
* `--genome`: Genome used ('human' or 'mouse').

---

## Output Explanation

* **Detailed TSV file** (`output_tsv`):

  * Rows represent spatial barcoded spots.
  * Columns represent genomic bins (5M base intervals).
  * Values are unique fragment counts per bin per spot.

* **Summary CSV file** (`output_sum_csv`):

  * Aggregates fragment counts across all spots.
  * Useful for quick quality control and data summaries.

---

## Performance Optimization

* The workflow uses parallelization (Python `multiprocessing`) for computational efficiency.
* Ensure your environment is configured to utilize multiple CPU cores effectively.

---

## Contact

For any issues or contributions, please create an issue or pull request on GitHub. We welcome community engagement and feedback to improve this pipeline.

