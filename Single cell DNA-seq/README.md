# Single-Cell DNA Sequencing Data Processing

This section handles high-resolution DNA sequencing data from individual cells. The barcode-aware preprocessing allows reconstruction of single-cell data from pooled reads.

## Step 1: Alignment and Preprocessing (`preprocess1.py`)

Aligns single-end reads, performs BAM conversion and deduplication, and outputs a cleaned SAM file.

### Usage

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

## Step 2: Barcode Extraction and SAM Splitting (`preprocess2.py`)

Extracts three barcodes per read, performs error correction (1 mismatch tolerated), filters by read count, and splits SAM file into per-cell files.

### Required Input

* **Barcode whitelist file**: Found in this repository under:

```
testData/barcode/barcode_list.txt
```

* **FASTQ file**: R2 read in `.fq.gz` format with barcodes encoded at positions 0–8, 38–46, and 76–84

### Usage

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

* **Per-cell SAM files (single-cell)**: Enable independent cell-level variant calling or fragment reconstruction.

---

## Performance Optimization

* Scripts utilize Python's `multiprocessing.Pool` to maximize CPU utilization.
* Ensure adequate system memory and thread capacity when processing large datasets.

---
