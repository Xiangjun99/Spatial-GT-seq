#!/usr/bin/env python3
import os
import argparse

# Software paths for easy modification
BWA_PATH = "/home/sw2448/project/software/BWA/bwa/bwa"
SAMTOOLS_PATH = "/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools"
SAMBAMBA_PATH = "/home/sw2448/.conda/envs/theSAMBAMBA/bin/sambamba"


def run_command(command):
    """
    Execute a system command and print it.
    """
    print(f"Running command: {command}")
    os.system(command)


def bwa_mem(fq1, fq2, ref_genome, output_dir):
    """
    Align reads to the reference genome using BWA MEM, producing a SAM file.

    Args:
        fq1 (str): Path to the first FASTQ file (required).
        fq2 (str): Path to the second FASTQ file (ignored for single-end mode).
        ref_genome (str): Path to the reference genome FASTA file.
        output_dir (str): Directory to store output SAM file.

    Returns:
        str: Path to the generated SAM file.
    """
    sample_name = os.path.basename(fq1).split('_')[0]
    sam_path = os.path.join(output_dir, f"{sample_name}.sam")
    # Construct BWA MEM command (single-end mode)
    bwa_cmd = f"{BWA_PATH} mem -t 10 {ref_genome} {fq1} > {sam_path}"
    run_command(bwa_cmd)
    return sam_path


def preprocess_sam(sam_path, output_dir):
    """
    Convert SAM to BAM, sort, mark and remove duplicates, and convert back to SAM.

    Args:
        sam_path (str): Path to the input SAM file.
        output_dir (str): Base directory for preprocessing outputs.

    Returns:
        str: Path to the cleaned SAM file.
    """
    sample_name = os.path.basename(sam_path).split('.')[0]
    preprocess_dir = os.path.join(output_dir, "preProcess")
    os.makedirs(preprocess_dir, exist_ok=True)

    bam_path = os.path.join(preprocess_dir, f"{sample_name}.bam")
    sort_bam_path = os.path.join(preprocess_dir, f"sort_{sample_name}.bam")
    clean_bam_path = os.path.join(preprocess_dir, f"sort_{sample_name}_clean.bam")
    clean_sam_path = os.path.join(preprocess_dir, f"sort_{sample_name}_clean.sam")

    # 1. Convert SAM to BAM
    cmd1 = f"{SAMTOOLS_PATH} view -bS {sam_path} > {bam_path}"
    # 2. Sort the BAM file
    cmd2 = f"{SAMTOOLS_PATH} sort {bam_path} > {sort_bam_path}"
    # 3. Index the sorted BAM
    cmd3 = f"{SAMTOOLS_PATH} index {sort_bam_path}"
    # 4. Mark and remove duplicates using Sambamba
    cmd4 = f"{SAMBAMBA_PATH} markdup -r -p -t 2 {sort_bam_path} {clean_bam_path}"
    # 5. Convert cleaned BAM back to SAM format
    cmd5 = f"{SAMTOOLS_PATH} view -h -o {clean_sam_path} {clean_bam_path}"

    for cmd in [cmd1, cmd2, cmd3, cmd4, cmd5]:
        run_command(cmd)

    return clean_sam_path


def main():
    """
    Parse arguments and run BWA alignment and SAM preprocessing steps.
    """
    parser = argparse.ArgumentParser(description="BWA MEM Alignment and SAM File Preprocessing.")
    parser.add_argument("--fq1", required=True, help="Path to the first FASTQ file.")
    parser.add_argument("--fq2", required=False, help="Path to the second FASTQ file (optional for single-end mode).")
    parser.add_argument("--ref", required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("--output", required=True, help="Directory to store output files.")

    args = parser.parse_args()
    fq1 = args.fq1
    fq2 = args.fq2
    ref_genome = args.ref
    output_dir = args.output

    os.makedirs(output_dir, exist_ok=True)

    print("Starting BWA MEM alignment...")
    sam_path = bwa_mem(fq1, fq2, ref_genome, output_dir)

    print("Starting preprocessing of SAM file...")
    preprocess_sam(sam_path, output_dir)

    print("Processing complete.")


if __name__ == "__main__":
    main()
