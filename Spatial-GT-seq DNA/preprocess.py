import os
import shutil
import gzip
import re
import pysam
import argparse
from multiprocessing import Pool, cpu_count

# Define software paths for easy updates 
BWA_PATH = "/home/sw2448/project/software/BWA/bwa/bwa"
SAMTOOLS_PATH = "/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools"
SAMBAMBA_PATH = "/home/sw2448/.conda/envs/theSAMBAMBA/bin/sambamba"


def run_command(command):
    print(f"Running command: {command}")
    os.system(command)


def bwa_mem(fq1, fq2, ref_genome, output_dir):
    sample_name = os.path.basename(fq1).split('_')[0]
    sam_path = os.path.join(output_dir, f"{sample_name}.sam")
    bwa_cmd = f"{BWA_PATH} mem -t {cpu_count()} {ref_genome} {fq1} {fq2} > {sam_path}"
    run_command(bwa_cmd)
    return sam_path


def preprocess_sam(sam_path, output_dir):
    sample_name = os.path.basename(sam_path).split('.')[0]
    preprocess_dir = os.path.join(output_dir, "preProcess")
    os.makedirs(preprocess_dir, exist_ok=True)

    bam_path = os.path.join(preprocess_dir, f"{sample_name}.bam")
    sort_bam_path = os.path.join(preprocess_dir, f"sort_{sample_name}.bam")
    clean_bam_path = os.path.join(preprocess_dir, f"sort_{sample_name}_clean.bam")
    clean_sam_path = os.path.join(preprocess_dir, f"sort_{sample_name}_clean.sam")

    # Execute commands sequentially
    cmds = [
        f"{SAMTOOLS_PATH} view -bS {sam_path} -@ {cpu_count()} > {bam_path}",
        f"{SAMTOOLS_PATH} sort -@ {cpu_count()} {bam_path} > {sort_bam_path}",
        f"{SAMTOOLS_PATH} index {sort_bam_path}",
        f"{SAMBAMBA_PATH} markdup --remove-duplicates -t {cpu_count()} {sort_bam_path} {clean_bam_path}",
        f"{SAMTOOLS_PATH} view -h {clean_bam_path} > {clean_sam_path}"
    ]

    for cmd in cmds:
        run_command(cmd)

    return clean_sam_path


def write_to_file(task):
    barcode, lines, header_path, split_dir = task
    sam_path = os.path.join(split_dir, f"{barcode}.sam")
    if not os.path.exists(sam_path):
        shutil.copy(header_path, sam_path)
    with open(sam_path, 'a') as f:
        f.writelines(lines)


def split_sam_optimized(clean_sam_path, fq2, output_dir):
    split_dir = os.path.join(output_dir, "separatedSam")
    os.makedirs(split_dir, exist_ok=True)

    readsID_barcode_dic = {}
    line_i = 0

    # Read FASTQ file for barcode extraction
    if fq2.endswith('.gz'):
        fq2_file = gzip.open(fq2, 'rt')
    else:
        fq2_file = open(fq2, 'r')

    oneBarcodeArr = [
        'AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT',
        'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA',
        'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA',
        'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA',
        'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA',
        'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA',
        'GATAGACA', 'GCCACATA'
    ]

    def if_barcode_match(a, b):
        # Require at least 7 matching bases
        match_base_amount = 7
        return sum(1 for i in range(len(a)) if a[i] == b[i]) >= match_base_amount

    # Build mapping from read IDs to barcodes
    while True:
        line = fq2_file.readline()
        if not line:
            break
        line_i += 1
        if line_i % 4 == 1:
            match = re.search(r'@(.*?)\s', line)
            thisKEY = match.group(1)
        elif line_i % 4 == 2:
            BARCODE_B = line[32:40]
            BARCODE_A = line[70:78]
            match_flagA = match_flagB = 0
            for barcode in oneBarcodeArr:
                if if_barcode_match(barcode, BARCODE_B):
                    match_flagB = 1
                    this_BARCODE_B = barcode
                if if_barcode_match(barcode, BARCODE_A):
                    match_flagA = 1
                    this_BARCODE_A = barcode
                if match_flagB and match_flagA:
                    break
            if match_flagB and match_flagA:
                readsID_barcode_dic[thisKEY] = f"{this_BARCODE_B}+{this_BARCODE_A}"
    fq2_file.close()

    print("ReadsID to barcode dictionary created.")

    header_path = os.path.join(split_dir, "header.sam")
    run_command(f"{SAMTOOLS_PATH} view -H {clean_sam_path} > {header_path}")

    grouped_lines = {}
    with open(clean_sam_path, 'r') as fq1sam:
        for line in fq1sam:
            if line.startswith("@"):
                continue
            read_id = line.split("\t")[0]
            barcode = readsID_barcode_dic.get(read_id)
            if barcode:
                grouped_lines.setdefault(barcode, []).append(line)

    tasks = [(barcode, lines, header_path, split_dir) for barcode, lines in grouped_lines.items()]
    with Pool(cpu_count()) as pool:
        pool.map(write_to_file, tasks)

    print("SAM file splitting completed.")
    return split_dir


def process_task(task):
    """
    Process a single barcode's SAM to BAM conversion task.
    """
    barcode_name, sam_path, bam_dir = task
    bam_path = os.path.join(bam_dir, f"{barcode_name}.bam")
    sorted_bam_path = os.path.join(bam_dir, f"sort_{barcode_name}.bam")

    # Execute view, sort, index sequentially
    cmds = [
        f"{SAMTOOLS_PATH} view -bS {sam_path} -@ {cpu_count()} > {bam_path}",
        f"{SAMTOOLS_PATH} sort -@ {cpu_count()} {bam_path} > {sorted_bam_path}",
        f"{SAMTOOLS_PATH} index {sorted_bam_path}"
    ]

    for cmd in cmds:
        run_command(cmd)


def convert_sam_to_bam(split_dir):
    bam_dir = os.path.join(os.path.dirname(split_dir), "separatedBam")
    os.makedirs(bam_dir, exist_ok=True)

    oneBarcodeArr = [
        'AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT',
        'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA',
        'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA',
        'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA',
        'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA',
        'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA',
        'GATAGACA', 'GCCACATA'
    ]

    tasks = []
    for BarcodeA in oneBarcodeArr:
        for BarcodeB in oneBarcodeArr:
            thisBarcodeName = f"{BarcodeA}+{BarcodeB}"
            thisSamPath = os.path.join(split_dir, f"{thisBarcodeName}.sam")
            if os.path.exists(thisSamPath):
                tasks.append((thisBarcodeName, thisSamPath, bam_dir))

    with Pool(cpu_count()) as pool:
        pool.map(process_task, tasks)

    return bam_dir


def filter_and_copy_bam(bam_dir, final_output_dir, threshold):
    # Copy BAM files that meet autosomal read threshold
    os.makedirs(final_output_dir, exist_ok=True)

    def calculate_autosomal_reads(bam_path):
        total_autosomal_reads = 0
        autosomes = {f"chr{i}" for i in range(1, 23)}
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for read in bam_file:
                if read.reference_name in autosomes:
                    total_autosomal_reads += 1
        return total_autosomal_reads

    def filter_task(bam_path):
        if calculate_autosomal_reads(bam_path) >= threshold:
            shutil.copy(bam_path, final_output_dir)
            bai_file = f"{bam_path}.bai"
            if os.path.exists(bai_file):
                shutil.copy(bai_file, final_output_dir)

    for root, _, files in os.walk(bam_dir):
        for bam_file in files:
            if bam_file.startswith("sort_") and bam_file.endswith(".bam"):
                filter_task(os.path.join(root, bam_file))


def main():
    parser = argparse.ArgumentParser(description="Automate sequencing analysis workflow up to BAM generation.")
    parser.add_argument("--fq1", required=True, help="Path to the first FASTQ file.")
    parser.add_argument("--fq2", required=True, help="Path to the second FASTQ file.")
    parser.add_argument("--ref", required=True, help="Path to the reference genome file.")
    parser.add_argument("--output", required=True, help="Path to the output directory.")
    parser.add_argument("--threshold", type=int, default=5, help="Threshold for autosomal reads.")

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)

    print("Starting BWA MEM alignment...")
    sam_path = bwa_mem(args.fq1, args.fq2, args.ref, args.output)

    print("Starting SAM file preprocessing...")
    clean_sam_path = preprocess_sam(sam_path, args.output)

    print("Splitting SAM file by barcode...")
    split_dir = split_sam_optimized(clean_sam_path, args.fq2, args.output)

    print("Converting split SAM to BAM files...")
    bam_dir = convert_sam_to_bam(split_dir)

    final_output_dir = os.path.join(args.output, "bamFile")
    print("Filtering and copying BAM files...")
    filter_and_copy_bam(bam_dir, final_output_dir, args.threshold)

    print("Analysis up to BAM generation complete.")


if __name__ == "__main__":
    main()
