#!/usr/bin/env python3
import gzip
import os
import argparse
import csv
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial
import sys
import numpy as np


def get_candidate_barcodes(barcode_file):
    """
    Read the list of candidate barcodes from a file, one per line.
    """
    candidate_barcodes = []
    with open(barcode_file, 'r') as f:
        for line in f:
            barcode = line.strip()
            if barcode:  # Skip empty lines
                candidate_barcodes.append(barcode)
    return candidate_barcodes


def is_one_char_different(s1, s2):
    """
    Check if two strings of equal length differ by exactly one character.
    """
    if len(s1) != len(s2):
        return False
    diff_count = sum(1 for a, b in zip(s1, s2) if a != b)
    return diff_count == 1


def process_fastq_chunk(lines, candidate_barcodes):
    """
    Process a chunk of FASTQ lines, extracting and correcting barcodes.

    Args:
        lines (list[str]): FASTQ lines chunk (multiples of 4).
        candidate_barcodes (list[str]): Valid barcode sequences.

    Returns:
        dict[str, list[str]]: Mapping from combined cell barcode to list of read names.
    """
    cell_barcodes = defaultdict(list)
    for i in range(0, len(lines), 4):
        try:
            read_name = lines[i].strip()
            sequence = lines[i + 1].strip()

            # Extract barcode segments from sequence
            barcode1 = sequence[76:84]  # Barcode1 at positions 76-83
            barcode2 = sequence[38:46]  # Barcode2 at positions 38-45
            barcode3 = sequence[0:8]    # Barcode3 at positions 0-7

            corrected = []
            for bc in (barcode1, barcode2, barcode3):
                if bc in candidate_barcodes:
                    corrected.append(bc)
                else:
                    # Attempt single-character correction
                    match = next((cb for cb in candidate_barcodes if is_one_char_different(bc, cb)), None)
                    corrected.append(match)

            if all(corrected):  # All three barcodes valid or corrected
                cell_bc = ''.join(corrected)
                # Remove '@' and any trailing metadata
                clean_name = read_name.lstrip('@').split(' ')[0]
                cell_barcodes[cell_bc].append(clean_name)
        except IndexError:
            break
    return cell_barcodes


def extract_barcodes(fq2_file, candidate_barcodes, num_processes):
    """
    Extract cell barcodes from a FASTQ.gz file using multiprocessing.
    """
    cell_barcodes = defaultdict(list)
    lines = []
    with gzip.open(fq2_file, 'rt') as fq2:
        for line in fq2:
            lines.append(line)
            if len(lines) >= 400000:  # Process in ~100k read chunks (400k lines)
                chunks = [lines[i:i+4000] for i in range(0, len(lines), 4000)]
                with Pool(num_processes) as pool:
                    results = pool.starmap(process_fastq_chunk, [(chunk, candidate_barcodes) for chunk in chunks])
                for res in results:
                    for bc, reads in res.items():
                        cell_barcodes[bc].extend(reads)
                lines = []
        # Process any remaining lines
        if lines:
            chunks = [lines[i:i+4000] for i in range(0, len(lines), 4000)]
            with Pool(num_processes) as pool:
                results = pool.starmap(process_fastq_chunk, [(chunk, candidate_barcodes) for chunk in chunks])
            for res in results:
                for bc, reads in res.items():
                    cell_barcodes[bc].extend(reads)
    return cell_barcodes


def build_read_to_barcode_mapping(cell_barcodes):
    """
    Build a dict mapping read names to their cell barcode.
    """
    return {read: bc for bc, reads in cell_barcodes.items() for read in reads}


def filter_barcodes_by_elbow(cell_barcodes, threshold=100):
    """
    Filter out barcodes with fewer reads than threshold.
    """
    print(f"Barcode count before filtering: {len(cell_barcodes)}")
    print("Filtering barcodes with low read counts...")
    filtered = {bc: reads for bc, reads in cell_barcodes.items() if len(reads) >= threshold}
    print(f"Barcode count after filtering: {len(filtered)}")
    return filtered


def process_sam_chunk(chunk, header_lines, read_to_bc):
    """
    Process a chunk of SAM lines, grouping by barcode for matched reads.
    """
    local_map = defaultdict(list)
    for line in chunk:
        if line.startswith('@'):
            continue
        read_name = line.split('\t')[0]
        if read_name in read_to_bc:
            local_map[read_to_bc[read_name]].append(line)
    return local_map


def decompress_file(output_dir, file):
    """
    Decompress a single .sam.gz file to .sam and remove the compressed file.
    """
    if file.endswith('.sam.gz'):
        infile = os.path.join(output_dir, file)
        outfile = infile[:-3]
        with gzip.open(infile, 'rt') as fin, open(outfile, 'w') as fout:
            fout.writelines(fin)
        os.remove(infile)


def decompress_all_files(output_dir, num_processes):
    """
    Decompress all .sam.gz files in the output directory using multiprocessing.
    """
    files = [f for f in os.listdir(output_dir) if f.endswith('.sam.gz')]
    with Pool(num_processes) as pool:
        pool.map(partial(decompress_file, output_dir), files)


def split_sam_by_barcodes(sam_file, cell_barcodes, output_dir, num_processes):
    """
    Split a SAM file into per-barcode SAM.gz files using multiprocessing.
    """
    os.makedirs(output_dir, exist_ok=True)
    header_lines = []
    chunks = []
    chunk_size = 100000
    with open(sam_file, 'r') as sf:
        buf = []
        for line in sf:
            if line.startswith('@'):
                header_lines.append(line)
            else:
                buf.append(line)
                if len(buf) >= chunk_size:
                    chunks.append(buf)
                    buf = []
        if buf:
            chunks.append(buf)
    read_to_bc = build_read_to_barcode_mapping(cell_barcodes)
    with Pool(num_processes) as pool:
        results = pool.starmap(process_sam_chunk, [(c, header_lines, read_to_bc) for c in chunks])
    final_map = defaultdict(list)
    for local in results:
        for bc, reads in local.items():
            final_map[bc].extend(reads)
    for bc, reads in final_map.items():
        out_file = os.path.join(output_dir, f"{bc}.sam.gz")
        with gzip.open(out_file, 'at') as out:
            if os.stat(out_file).st_size == 0:
                out.writelines(header_lines)
            out.writelines(reads)
    print("Decompressing .sam.gz files...")
    decompress_all_files(output_dir, num_processes)


def main():
    parser = argparse.ArgumentParser(
        description="Split SAM file by cell barcodes extracted from FASTQ file."
    )
    parser.add_argument("fq2_file", help="Path to input FASTQ.gz file containing barcode information.")
    parser.add_argument("sam_file", help="Path to input SAM file to split by barcode.")
    parser.add_argument("barcode_file", help="Path to file listing candidate barcodes, one per line.")
    parser.add_argument("output_dir", help="Directory to store split SAM files.")
    args = parser.parse_args()

    fq2 = args.fq2_file
    sam = args.sam_file
    bc_file = args.barcode_file
    out_dir = args.output_dir

    num_proc = cpu_count()

    print("Reading candidate barcodes...")
    candidates = get_candidate_barcodes(bc_file)

    print("Extracting cell barcodes from FASTQ...")
    cell_bcs = extract_barcodes(fq2, candidates, num_proc)

    print("Filtering low-read barcodes...")
    cell_bcs = filter_barcodes_by_elbow(cell_bcs)

    print("Splitting SAM by cell barcodes...")
    split_sam_by_barcodes(sam, cell_bcs, out_dir, num_proc)

    print(f"Done! Split SAM files saved in directory: {out_dir}")


if __name__ == "__main__":
    main()
