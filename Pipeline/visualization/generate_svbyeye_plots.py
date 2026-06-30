#!/usr/bin/env python3
import os
import sys
import glob
import pandas as pd
import argparse
import tempfile
import subprocess
import uuid

import multiprocessing
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(description="Extract original PAF coordinates and generate SVbyEye SVG plots.")
    parser.add_argument("--results-dir", required=True, help="Path to the 'results' directory of the pipeline.")
    parser.add_argument("--out-dir", required=True, help="Output directory for plots.")
    parser.add_argument("--r-script", required=True, help="Path to plot_svbyeye_single.R")
    parser.add_argument("--padding", type=int, default=0, help="Padding in bp around the original alignment.")
    parser.add_argument("--cores", type=int, default=1, help="Number of CPU cores to use for parallel plotting.")
    return parser.parse_args()

def extract_paf_slice(original_paf_path, target_contig, target_start, target_end, out_temp_paf, padding=0):
    """
    Extracts PAF lines for the target contig within the expanded (target_start - padding, target_end + padding) range.
    """
    matched = []
    with open(original_paf_path, 'r') as fin:
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            
            q_name = parts[0]
            q_start = int(parts[2])
            q_end = int(parts[3])
            aln_len = int(parts[10])
            
            # Check overlap and target sequence
            if target_contig in q_name and parts[5] == "chrM":
                if aln_len >= 200:
                    if (q_start <= target_end + padding) and (q_end >= target_start - padding):
                        matched.append(line)
                    
    with open(out_temp_paf, 'w') as fout:
        for line in matched:
            fout.write(line)
            
    return len(matched)

def process_single_locus(row_data, args, sample, haplotype):
    hg38_chr, hg38_start, original_contig, original_start, original_end, category, is_ref = row_data
    
    if pd.isna(original_contig) or original_contig is None:
        return
        
    locus_name = f"{hg38_chr}_{hg38_start}"
    folder_class = "is_ref" if is_ref == "YES" else "non_ref"
    
    # Determine target folder: out_dir / class / locus /
    locus_dir = os.path.join(args.out_dir, folder_class, locus_name)
    os.makedirs(locus_dir, exist_ok=True)
    
    pdf_filename = f"{hg38_chr}_{hg38_start}_{sample}_{haplotype}_Cat{category}.pdf"
    pdf_path = os.path.join(locus_dir, pdf_filename)
    
    if os.path.exists(pdf_path):
        print(f"  [SKIP] {pdf_filename} already exists.")
        return
        
    # Locate original PAF file
    paf_file = os.path.join(args.results_dir, "01_classification", f"{sample}.numt_classification.annotated.tsv")
    paf_file = paf_file.replace("01_classification", "data/paf").replace(".numt_classification.annotated.tsv", "_combined.paf")
    
    if not os.path.exists(paf_file):
        print(f"  [WARN] Missing PAF: {paf_file}")
        return
    
    # Extract relevant alignments
    tmp_paf_name = os.path.join(locus_dir, f"tmp_{uuid.uuid4().hex[:8]}.paf")
    written = extract_paf_slice(paf_file, original_contig, original_start, original_end, tmp_paf_name, padding=args.padding)
    
    if written > 0:
        title = f"{locus_name} | {sample} ({haplotype})"
        subtitle = f"Category: {category} | Orig Contig: {original_contig}:{original_start}-{original_end}"
        
        rscript_bin = "Rscript"
        if sys.executable:
            rscript_path = os.path.join(os.path.dirname(sys.executable), "Rscript")
            if os.path.exists(rscript_path):
                rscript_bin = rscript_path
        
        cmd = [
            rscript_bin, args.r_script,
            "-p", tmp_paf_name,
            "-o", pdf_path,
            "-t", title,
            "-s", subtitle
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  [ERROR] R script failed for {pdf_filename}:\n{result.stderr}")
        else:
            print(f"  [OK] Generated: {pdf_filename}")
    else:
        print(f"  [WARN] No alignment lines extracted for {original_contig}:{original_start}-{original_end}")
        
    # Cleanup temp
    if os.path.exists(tmp_paf_name):
        os.remove(tmp_paf_name)

def main():
    args = parse_args()
    
    # Locate liftover TSVs
    liftover_files = glob.glob(os.path.join(args.results_dir, "02_liftover", "*.classification_hg38.tsv"))
    if not liftover_files:
        print(f"[ERROR] No liftover TSVs found in {args.results_dir}/02_liftover")
        sys.exit(1)
        
    print(f"Found {len(liftover_files)} liftover files.")
    
    tasks = []
    
    # Process each file to collect tasks
    for liftover_file in liftover_files:
        filename = os.path.basename(liftover_file)
        # e.g., HG002.maternal.classification_hg38.tsv
        parts = filename.split('.')
        sample = parts[0] if len(parts) > 0 else "Unknown"
        haplotype = parts[1] if len(parts) > 1 else "Unknown"
        
        df = pd.read_csv(liftover_file, sep='\t')
        
        # Filter for validated NUMTs (if your DSA pipeline has this, otherwise just use all)
        if 'valid_blastn' in df.columns and 'valid_blastx' in df.columns:
            df = df[(df['valid_blastn'] == 'YES') & (df['valid_blastx'] == 'YES')]
            
        for idx, row in df.iterrows():
            # Get key info from DSA-NUMT-pipeline columns
            hg38_chr = row.get('hg38_chr', 'chrUn')
            hg38_start = row.get('hg38_start', 0)
            
            original_contig = row.get('contig', None)
            original_start = row.get('locus_start', 0)
            original_end = row.get('locus_end', 0)
            
            category = row.get('category', 'Unknown')
            is_ref = str(row.get('is_ref_numt', 'NO')).upper()
            
            row_data = (hg38_chr, hg38_start, original_contig, original_start, original_end, category, is_ref)
            tasks.append((row_data, args, sample, haplotype))
            
    print(f"Total plots to generate: {len(tasks)}")
    
    if args.cores > 1:
        print(f"Running in parallel with {args.cores} cores...")
        with multiprocessing.Pool(args.cores) as pool:
            pool.starmap(process_single_locus, tasks)
    else:
        for task in tasks:
            process_single_locus(*task)

if __name__ == "__main__":
    main()
