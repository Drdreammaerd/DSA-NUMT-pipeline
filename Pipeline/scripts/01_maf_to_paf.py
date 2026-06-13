#!/usr/bin/env python3
###################################
# Purpose: Convert LAST MAF → PAF format WITHOUT requiring LAST/Docker.
#          This is a standalone Python fallback if maf-convert is unavailable.
#          Also extracts specific loci for SVbyEye visualization.
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# AI Assistant: Gemini
# Usage: python3 02b_maf_to_paf_python.py
###################################

import os
import sys

# =========== CONFIG ===========
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
# ==============================


def maf_to_paf(maf_path, paf_path):
    """
    Convert LAST MAF to PAF format.
    
    PAF columns:
    1. query_name  2. query_len  3. query_start  4. query_end  5. strand
    6. target_name 7. target_len 8. target_start  9. target_end
    10. num_matches 11. alignment_block_len 12. mapping_quality
    """
    count = 0
    with open(maf_path, 'r') as fin, open(paf_path, 'w') as fout:
        score = 0
        s_lines = []
        
        for line in fin:
            line = line.rstrip('\n')
            
            if line.startswith('#') or line == '':
                if len(s_lines) >= 2:
                    paf_line = _convert_block(s_lines, score)
                    if paf_line:
                        fout.write(paf_line + '\n')
                        count += 1
                s_lines = []
                continue
            
            if line.startswith('a '):
                if len(s_lines) >= 2:
                    paf_line = _convert_block(s_lines, score)
                    if paf_line:
                        fout.write(paf_line + '\n')
                        count += 1
                s_lines = []
                # Parse score
                score = 0
                for part in line.split():
                    if part.startswith('score='):
                        try:
                            score = int(float(part.split('=')[1]))
                        except:
                            score = 0
                continue
            
            if line.startswith('s '):
                s_lines.append(line)
        
        # Last block
        if len(s_lines) >= 2:
            paf_line = _convert_block(s_lines, score)
            if paf_line:
                fout.write(paf_line + '\n')
                count += 1
    
    return count


def _convert_block(s_lines, score):
    """Convert a pair of MAF 's' lines to PAF format."""
    # In LAST MAF: first s-line = reference (chrM), second = query (contig)
    ref_parts = s_lines[0].split()
    qry_parts = s_lines[1].split()
    
    if len(ref_parts) < 7 or len(qry_parts) < 7:
        return None
    
    # Parse reference (target in PAF terms)
    t_name = ref_parts[1]
    t_start_raw = int(ref_parts[2])
    t_aln_size = int(ref_parts[3])
    t_strand = ref_parts[4]
    t_srcsize = int(ref_parts[5])
    t_seq = ref_parts[6]
    
    # Parse query
    q_name = qry_parts[1]
    q_start_raw = int(qry_parts[2])
    q_aln_size = int(qry_parts[3])
    q_strand = qry_parts[4]
    q_srcsize = int(qry_parts[5])
    q_seq = qry_parts[6]
    
    # Determine relative strand
    if t_strand == q_strand:
        strand = '+'
    else:
        strand = '-'
    
    # Convert coordinates to forward strand
    if t_strand == '+':
        t_start = t_start_raw
        t_end = t_start_raw + t_aln_size
    else:
        t_start = t_srcsize - t_start_raw - t_aln_size
        t_end = t_srcsize - t_start_raw
    
    if q_strand == '+':
        q_start = q_start_raw
        q_end = q_start_raw + q_aln_size
    else:
        q_start = q_srcsize - q_start_raw - q_aln_size
        q_end = q_srcsize - q_start_raw
    
    # Count matches and block length
    matches = 0
    block_len = 0
    for tc, qc in zip(t_seq, q_seq):
        if tc != '-' or qc != '-':
            block_len += 1
        if tc != '-' and qc != '-' and tc.upper() == qc.upper():
            matches += 1
    
    # Mapping quality (use score as proxy, cap at 255)
    mapq = min(255, max(0, score // 10))
    
    # PAF format (12 mandatory columns + optional)
    paf = (
        f"{q_name}\t{q_srcsize}\t{q_start}\t{q_end}\t{strand}\t"
        f"{t_name}\t{t_srcsize}\t{t_start}\t{t_end}\t"
        f"{matches}\t{block_len}\t{mapq}\t"
        f"AS:i:{score}"
    )
    
    return paf


def extract_cases(paf_path, cases, case_dir):
    """Extract specific loci from PAF into individual files."""
    os.makedirs(case_dir, exist_ok=True)
    
    # Read all PAF lines
    with open(paf_path, 'r') as f:
        lines = f.readlines()
    
    for sample, contig, start, end, label in cases:
        case_paf = os.path.join(case_dir, f"{label}.paf")
        matched = []
        
        for line in lines:
            parts = line.split('\t')
            if len(parts) < 12:
                continue
            q_name = parts[0]
            q_start = int(parts[2])
            q_end = int(parts[3])
            
            if q_name == contig and q_start < end and q_end > start:
                matched.append(line)
        
        with open(case_paf, 'w') as f:
            f.writelines(matched)
        
        print(f"  [OK] {label}: {len(matched)} alignment blocks")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Convert LAST MAF to PAF format.')
    parser.add_argument('--maf', help='Input MAF file path')
    parser.add_argument('--out', help='Output PAF file path')
    args = parser.parse_args()
    
    # CLI mode: convert single MAF → PAF
    if args.maf and args.out:
        if not os.path.exists(args.maf):
            print(f"[ERROR] MAF not found: {args.maf}", file=sys.stderr)
            sys.exit(1)
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        count = maf_to_paf(args.maf, args.out)
        print(f"{count} alignment records → {args.out}")
        return
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

