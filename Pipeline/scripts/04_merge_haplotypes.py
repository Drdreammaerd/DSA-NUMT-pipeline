#!/usr/bin/env python3
###################################
# Purpose: Merge maternal + paternal haplotypes per individual.
#          Determines zygosity (HOM/HET) and structural concordance.
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# AI Assistant: Gemini
# Usage: python3 scripts/05_merge_haplotypes.py
###################################

import os
import sys
import re
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

MATCH_DIST = 2000  # bp tolerance for position matching
MT_COORD_TOLERANCE = 100  # bp tolerance for mt_structure segment matching
COMPLEX_CATS = {"C_inversion", "D_tandem_repeat", "E_complex_chimeric"}
CAT_RANK = {"A_single_block": 0, "B_divergence_gap": 1, "C_inversion": 2,
            "D_tandem_repeat": 3, "E_complex_chimeric": 4}


# ── Parsing ──

def parse_tsv(path):
    """Parse liftover classification TSV."""
    loci = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            while len(fields) < len(header):
                fields.append('')
            row = dict(zip(header, fields))
            row['chr'] = row['hg38_chr']
            row['start'] = int(row['hg38_start'])
            row['end'] = int(row['hg38_end'])
            loci.append(row)
    return loci


def parse_mt_structure(mt_str):
    """Parse mt_structure string into list of (start, end, strand) tuples.
    e.g. '992-1514(-);5603-6982(+)' → [(992, 1514, '-'), (5603, 6982, '+')]
    """
    if not mt_str or mt_str == '.':
        return []
    segments = []
    for part in mt_str.split(';'):
        m = re.match(r'(\d+)-(\d+)\(([+-])\)', part.strip())
        if m:
            segments.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    return segments


# ── Structure comparison ──

def _flip_strand(strand):
    return '-' if strand == '+' else '+'


def _segments_match(segs_a, segs_b, tol):
    """Check if two segment lists match (same order). Returns (coords_match, strand_match)."""
    if len(segs_a) != len(segs_b):
        return False, False
    coords_match = True
    strand_match = True
    for (s_a, e_a, st_a), (s_b, e_b, st_b) in zip(segs_a, segs_b):
        if abs(s_a - s_b) > tol or abs(e_a - e_b) > tol:
            coords_match = False
            break
        if st_a != st_b:
            strand_match = False
    return coords_match, strand_match


def compare_mt_structure(struct_a, struct_b, tol=MT_COORD_TOLERANCE):
    """Compare two mt_structure strings. Returns similarity level.
    
    Accounts for contig orientation: if one assembly reverses the contig,
    segments appear in reverse order with flipped strands. This is NOT
    a structural difference.
    
    Returns:
        'IDENTICAL'    - exact same segments
        'SIMILAR'      - same segments, coords within tolerance
        'STRAND_FLIP'  - same coords, strand(s) differ (forward order)
        'REV_COMP'     - same segments but reversed order + strand flip (contig orientation)
        'DIFFERENT'    - truly different structure
    """
    segs_a = parse_mt_structure(struct_a)
    segs_b = parse_mt_structure(struct_b)
    
    if not segs_a or not segs_b:
        return 'UNKNOWN'
    
    # Exact match
    if segs_a == segs_b:
        return 'IDENTICAL'
    
    # Different number of segments → truly different
    if len(segs_a) != len(segs_b):
        return 'DIFFERENT'
    
    # Forward comparison
    fwd_coords, fwd_strand = _segments_match(segs_a, segs_b, tol)
    
    if fwd_coords and fwd_strand:
        return 'SIMILAR'
    elif fwd_coords and not fwd_strand:
        return 'STRAND_FLIP'
    
    # Reverse-complement comparison:
    # Reverse order of segs_b and flip strands
    segs_b_rc = [(s, e, _flip_strand(st)) for s, e, st in reversed(segs_b)]
    rev_coords, rev_strand = _segments_match(segs_a, segs_b_rc, tol)
    
    if rev_coords and rev_strand:
        return 'REV_COMP'
    elif rev_coords and not rev_strand:
        return 'REV_COMP'  # coords match in rev order (strand may partially match)
    
    return 'DIFFERENT'


def consensus_category(cat_a, cat_b):
    """Return the more complex category."""
    rank_a = CAT_RANK.get(cat_a, -1)
    rank_b = CAT_RANK.get(cat_b, -1)
    return cat_a if rank_a >= rank_b else cat_b


# ── Merging ──

def merge_haplotypes(mat_loci, pat_loci, dist=MATCH_DIST):
    """Merge maternal and paternal loci. Returns list of merged records."""
    
    # Index paternal by chr
    pat_by_chr = defaultdict(list)
    for p in pat_loci:
        pat_by_chr[p['chr']].append(p)
    
    pat_matched = set()  # Track which paternal loci were matched
    merged = []
    
    # For each maternal locus, find matching paternal
    for m in mat_loci:
        best_match = None
        best_overlap = -1
        
        for j, p in enumerate(pat_by_chr.get(m['chr'], [])):
            if id(p) in pat_matched:
                continue
            # Check overlap with tolerance
            overlap_start = max(m['start'], p['start'])
            overlap_end = min(m['end'], p['end'])
            overlap = overlap_end - overlap_start
            
            if (m['start'] - dist <= p['end']) and (m['end'] + dist >= p['start']):
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_match = p
        
        if best_match is not None:
            pat_matched.add(id(best_match))
            # HOM: present in both
            struct_comp = compare_mt_structure(
                m.get('mt_structure', ''), 
                best_match.get('mt_structure', '')
            )
            cons_cat = consensus_category(
                m.get('category', ''), 
                best_match.get('category', '')
            )
            
            merged.append({
                'hg38_chr': m['chr'],
                'hg38_start': min(m['start'], best_match['start']),
                'hg38_end': max(m['end'], best_match['end']),
                'zygosity': 'HOM',
                'consensus_category': cons_cat,
                'mat_category': m.get('category', ''),
                'pat_category': best_match.get('category', ''),
                'structure_concordance': struct_comp,
                'mat_mt_structure': m.get('mt_structure', ''),
                'pat_mt_structure': best_match.get('mt_structure', ''),
                'mat_n_blocks': m.get('n_blocks', ''),
                'pat_n_blocks': best_match.get('n_blocks', ''),
                'mat_contig': m.get('contig', ''),
                'pat_contig': best_match.get('contig', ''),
                'Blast_Status': 'VALIDATED' if m.get('Blast_Status') == 'VALIDATED' or best_match.get('Blast_Status') == 'VALIDATED' else m.get('Blast_Status', ''),
            })
        else:
            # HET: maternal only
            merged.append({
                'hg38_chr': m['chr'],
                'hg38_start': m['start'],
                'hg38_end': m['end'],
                'zygosity': 'HET_MAT',
                'consensus_category': m.get('category', ''),
                'mat_category': m.get('category', ''),
                'pat_category': '.',
                'structure_concordance': '.',
                'mat_mt_structure': m.get('mt_structure', ''),
                'pat_mt_structure': '.',
                'mat_n_blocks': m.get('n_blocks', ''),
                'pat_n_blocks': '.',
                'mat_contig': m.get('contig', ''),
                'pat_contig': '.',
                'Blast_Status': m.get('Blast_Status', ''),
            })
    
    # Paternal-only (HET)
    for p in pat_loci:
        if id(p) not in pat_matched:
            merged.append({
                'hg38_chr': p['chr'],
                'hg38_start': p['start'],
                'hg38_end': p['end'],
                'zygosity': 'HET_PAT',
                'consensus_category': p.get('category', ''),
                'mat_category': '.',
                'pat_category': p.get('category', ''),
                'structure_concordance': '.',
                'mat_mt_structure': '.',
                'pat_mt_structure': p.get('mt_structure', ''),
                'mat_n_blocks': '.',
                'pat_n_blocks': p.get('n_blocks', ''),
                'mat_contig': '.',
                'pat_contig': p.get('contig', ''),
                'Blast_Status': p.get('Blast_Status', ''),
            })
    
    # Sort by chr, start
    chr_order = {f'chr{i}': i for i in range(1, 23)}
    chr_order['chrX'] = 23
    chr_order['chrY'] = 24
    chr_order['chrM'] = 25
    merged.sort(key=lambda x: (chr_order.get(x['hg38_chr'], 99), x['hg38_start']))
    
    return merged


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Merge maternal and paternal haplotypes.')
    parser.add_argument('--mat', required=True, help='Maternal liftover TSV')
    parser.add_argument('--pat', required=True, help='Paternal liftover TSV')
    parser.add_argument('--out', required=True, help='Output TSV')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    
    HEADER = [
        'hg38_chr', 'hg38_start', 'hg38_end', 'zygosity',
        'consensus_category', 'mat_category', 'pat_category',
        'structure_concordance',
        'mat_mt_structure', 'pat_mt_structure',
        'mat_n_blocks', 'pat_n_blocks',
        'mat_contig', 'pat_contig',
        'Blast_Status',
    ]
    
    mat_loci = parse_tsv(args.mat) if os.path.exists(args.mat) else []
    pat_loci = parse_tsv(args.pat) if os.path.exists(args.pat) else []
    
    merged = merge_haplotypes(mat_loci, pat_loci)
    
    with open(args.out, 'w') as f:
        f.write('\t'.join(HEADER) + '\n')
        for m in merged:
            f.write('\t'.join(str(m.get(h, '')) for h in HEADER) + '\n')


if __name__ == "__main__":
    main()
