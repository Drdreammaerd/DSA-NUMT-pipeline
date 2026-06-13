#!/usr/bin/env python3
"""
Classify ALL NUMT loci from LAST MAF/PAF output.

Reads the full PAF file and classifies every NUMT locus into:
  A. Single-block       — 1 alignment block, simplest NUMT
  B. Divergence-gap     — multi-block, same strand, continuous chrM
  C. Inversion          — mixed +/- strands
  D. Tandem-repeat      — same chrM coordinates repeated
  E. Complex-chimeric   — different chrM regions, same strand

Output: TSV with one row per locus + summary statistics.

Author: Yung-Chun Wang <yung-chun@wustl.edu>
AI Assistant: Gemini
Usage:
  python3 01_classify_all_numts.py                       # default HG002.maternal
  python3 01_classify_all_numts.py --sample HG005.paternal
  python3 01_classify_all_numts.py --paf data/paf/custom.paf
"""

import os
import sys
import argparse
from collections import defaultdict, Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)


def parse_bed(bed_path):
    """Read BED file and return dict mapping contig to list of (start, end)."""
    bed_intervals = defaultdict(list)
    if bed_path and os.path.exists(bed_path):
        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    bed_intervals[parts[0]].append((int(parts[1]), int(parts[2])))
    return bed_intervals


def parse_paf(paf_path):
    """Read PAF file and return list of alignment dicts."""
    alns = []
    with open(paf_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            alns.append({
                'contig': parts[0],
                'q_len': int(parts[1]),
                'q_start': int(parts[2]),
                'q_end': int(parts[3]),
                'strand': parts[4],
                't_name': parts[5],
                't_len': int(parts[6]),
                't_start': int(parts[7]),
                't_end': int(parts[8]),
                'n_match': int(parts[9]),
                'aln_len': int(parts[10]),
                'mapq': int(parts[11]),
                'mt_size': int(parts[8]) - int(parts[7]),
            })
    return alns


def cluster_loci(alns, gap=2000):
    """Cluster alignments on the same contig within 'gap' bp."""
    by_contig = defaultdict(list)
    for a in alns:
        by_contig[a['contig']].append(a)

    clusters = []
    for contig, hits in by_contig.items():
        hits.sort(key=lambda x: x['q_start'])
        cluster = [hits[0]]
        for h in hits[1:]:
            if h['q_start'] - cluster[-1]['q_end'] < gap:
                cluster.append(h)
            else:
                clusters.append(cluster)
                cluster = [h]
        clusters.append(cluster)

    return clusters


def classify_locus(cluster):
    """
    Classify a single NUMT locus.

    Returns: (category, detail_dict)
    """
    n = len(cluster)
    strands = set(h['strand'] for h in cluster)
    total_mt = sum(h['mt_size'] for h in cluster)
    span = cluster[-1]['q_end'] - cluster[0]['q_start']

    # Check for tandem repeats (identical chrM coordinates)
    mt_coords = [(h['t_start'], h['t_end']) for h in cluster]
    coord_counter = Counter(mt_coords)
    has_tandem = any(v > 1 for v in coord_counter.values())
    max_repeat = max(coord_counter.values()) if coord_counter else 1

    # Check chrM continuity
    mt_sorted = sorted(mt_coords)
    is_continuous = True
    max_mt_gap = 0
    if n >= 2:
        for i in range(len(mt_sorted) - 1):
            gap = mt_sorted[i + 1][0] - mt_sorted[i][1]
            max_mt_gap = max(max_mt_gap, abs(gap))
            if abs(gap) > 500:
                is_continuous = False

    # Nuclear gaps
    sorted_by_q = sorted(cluster, key=lambda x: x['q_start'])
    nuclear_gaps = []
    for i in range(len(sorted_by_q) - 1):
        nuclear_gaps.append(sorted_by_q[i + 1]['q_start'] - sorted_by_q[i]['q_end'])

    # Classify
    if n == 1:
        category = 'A_single_block'
    elif has_tandem:
        category = 'D_tandem_repeat'
    elif len(strands) > 1:
        category = 'C_inversion'
    elif n >= 2 and is_continuous:
        category = 'B_divergence_gap'
    else:
        category = 'E_complex_chimeric'

    # Build chrM structure string
    mt_structure_parts = []
    for h in sorted_by_q:
        mt_structure_parts.append(f"{h['t_start']}-{h['t_end']}({h['strand']})")
    mt_structure = ';'.join(mt_structure_parts)

    # Tandem detail
    tandem_detail = ''
    if has_tandem:
        for coord, count in coord_counter.items():
            if count > 1:
                tandem_detail = f"chrM:{coord[0]}-{coord[1]}x{count}"
                break

    detail = {
        'contig': cluster[0]['contig'],
        'locus_start': cluster[0]['q_start'],
        'locus_end': cluster[-1]['q_end'],
        'locus_span': span,
        'n_blocks': n,
        'strands': '/'.join(sorted(strands)),
        'total_mt_coverage': total_mt,
        'mt_structure': mt_structure,
        'max_mt_gap': max_mt_gap,
        'nuclear_gaps': ','.join(str(g) for g in nuclear_gaps) if nuclear_gaps else '.',
        'max_score': max(h['mapq'] for h in cluster),
        'tandem_detail': tandem_detail,
        'category': category,
    }

    return category, detail


def main():
    parser = argparse.ArgumentParser(description='Classify all NUMT loci from PAF file.')
    parser.add_argument('--sample', default='HG002.maternal',
                        help='Sample name (default: HG002.maternal)')
    parser.add_argument('--paf', default=None,
                        help='Path to PAF file (overrides --sample)')
    parser.add_argument('--bed', default=None,
                        help='Path to BED file for filtering true NUMTs')
    parser.add_argument('--gap', type=int, default=2000,
                        help='Clustering gap in bp (default: 2000)')
    parser.add_argument('--outdir', default=None,
                        help='Output directory (default: PROJECT_DIR/results/classification/)')
    args = parser.parse_args()

    # Input PAF
    if args.paf:
        paf_path = args.paf
    else:
        paf_path = os.path.join(PROJECT_DIR, 'data', 'paf', f'{args.sample}.nu2mitogeno.paf')

    if not os.path.exists(paf_path):
        print(f"[ERROR] PAF not found: {paf_path}")
        print(f"  Run 02b_maf_to_paf_python.py first to generate PAF files.")
        sys.exit(1)

    # Output dir
    out_dir = args.outdir or os.path.join(PROJECT_DIR, 'results', 'classification')
    os.makedirs(out_dir, exist_ok=True)

    sample_name = args.sample

    print("=" * 70)
    print(f"  NUMT Locus Classification: {sample_name}")
    print(f"  Input: {paf_path}")
    print("=" * 70)
    print()

    # Parse and cluster
    alns = parse_paf(paf_path)
    print(f"Total alignment blocks: {len(alns)}")

    clusters = cluster_loci(alns, gap=args.gap)
    print(f"Clustered NUMT loci (gap={args.gap}bp): {len(clusters)}")

    # Filter by BED
    if args.bed:
        bed_intervals = parse_bed(args.bed)
        filtered_clusters = []
        for cluster in clusters:
            contig = cluster[0]['contig']
            c_start = cluster[0]['q_start']
            c_end = cluster[-1]['q_end']
            
            overlaps = False
            for b_start, b_end in bed_intervals.get(contig, []):
                if max(c_start, b_start) < min(c_end, b_end):
                    overlaps = True
                    break
            if overlaps:
                filtered_clusters.append(cluster)
        print(f"Filtered NUMT loci (overlapping BED): {len(filtered_clusters)}")
        clusters = filtered_clusters
    print()

    # Classify all
    results = []
    category_counts = Counter()

    for cluster in clusters:
        category, detail = classify_locus(cluster)
        category_counts[category] += 1
        results.append(detail)

    # Print summary
    total = len(results)
    cat_labels = {
        'A_single_block': 'A. Single-block',
        'B_divergence_gap': 'B. Divergence-gap (same strand, continuous chrM)',
        'C_inversion': 'C. Inversion (mixed strands)',
        'D_tandem_repeat': 'D. Tandem repeat (same chrM repeated)',
        'E_complex_chimeric': 'E. Complex chimeric (different chrM regions)',
    }

    print("-" * 70)
    for key in ['A_single_block', 'B_divergence_gap', 'C_inversion', 'D_tandem_repeat', 'E_complex_chimeric']:
        count = category_counts.get(key, 0)
        pct = count / total * 100 if total > 0 else 0
        print(f"  {cat_labels[key]:<55} {count:>5}  ({pct:.1f}%)")
    print("-" * 70)

    simple = category_counts.get('A_single_block', 0) + category_counts.get('B_divergence_gap', 0)
    complex_sv = total - simple
    print(f"  Simple/Normal (A+B): {simple:>5}  ({simple / total * 100:.1f}%)")
    print(f"  Complex SV (C+D+E): {complex_sv:>5}  ({complex_sv / total * 100:.1f}%)")
    print(f"  {'Total:':<55} {total:>5}")
    print()

    # Write TSV
    tsv_path = os.path.join(out_dir, f'{sample_name}.numt_classification.tsv')
    cols = ['contig', 'locus_start', 'locus_end', 'locus_span', 'n_blocks',
            'strands', 'total_mt_coverage', 'category', 'mt_structure',
            'max_mt_gap', 'nuclear_gaps', 'max_score', 'tandem_detail']

    with open(tsv_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in sorted(results, key=lambda x: (x['category'], -x['total_mt_coverage'])):
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    print(f"Full classification: {tsv_path}")

    # Write complex-only TSV
    complex_path = os.path.join(out_dir, f'{sample_name}.complex_sv_numts.tsv')
    with open(complex_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in sorted(results, key=lambda x: (x['category'], -x['total_mt_coverage'])):
            if r['category'] in ('C_inversion', 'D_tandem_repeat', 'E_complex_chimeric'):
                f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    print(f"Complex SV only:     {complex_path}")
    print()
    print("Done!")


if __name__ == "__main__":
    main()
