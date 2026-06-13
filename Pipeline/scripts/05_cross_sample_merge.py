#!/usr/bin/env python3
###################################
# Purpose: Cross-sample merge. Combine 6 per-individual callsets into a
#          unified NUMT catalogue with population frequency and structural
#          concordance across individuals.
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# AI Assistant: Gemini
# Usage: python3 scripts/06_cross_sample_merge.py
###################################

import os
import re
import argparse
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

INDIVIDUALS = []
MATCH_DIST = 2000
MT_COORD_TOLERANCE = 100
COMPLEX_CATS = {"C_inversion", "D_tandem_repeat", "E_complex_chimeric"}
CAT_RANK = {"A_single_block": 0, "B_divergence_gap": 1, "C_inversion": 2,
            "D_tandem_repeat": 3, "E_complex_chimeric": 4}


# ── Parsing ──

def parse_merged_tsv(path, indiv_name):
    """Parse per-individual merged TSV."""
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
            row['individual'] = indiv_name
            loci.append(row)
    return loci


def parse_ref_tsv(path):
    """Parse hg38 reference classification TSV."""
    loci = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            while len(fields) < len(header):
                fields.append('')
            row = dict(zip(header, fields))
            row['chr'] = row['contig']
            row['start'] = int(row['locus_start'])
            row['end'] = int(row['locus_end'])
            loci.append(row)
    return loci


def parse_mt_structure(mt_str):
    if not mt_str or mt_str == '.':
        return []
    segs = []
    for part in mt_str.split(';'):
        m = re.match(r'(\d+)-(\d+)\(([+-])\)', part.strip())
        if m:
            segs.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    return segs


# ── Clustering ──

def cluster_loci(all_loci, dist=MATCH_DIST):
    """Cluster loci from all individuals by genomic position.
    
    Uses single-linkage clustering: two loci are in the same cluster
    if their genomic intervals overlap (with tolerance).
    
    Returns list of clusters, each cluster is a list of loci dicts.
    """
    # Sort by chr, start
    chr_order = {f'chr{i}': i for i in range(1, 23)}
    chr_order['chrX'] = 23
    chr_order['chrY'] = 24
    chr_order['chrM'] = 25
    all_loci.sort(key=lambda x: (chr_order.get(x['chr'], 99), x['start']))
    
    clusters = []
    current_cluster = []
    cluster_chr = None
    cluster_end = -1
    
    for locus in all_loci:
        if locus['chr'] == cluster_chr and locus['start'] <= cluster_end + dist:
            # Extend current cluster
            current_cluster.append(locus)
            cluster_end = max(cluster_end, locus['end'])
        else:
            # Start new cluster
            if current_cluster:
                clusters.append(current_cluster)
            current_cluster = [locus]
            cluster_chr = locus['chr']
            cluster_end = locus['end']
    
    if current_cluster:
        clusters.append(current_cluster)
    
    return clusters


def match_to_ref(cluster_chr, cluster_start, cluster_end, ref_by_chr, dist=MATCH_DIST):
    """Find matching reference locus for a cluster."""
    for r in ref_by_chr.get(cluster_chr, []):
        if (cluster_start - dist <= r['end']) and (cluster_end + dist >= r['start']):
            return r
    return None


# ── Consensus ──

def get_consensus_category(categories):
    """Return the most complex category from a list."""
    best = None
    best_rank = -1
    for cat in categories:
        rank = CAT_RANK.get(cat, -1)
        if rank > best_rank:
            best_rank = rank
            best = cat
    return best


def normalize_mt_segments(mt_str):
    """Normalize mt_structure to a canonical form (sorted coords, ignore strand).
    Used for cross-sample comparison."""
    segs = parse_mt_structure(mt_str)
    if not segs:
        return ""
    # Sort by start coord, ignore strand for comparison
    coords = sorted([(s, e) for s, e, _ in segs])
    return ";".join(f"{s}-{e}" for s, e in coords)


def summarize_cluster(cluster, ref_match):
    """Summarize a cluster of loci from multiple individuals."""
    
    SIMPLE_TIER = {"A_single_block", "B_divergence_gap"}
    
    # Genomic coordinates (union)
    chrom = cluster[0]['chr']
    start = min(l['start'] for l in cluster)
    end = max(l['end'] for l in cluster)
    
    # Which individuals have this locus?
    indiv_data = {}
    for l in cluster:
        indiv = l['individual']
        if indiv not in indiv_data:
            indiv_data[indiv] = l
    
    n_individuals = len(indiv_data)
    present_in = ",".join(sorted(indiv_data.keys()))
    
    # Zygosity per individual
    zygosity_list = []
    for indiv in INDIVIDUALS:
        if indiv in indiv_data:
            zygosity_list.append(indiv_data[indiv].get('zygosity', '?'))
        else:
            zygosity_list.append('ABSENT')
    
    # Category per individual
    cat_list = []
    for indiv in INDIVIDUALS:
        if indiv in indiv_data:
            cat_list.append(indiv_data[indiv].get('consensus_category', '?'))
        else:
            cat_list.append('.')
    
    # Reference annotation
    if ref_match:
        is_ref = 'YES'
        ref_category = ref_match.get('category', '.')
        ref_mt_structure = ref_match.get('mt_structure', '.')
    else:
        is_ref = 'NO'
        ref_category = '.'
        ref_mt_structure = '.'
    
    # ── Category vote (includes ref) ──
    # Build vote: count occurrences of each category
    present_cats = [c for c in cat_list if c != '.']
    all_cats_for_vote = list(present_cats)
    if ref_category and ref_category != '.':
        all_cats_for_vote.append(ref_category)
    
    vote_counts = defaultdict(int)
    for c in all_cats_for_vote:
        vote_counts[c] += 1
    # Format: "5:B,1:A,ref:A" or "6:A,ref:A"
    vote_parts = []
    for cat in sorted(vote_counts.keys(), key=lambda x: -vote_counts[x]):
        # Count how many are from HPRC individuals vs ref
        hprc_count = sum(1 for c in present_cats if c == cat)
        ref_has = (ref_category == cat)
        parts = []
        if hprc_count > 0:
            parts.append(f"{hprc_count}:{cat}")
        if ref_has:
            parts.append(f"ref:{cat}")
        vote_parts.extend(parts)
    category_vote = ",".join(vote_parts) if vote_parts else "."
    
    # ── Structural polymorphism flag ──
    # Check if categories span the simple/complex boundary (including ref)
    all_present_cats = set(present_cats)
    if ref_category and ref_category != '.':
        all_present_cats.add(ref_category)
    
    has_simple = bool(all_present_cats & SIMPLE_TIER)
    has_complex = bool(all_present_cats & COMPLEX_CATS)
    
    if has_simple and has_complex:
        structural_polymorphism = 'YES'
    elif len(all_present_cats) > 1:
        structural_polymorphism = 'WITHIN_TIER'  # e.g., A vs B or C vs E
    else:
        structural_polymorphism = 'NO'
    
    # Overall max category (includes ref)
    max_cat = get_consensus_category(list(all_present_cats)) if all_present_cats else '.'
    
    # Is complex? (based on max)
    is_complex = max_cat in COMPLEX_CATS
    
    # Structural concordance across individuals (using mat_mt_structure as representative)
    mt_structures = []
    for indiv in INDIVIDUALS:
        if indiv in indiv_data:
            mat_mt = indiv_data[indiv].get('mat_mt_structure', '.')
            pat_mt = indiv_data[indiv].get('pat_mt_structure', '.')
            # Use whichever is non-empty, preferring mat
            mt = mat_mt if mat_mt and mat_mt != '.' else pat_mt
            mt_structures.append(mt)
        else:
            mt_structures.append('.')
    
    # Normalize and compare
    present_normalized = [normalize_mt_segments(m) for m in mt_structures if m and m != '.']
    unique_structures = set(present_normalized)
    if len(unique_structures) <= 1:
        cross_sample_concordance = 'IDENTICAL'
    elif len(unique_structures) <= 2:
        cross_sample_concordance = 'SIMILAR'
    else:
        cross_sample_concordance = 'VARIABLE'
    
    locus_size = end - start
    
    # Blast Status summary (VALIDATED if any individual has it)
    blast_statuses = [l.get('Blast_Status', 'FAILED') for l in cluster]
    if 'VALIDATED' in blast_statuses:
        cluster_blast_status = 'VALIDATED'
    else:
        cluster_blast_status = 'FAILED'
    
    return {
        'hg38_chr': chrom,
        'hg38_start': start,
        'hg38_end': end,
        'locus_size': locus_size,
        'n_individuals': n_individuals,
        'present_in': present_in,
        'max_category': max_cat,
        'is_complex': 'YES' if is_complex else 'NO',
        'category_vote': category_vote,
        'structural_polymorphism': structural_polymorphism,
        'cross_sample_concordance': cross_sample_concordance,
        'blast_validation': cluster_blast_status,
        'is_ref_numt': is_ref,
        'ref_category': ref_category,
        # Per-individual zygosity
        **{f'{indiv}_zygosity': z for indiv, z in zip(INDIVIDUALS, zygosity_list)},
        # Per-individual category
        **{f'{indiv}_category': c for indiv, c in zip(INDIVIDUALS, cat_list)},
        # Per-individual mt_structure
        **{f'{indiv}_mt_structure': m for indiv, m in zip(INDIVIDUALS, mt_structures)},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_sheet", required=True)
    parser.add_argument("--ref_tsv", required=True)
    parser.add_argument("--mode", default="annotate")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--indiv_dir", required=True)
    args = parser.parse_args()

    global INDIVIDUALS
    with open(args.sample_sheet) as f:
        header = f.readline()
        individuals = set()
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                individuals.add(parts[0])
    INDIVIDUALS.extend(sorted(list(individuals)))

    os.makedirs(args.outdir, exist_ok=True)
    
    print("=" * 70)
    print(f"  Phase 3b: Cross-Sample Merge ({len(INDIVIDUALS)} individuals → unified callset)")
    print("=" * 70)
    
    # Load all per-individual data
    all_loci = []
    for indiv in INDIVIDUALS:
        path = os.path.join(args.indiv_dir, f"{indiv}.merged_haplotypes.tsv")
        if not os.path.exists(path):
            print(f"[SKIP] {indiv}: file not found")
            continue
        loci = parse_merged_tsv(path, indiv)
        print(f"  Loaded {indiv}: {len(loci)} loci")
        all_loci.extend(loci)
    
    print(f"\n  Total loci across all individuals: {len(all_loci)}")
    
    # Load reference
    ref_loci = parse_ref_tsv(args.ref_tsv)
    ref_by_chr = defaultdict(list)
    for r in ref_loci:
        ref_by_chr[r['chr']].append(r)
    print(f"  hg38 reference: {len(ref_loci)} loci")
    
    # Cluster
    clusters = cluster_loci(all_loci)
    print(f"\n  Clustered into {len(clusters)} unified loci")
    
    # Summarize each cluster
    unified = []
    for cluster in clusters:
        chrom = cluster[0]['chr']
        start = min(l['start'] for l in cluster)
        end = max(l['end'] for l in cluster)
        ref_match = match_to_ref(chrom, start, end, ref_by_chr)
        summary = summarize_cluster(cluster, ref_match)
        unified.append(summary)
    
    # Filter mode
    if args.mode == "filter":
        filtered_unified = []
        for u in unified:
            # Exclude if it perfectly matches the reference and has no structural polymorphism across individuals
            if u['is_ref_numt'] == 'YES' and u['structural_polymorphism'] != 'YES':
                continue
            filtered_unified.append(u)
        print(f"  Filtered {len(unified) - len(filtered_unified)} reference NUMTs. {len(filtered_unified)} remain.")
        unified = filtered_unified
    
    # Stats
    print(f"\n{'─' * 70}")
    print(f"  Unified Callset Statistics")
    print(f"{'─' * 70}")
    
    # Frequency distribution
    freq_dist = defaultdict(int)
    for u in unified:
        freq_dist[u['n_individuals']] += 1
    
    print(f"\n  Population frequency:")
    for n in sorted(freq_dist.keys()):
        pct = freq_dist[n] * 100 / len(unified)
        bar = '█' * int(pct / 2)
        print(f"    {n}/{len(INDIVIDUALS)} individuals: {freq_dist[n]:>4} loci ({pct:5.1f}%) {bar}")
    
    # Category breakdown
    print(f"\n  Category breakdown (max_category):")
    cat_counts = defaultdict(int)
    for u in unified:
        cat_counts[u['max_category']] += 1
    for cat in sorted(cat_counts.keys(), key=lambda x: CAT_RANK.get(x, 99)):
        print(f"    {cat}: {cat_counts[cat]}")
    
    # Complex NUMTs
    complex_loci = [u for u in unified if u['is_complex'] == 'YES']
    print(f"\n  Complex NUMTs (C+D+E): {len(complex_loci)}")
    
    # Complex frequency
    complex_freq = defaultdict(int)
    for u in complex_loci:
        complex_freq[u['n_individuals']] += 1
    total_indivs = len(INDIVIDUALS)
    print(f"    Shared by {total_indivs}/{total_indivs}: {complex_freq.get(total_indivs, 0)}")
    print(f"    Shared by {total_indivs-1}/{total_indivs}: {complex_freq.get(total_indivs-1, 0)}")
    print(f"    Shared by <{total_indivs-1}:  {sum(v for k,v in complex_freq.items() if k < total_indivs-1)}")
    
    # Reference overlap
    n_ref = sum(1 for u in unified if u['is_ref_numt'] == 'YES')
    n_novel = sum(1 for u in unified if u['is_ref_numt'] == 'NO')
    print(f"\n  Reference overlap:")
    print(f"    Reference-known: {n_ref}")
    print(f"    Novel (not in hg38): {n_novel}")
    
    # Novel complex
    novel_complex = [u for u in unified if u['is_complex'] == 'YES' and u['is_ref_numt'] == 'NO']
    print(f"    Novel complex: {len(novel_complex)}")
    for nc in novel_complex:
        print(f"      {nc['hg38_chr']}:{nc['hg38_start']}-{nc['hg38_end']} "
              f"cat={nc['max_category']} n={nc['n_individuals']} "
              f"vote={nc['category_vote']}")
    
    # Structural polymorphism
    n_poly_yes = sum(1 for u in unified if u['structural_polymorphism'] == 'YES')
    n_poly_within = sum(1 for u in unified if u['structural_polymorphism'] == 'WITHIN_TIER')
    n_poly_no = sum(1 for u in unified if u['structural_polymorphism'] == 'NO')
    print(f"\n  Structural polymorphism (category vote including ref):")
    print(f"    NO:          {n_poly_no:>4} ({n_poly_no*100/len(unified):.1f}%) — all agree")
    print(f"    WITHIN_TIER: {n_poly_within:>4} ({n_poly_within*100/len(unified):.1f}%) — A↔B or C↔D↔E")
    print(f"    YES:         {n_poly_yes:>4} ({n_poly_yes*100/len(unified):.1f}%) — simple↔complex boundary")
    
    if n_poly_yes > 0:
        print(f"\n    ▶ Cross-tier polymorphism cases (simple↔complex):")
        for u in unified:
            if u['structural_polymorphism'] == 'YES':
                print(f"      {u['hg38_chr']}:{u['hg38_start']}-{u['hg38_end']} "
                      f"vote={u['category_vote']}")
    
    # Write output
    header_cols = [
        'hg38_chr', 'hg38_start', 'hg38_end', 'locus_size',
        'n_individuals', 'present_in',
        'max_category', 'is_complex', 'category_vote',
        'structural_polymorphism',
        'cross_sample_concordance',
        'blast_validation',
        'is_ref_numt', 'ref_category',
    ]
    # Per-individual columns
    for indiv in INDIVIDUALS:
        header_cols.extend([
            f'{indiv}_zygosity',
            f'{indiv}_category',
        ])
    # mt_structure columns last (they're wide)
    for indiv in INDIVIDUALS:
        header_cols.append(f'{indiv}_mt_structure')
    
    out_path = os.path.join(args.outdir, "unified_callset.tsv")
    with open(out_path, 'w') as f:
        f.write('\t'.join(header_cols) + '\n')
        for u in unified:
            f.write('\t'.join(str(u.get(h, '')) for h in header_cols) + '\n')
    
    # Also write a compact summary (no mt_structure)
    compact_cols = header_cols[:14] + [f'{i}_zygosity' for i in INDIVIDUALS] + [f'{i}_category' for i in INDIVIDUALS]
    compact_path = os.path.join(args.outdir, "unified_callset_compact.tsv")
    with open(compact_path, 'w') as f:
        f.write('\t'.join(compact_cols) + '\n')
        for u in unified:
            f.write('\t'.join(str(u.get(h, '')) for h in compact_cols) + '\n')
    
    print(f"\n  Output files:")
    print(f"    Full:    {out_path}")
    print(f"    Compact: {compact_path}")
    print(f"\n  Done!")


if __name__ == "__main__":
    main()

