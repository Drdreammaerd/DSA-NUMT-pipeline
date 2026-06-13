#!/usr/bin/env python3
###################################
# Purpose: Generate final callset from unified cross-sample merge.
#          Adds AF, QUAL flag, and outputs BED + VCF formats.
# Author: Yung-Chun Wang <yung-chun@wustl.edu>
# AI Assistant: Gemini
# Usage: python3 scripts/07_final_callset.py
###################################

import os
import re
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)

INDIVIDUALS = []
DONOR_FRACTIONS = {}
N_HAPLOTYPES = 0


# ── Load ──

def load_unified(path):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            while len(fields) < len(header):
                fields.append('')
            rows.append(dict(zip(header, fields)))
    return rows


# ── AF Calculation ──

def calc_af(row):
    """Calculate weighted allele frequency.
    
    If custom donor fractions are provided in the sample sheet, AF is weighted:
    HOM = donor_fraction, HET = donor_fraction / 2, ABSENT = 0.
    If no fractions are provided, falls back to uniform AF: 1 / N_INDIVIDUALS per donor.
    """
    allele_count = 0
    weighted_af = 0.0
    
    for indiv in INDIVIDUALS:
        zyg = row.get(f'{indiv}_zygosity', 'ABSENT')
        donor_fraction = DONOR_FRACTIONS.get(indiv, 0.0)
        
        if zyg == 'HOM':
            allele_count += 2
            weighted_af += donor_fraction
        elif zyg in ('HET_MAT', 'HET_PAT'):
            allele_count += 1
            weighted_af += donor_fraction / 2.0
    
    return allele_count, weighted_af


# ── Evidence Metrics ──

def calc_evidence(row):
    """Calculate donor-level evidence metrics.
    
    Returns dict with:
        n_hom:  Number of donors where NUMT is HOM (both haplotypes)
        n_het:  Number of donors where NUMT is HET (one haplotype only)
    
    These are raw metrics at the donor level.
    AF (allele frequency) is calculated separately at the haplotype level.
    """
    n_hom = 0
    n_het = 0
    for indiv in INDIVIDUALS:
        zyg = row.get(f'{indiv}_zygosity', 'ABSENT')
        if zyg == 'HOM':
            n_hom += 1
        elif zyg in ('HET_MAT', 'HET_PAT'):
            n_het += 1
    
    return {'n_hom': n_hom, 'n_het': n_het}


# ── HG005 Status ──

def calc_hg005_status(row):
    """Classify each NUMT relative to HG005 (HapMap reference sample).
    
    From the HapMap validation perspective:
      - HG005_germline:  Present in HG005 → detectable in HG005 short reads
      - other_only:      Absent from HG005 → germline in other donors,
                         not expected in HG005 data
    
    Note: ALL NUMTs in HPRC assemblies are germline. 'other_only' means
    germline in other individuals, NOT somatic.
    """
    hg005_zyg = row.get('HG005_zygosity', 'ABSENT')
    if hg005_zyg in ('HOM', 'HET_MAT', 'HET_PAT'):
        return 'HG005_germline'
    else:
        return 'other_only'


# ── chrM coverage range ──

def get_chrm_range(row):
    """Extract the overall chrM coordinate range from mt_structures across individuals."""
    mt_min = float('inf')
    mt_max = 0
    
    for indiv in INDIVIDUALS:
        mt_str = row.get(f'{indiv}_mt_structure', '.')
        if not mt_str or mt_str == '.':
            continue
        for part in mt_str.split(';'):
            m = re.match(r'(\d+)-(\d+)', part.strip())
            if m:
                mt_min = min(mt_min, int(m.group(1)))
                mt_max = max(mt_max, int(m.group(2)))
    
    if mt_min == float('inf'):
        return '.', '.', 0
    return str(mt_min), str(mt_max), mt_max - mt_min


# ── Output Writers ──

def write_final_tsv(rows, path):
    """Write final annotated TSV."""
    header = [
        'hg38_chr', 'hg38_start', 'hg38_end', 'locus_size',
        'allele_count', 'allele_freq',
        'n_individuals', 'n_hom', 'n_het',
        'hg005_status',
        'present_in', 'blast_validation',
        'max_category', 'is_complex', 'category_vote',
        'structural_polymorphism', 'cross_sample_concordance',
        'is_ref_numt', 'ref_category',
        'chrM_start', 'chrM_end', 'chrM_span',
    ]
    for indiv in INDIVIDUALS:
        header.extend([f'{indiv}_zygosity', f'{indiv}_category'])
    
    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(h, '.')) for h in header) + '\n')
    
    return header


def write_bed(rows, path):
    """Write BED6 format (0-based coords)."""
    with open(path, 'w') as f:
        f.write('#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n')
        for i, r in enumerate(rows):
            chrom = r['hg38_chr']
            start = int(r['hg38_start'])  # BED is 0-based
            end = int(r['hg38_end'])
            name = f"NUMT_{i+1:04d}_{r['max_category']}"
            # Score: map AF to 0-1000
            score = int(float(r['allele_freq']) * 1000)
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\n")


def write_vcf(rows, path):
    """Write VCF 4.2 format."""
    with open(path, 'w') as f:
        # Header
        f.write('##fileformat=VCFv4.2\n')
        f.write('##source=NUMT-CNV_pipeline\n')
        f.write('##reference=GRCh38/hg38\n')
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n')
        f.write(f'##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency (weighted by sample fraction if provided)">\n')
        f.write(f'##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count (out of {N_HAPLOTYPES} haplotypes)">\n')
        f.write('##INFO=<ID=NHOM,Number=1,Type=Integer,Description="Number of donors homozygous for NUMT">\n')
        f.write('##INFO=<ID=NHET,Number=1,Type=Integer,Description="Number of donors heterozygous for NUMT">\n')
        f.write('##INFO=<ID=NUMT_CAT,Number=1,Type=String,Description="NUMT structural category (A-E)">\n')
        f.write('##INFO=<ID=IS_COMPLEX,Number=0,Type=Flag,Description="Complex NUMT (C/D/E)">\n')
        f.write('##INFO=<ID=IS_REF,Number=0,Type=Flag,Description="Present in hg38 reference">\n')
        f.write('##INFO=<ID=BLAST_VAL,Number=0,Type=Flag,Description="Validated by BLAST">\n')
        f.write('##INFO=<ID=STRUCT_POLY,Number=1,Type=String,Description="Structural polymorphism flag (NO/WITHIN_TIER/YES)">\n')
        f.write('##INFO=<ID=CHRM_RANGE,Number=1,Type=String,Description="chrM coordinate range">\n')
        f.write('##INFO=<ID=N_INDIV,Number=1,Type=Integer,Description="Number of individuals with this NUMT">\n')
        f.write('##INFO=<ID=HG005_STATUS,Number=1,Type=String,Description="HG005_germline or other_only">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        # Column header
        sample_cols = '\t'.join(INDIVIDUALS)
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_cols}\n")
        
        # Data rows
        for i, r in enumerate(rows):
            chrom = r['hg38_chr']
            pos = r['hg38_start']
            vid = f"NUMT_{i+1:04d}"
            ref = 'N'
            alt = '<NUMT>'
            qual = '.'
            filt = '.'
            
            # INFO
            svlen = r['locus_size']
            info_parts = [
                'SVTYPE=NUMT',
                f"END={r['hg38_end']}",
                f"SVLEN={svlen}",
                f"AF={r['allele_freq']}",
                f"AC={r['allele_count']}",
                f"NHOM={r['n_hom']}",
                f"NHET={r['n_het']}",
                f"NUMT_CAT={r['max_category']}",
                f"STRUCT_POLY={r['structural_polymorphism']}",
                f"CHRM_RANGE={r['chrM_start']}-{r['chrM_end']}",
                f"N_INDIV={r['n_individuals']}",
                f"HG005_STATUS={r['hg005_status']}",
            ]
            if r['is_complex'] == 'YES':
                info_parts.append('IS_COMPLEX')
            if r['is_ref_numt'] == 'YES':
                info_parts.append('IS_REF')
            if r.get('blast_validation') == 'VALIDATED':
                info_parts.append('BLAST_VAL')
            info = ';'.join(info_parts)
            
            # Genotypes per individual
            gts = []
            for indiv in INDIVIDUALS:
                zyg = r.get(f'{indiv}_zygosity', 'ABSENT')
                if zyg == 'HOM':
                    gts.append('1/1')
                elif zyg in ('HET_MAT', 'HET_PAT'):
                    gts.append('0/1')
                else:
                    gts.append('0/0')
            gt_str = '\t'.join(gts)
            
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\tGT\t{gt_str}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_sheet", required=True)
    parser.add_argument("--input_tsv", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    global INDIVIDUALS
    global DONOR_FRACTIONS
    global N_HAPLOTYPES
    
    with open(args.sample_sheet) as f:
        header = f.readline().strip().split('\t')
        fraction_idx = header.index('Fraction') if 'Fraction' in header else -1
        
        individuals = set()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 0 and parts[0]:
                donor = parts[0]
                individuals.add(donor)
                if fraction_idx != -1 and len(parts) > fraction_idx:
                    try:
                        DONOR_FRACTIONS[donor] = float(parts[fraction_idx])
                    except ValueError:
                        pass

    INDIVIDUALS.extend(sorted(list(individuals)))
    N_HAPLOTYPES = len(INDIVIDUALS) * 2
    
    # If no valid fractions were provided, fallback to uniform distribution
    if not DONOR_FRACTIONS:
        uniform_fraction = 1.0 / len(INDIVIDUALS) if INDIVIDUALS else 0
        for indiv in INDIVIDUALS:
            DONOR_FRACTIONS[indiv] = uniform_fraction

    os.makedirs(args.outdir, exist_ok=True)
    
    print("=" * 70)
    print("  Step 5: Final Callset Generation")
    print("=" * 70)
    
    # Load unified callset
    rows = load_unified(args.input_tsv)
    print(f"  Loaded {len(rows)} unified loci")
    
    # Annotate each locus
    for r in rows:
        # AF
        ac, af = calc_af(r)
        r['allele_count'] = ac
        r['allele_freq'] = f"{af:.4f}"
        
        # Evidence metrics (replaces subjective QUAL)
        evidence = calc_evidence(r)
        r.update(evidence)
        
        # HG005 status
        r['hg005_status'] = calc_hg005_status(r)
        
        # chrM range
        mt_start, mt_end, mt_span = get_chrm_range(r)
        r['chrM_start'] = mt_start
        r['chrM_end'] = mt_end
        r['chrM_span'] = mt_span
    
    # Stats
    print(f"\n{'─' * 70}")
    print(f"  Final Callset Statistics")
    print(f"{'─' * 70}")
    
    # AF distribution
    af_buckets = [(0, 0.1), (0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
    print(f"\n  Allele Frequency distribution:")
    for lo, hi in af_buckets:
        n = sum(1 for r in rows if lo <= float(r['allele_freq']) < hi)
        bar = '█' * (n // 5)
        label = f"{lo:.1f}-{hi:.1f}" if hi <= 1 else f"{lo:.1f}-1.0"
        print(f"    AF {label}: {n:>4} {bar}")
    
    # Evidence metrics
    print(f"\n  Donor-level evidence:")
    n_indiv = len(INDIVIDUALS)
    for label, condition in [
        (f'{n_indiv}/{n_indiv} donors, all HOM', lambda r: int(r['n_individuals']) == n_indiv and int(r['n_hom']) == n_indiv),
        (f'{n_indiv}/{n_indiv} donors, some HET', lambda r: int(r['n_individuals']) == n_indiv and int(r['n_hom']) < n_indiv),
        (f'2-{n_indiv-1} donors', lambda r: 2 <= int(r['n_individuals']) <= n_indiv - 1),
        ('1 donor, HOM', lambda r: int(r['n_individuals']) == 1 and int(r['n_hom']) == 1),
        ('1 donor, HET only', lambda r: int(r['n_individuals']) == 1 and int(r['n_hom']) == 0),
    ]:
        n = sum(1 for r in rows if condition(r))
        print(f"    {label:<25}: {n:>4} ({n*100/len(rows):.1f}%)")
    
    # Complex summary
    complex_rows = [r for r in rows if r['is_complex'] == 'YES']
    print(f"\n  Complex NUMTs: {len(complex_rows)}")
    print(f"    ≥1 HOM donor: {sum(1 for r in complex_rows if int(r['n_hom']) >= 1)}")
    print(f"    HET only:    {sum(1 for r in complex_rows if int(r['n_hom']) == 0)}")
    
    # HG005 status breakdown
    n_hg005 = sum(1 for r in rows if r['hg005_status'] == 'HG005_germline')
    n_other = sum(1 for r in rows if r['hg005_status'] == 'other_only')
    print(f"\n  HG005 status (HapMap perspective):")
    print(f"    HG005_germline: {n_hg005:>4} ({n_hg005*100/len(rows):.1f}%) — present in HG005")
    print(f"    other_only:     {n_other:>4} ({n_other*100/len(rows):.1f}%) — absent from HG005")
    
    # Cross-tabulate with complexity
    print(f"\n    HG005_germline breakdown:")
    hg005_rows = [r for r in rows if r['hg005_status'] == 'HG005_germline']
    other_rows = [r for r in rows if r['hg005_status'] == 'other_only']
    for cat_label, subset in [('HG005_germline', hg005_rows), ('other_only', other_rows)]:
        n_simple = sum(1 for r in subset if r['is_complex'] == 'NO')
        n_complex = sum(1 for r in subset if r['is_complex'] == 'YES')
        n_ref = sum(1 for r in subset if r['is_ref_numt'] == 'YES')
        n_novel = sum(1 for r in subset if r['is_ref_numt'] == 'NO')
        print(f"      {cat_label:<18}: {len(subset):>4} total | "
              f"simple={n_simple} complex={n_complex} | ref={n_ref} novel={n_novel}")
    
    # Write outputs
    tsv_path = os.path.join(args.outdir, "HPRC_NUMT_callset.tsv")
    bed_path = os.path.join(args.outdir, "HPRC_NUMT_callset.bed")
    vcf_path = os.path.join(args.outdir, "HPRC_NUMT_callset.vcf")
    
    write_final_tsv(rows, tsv_path)
    write_bed(rows, bed_path)
    write_vcf(rows, vcf_path)
    
    print(f"\n  Output files:")
    print(f"    TSV: {tsv_path}")
    print(f"    BED: {bed_path}")
    print(f"    VCF: {vcf_path}")
    
    # Summary line counts
    print(f"\n  File sizes:")
    for p in [tsv_path, bed_path, vcf_path]:
        with open(p) as f:
            n_lines = sum(1 for _ in f)
        size = os.path.getsize(p)
        print(f"    {os.path.basename(p)}: {n_lines} lines, {size/1024:.1f} KB")
    
    print(f"\n  Done!")


if __name__ == "__main__":
    main()
