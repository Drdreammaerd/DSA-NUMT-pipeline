#!/usr/bin/env python3
"""
Biological Validation of NUMT loci using BLAST.
Appends a 'Blast_Status' column (VALIDATED/FAILED) to the numt_classification.tsv.
"""
import os
import sys
import argparse
import subprocess
import pysam

def parse_tsv(tsv_path):
    rows = []
    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            while len(fields) < len(header):
                fields.append('')
            rows.append(dict(zip(header, fields)))
    return header, rows

def check_and_index(db_path, db_type):
    """Checks if BLAST database exists; if not, creates it."""
    ext = ".nsq" if db_type == "nucl" else ".psq"
    if not os.path.exists(db_path + ext):
        print(f"[*] BLAST index not found for {db_path}. Running makeblastdb...")
        cmd = ["makeblastdb", "-in", db_path, "-dbtype", db_type]
        subprocess.run(cmd, check=True)
    else:
        print(f"[*] Found BLAST index for {os.path.basename(db_path)}.")

def main():
    parser = argparse.ArgumentParser(description="BLAST Validation for NUMTs")
    parser.add_argument("--input_tsv", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--mode", default="dna", choices=["dna", "protein", "both"])
    args = parser.parse_args()

    # Fixed paths based on Docker / resources dir
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    chrm_db = os.path.join(project_dir, "resources", "chrM.fa")
    chrm_prot_db = os.path.join(project_dir, "resources", "mito_proteins.fa")
    
    # Auto-index
    if args.mode in ['dna', 'both']:
        check_and_index(chrm_db, "nucl")
    if args.mode in ['protein', 'both']:
        check_and_index(chrm_prot_db, "prot")

    header, rows = parse_tsv(args.input_tsv)
    
    if not rows:
        # Empty TSV pass-through
        with open(args.out_tsv, 'w') as f:
            if 'Blast_Status' not in header:
                header.append('Blast_Status')
            f.write('\t'.join(header) + '\n')
        print("[*] Empty TSV, passing through.")
        return

    # Extract sequences
    temp_dir = os.path.dirname(args.out_tsv)
    base_name = os.path.basename(args.out_tsv).replace(".tsv", "")
    temp_fa = os.path.join(temp_dir, f"{base_name}_extracted.fasta")
    
    print(f"[*] Extracting sequences from {os.path.basename(args.fasta)}...")
    with pysam.FastaFile(args.fasta) as fasta_ref, open(temp_fa, "w") as out_fa:
        for i, row in enumerate(rows):
            contig = row['contig']
            start = int(row['locus_start'])
            end = int(row['locus_end'])
            region_id = f"LOCUS_{i}"
            try:
                seq = fasta_ref.fetch(contig, start, end)
                out_fa.write(f">{region_id}\n{seq}\n")
            except KeyError:
                print(f"[WARNING] Contig {contig} not found in FASTA.")

    # Run BLAST
    dna_hits = set()
    prot_hits = set()
    
    if args.mode in ['dna', 'both']:
        out_dna = os.path.join(temp_dir, f"{base_name}_vs_chrM.blastn.txt")
        print(f"[*] Running blastn (DNA-DNA) with {args.threads} threads...")
        cmd_dna = [
            "blastn", "-query", temp_fa, "-db", chrm_db,
            "-out", out_dna, "-evalue", "1e-10", "-outfmt", "6",
            "-num_threads", str(args.threads)
        ]
        result = subprocess.run(cmd_dna, stdout=subprocess.DEVNULL)
        if result.returncode != 0:
            print(f"[WARNING] blastn exited with code {result.returncode}. All loci will be marked FAILED for DNA.")
        if os.path.exists(out_dna) and os.path.getsize(out_dna) > 0:
            with open(out_dna) as f:
                for line in f:
                    dna_hits.add(line.split()[0])
                    
    if args.mode in ['protein', 'both']:
        out_prot = os.path.join(temp_dir, f"{base_name}_vs_mitoProt.blastx.txt")
        print(f"[*] Running blastx (DNA-Protein) with {args.threads} threads...")
        cmd_prot = [
            "blastx", "-query", temp_fa, "-db", chrm_prot_db,
            "-query_gencode", "2",
            "-out", out_prot, "-evalue", "1e-5", "-outfmt", "6",
            "-num_threads", str(args.threads)
        ]
        result = subprocess.run(cmd_prot, stdout=subprocess.DEVNULL)
        if result.returncode != 0:
            print(f"[WARNING] blastx exited with code {result.returncode}. All loci will be marked FAILED for Protein.")
        if os.path.exists(out_prot) and os.path.getsize(out_prot) > 0:
            with open(out_prot) as f:
                for line in f:
                    prot_hits.add(line.split()[0])

    # Annotate TSV
    if 'Blast_Status' not in header:
        header.append('Blast_Status')
        
    passed = 0
    with open(args.out_tsv, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for i, row in enumerate(rows):
            region_id = f"LOCUS_{i}"
            is_dna = region_id in dna_hits
            is_prot = region_id in prot_hits
            
            if (args.mode == 'dna' and is_dna) or \
               (args.mode == 'protein' and is_prot) or \
               (args.mode == 'both' and (is_dna or is_prot)):
                status = "VALIDATED"
                passed += 1
            else:
                status = "FAILED"
                
            row['Blast_Status'] = status
            f.write('\t'.join(str(row.get(h, '')) for h in header) + '\n')

    # Cleanup temp fasta and blast output to save space
    try:
        os.remove(temp_fa)
        if args.mode in ['dna', 'both'] and os.path.exists(out_dna):
            os.remove(out_dna)
        if args.mode in ['protein', 'both'] and os.path.exists(out_prot):
            os.remove(out_prot)
    except:
        pass

    print(f"\n[*] BLAST Validation Complete:")
    print(f"    Total loci: {len(rows)}")
    print(f"    Validated:  {passed} ({(passed/len(rows))*100 if len(rows) else 0:.1f}%)")
    print(f"    Failed:     {len(rows) - passed}")

if __name__ == "__main__":
    main()
