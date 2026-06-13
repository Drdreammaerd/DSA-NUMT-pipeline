# DSA-NUMT-pipeline (v1.0)

## Overview

An end-to-end Snakemake pipeline for detecting polymorphic Nuclear Mitochondrial DNA insertions (NUMTs) directly from genome assemblies using Diploid Sequence Alignment (DSA). This pipeline leverages **Huang's LAST Method** ([Huang & Frith, 2026](#references)) for accurate sequence alignment and NUMT classification, and is capable of processing multi-sample and multi-haplotype cohorts (e.g., HPRC/SMAHT datasets).

The pipeline autonomously executes NUMT discovery, LAST alignments, biological validation (via BLAST), cross-haplotype merging, and cross-sample population catalog generation. It is fully containerized via Docker and optimized for both Cluster environments (LSF) and Local servers.

## Pipeline Architecture

```
Assembly FASTA ──→ Alignment ──→ Classification ──→ Validation ──→ Merge (Sample) ──→ Merge (Population)
                   (LAST)        (Huang's Method)     (BLAST)        (Maternal+Paternal)
```

| Stage | Output | Description |
|-------|--------|-------------|
| **1: Alignment** | `.maf` & `.paf` | Aligns genome assemblies against human `chrM` using LAST (Huang's Method). |
| **2: Classification** | `numt_classification.tsv` | Classifies alignments into NUMT structural categories (e.g., single block, divergence gap, complex) using custom structural logic. |
| **3: Validation** | `numt_classification.annotated.tsv` | Biologically validates NUMT sequences against `chrM` (blastn) and mitochondrial proteins (blastx). |
| **4: Liftover** | `numt_liftover.tsv` | Lifts NUMT loci from the specific assembly coordinates back to GRCh38. |
| **5: Haplotype Merge** | `sample.merged.tsv` | Merges maternal and paternal NUMT calls for a single donor to determine zygosity (HOM/HET). |
| **6: Callset** | `HPRC_NUMT_callset.tsv/vcf/bed`| Cross-sample integration to produce the final population frequency catalog. |

> [!NOTE]
> For a detailed explanation of all output columns, Huang's structural categories (A-E), and VCF attributes across the pipeline stages, please see the **[DSA-NUMT Data Dictionary](docs/NUMT_DSA_Dictionary.md)**.

## Quick Start

The pipeline relies on a single configuration file: `config.yaml`.

### 1. Configure the Environment

Copy the provided template to create your configuration file:
```bash
cp Pipeline/config_template/config_template.yaml config.yaml
```

#### `config.yaml` Settings
Open `config.yaml` and configure your paths. 

```yaml
# ==========================================
# 1. Project Specific Paths (Modify these)
# ==========================================
work_dir: "/path/to/your/working/directory"

# ==========================================
# 2. Input Datasets
# ==========================================
# Define each sample, providing paths to paternal and maternal assemblies
# as well as the paths to the chain files for liftover to GRCh38.
samples:
  HG002:
    paternal: "/path/to/HG002.paternal.fa"
    paternal_chain: "/path/to/HG002_paternal_to_hg38.chain"
    maternal: "/path/to/HG002.maternal.fa"
    maternal_chain: "/path/to/HG002_maternal_to_hg38.chain"

# ==========================================
# 3. Reference Catalog (Auto-Generation)
# ==========================================
# If 'prebuilt_numts' is empty, the pipeline will automatically run DSA 
# on 'target_fasta' to generate the reference NUMT catalog from scratch.
paths:
  target_fasta: "/path/to/hg38.fa"
  chain_dir: "/path/to/chains"

reference_comparison:
  prebuilt_numts: "" 
  mode: "annotate"

# ==========================================
# 4. System Settings (Do not change)
# ==========================================
docker_image: "dreammaerd/numt_dsa:v1"
```

### 2. Execution

The pipeline supports two modes of execution: **Cluster Mode (LSF)** and **Local Mode (Standalone Server)**.

#### Option A: WashU RIS Cluster (LSF)
If you are on the WashU RIS cluster, submit the Snakemake orchestrator. It will run in the background and spawn LSF sub-jobs automatically:

```bash
# Run the full end-to-end pipeline
bash Pipeline/run_numt_dsa_pipeline.sh --configfile config.yaml --subscription
```

#### Option B: Local Machine / Standalone Server
If you are running on a local workstation, a powerful server, or a Mac, you do not need LSF. You can run the entire pipeline directly via Docker using your local CPU cores.

```bash
# Run the full pipeline locally
bash Pipeline/run_local.sh --configfile config.yaml
```
*Note: If your data is located on custom drives, they will be accessible as long as you launch the script from a directory that shares the mount path, or you manually add `-v` volume mounts to `run_local.sh`.*

## Dependencies & Architecture

- **Containerized:** The pipeline is entirely encapsulated within the `dreammaerd/numt_dsa:v1` Docker image. This includes Snakemake, Python, R, LAST, and BLAST.
- **Pre-built Indices:** The Docker image has `chrM.fa` and `mito_proteins.fa` indices permanently baked in, preventing read-only filesystem errors during execution.
- **No Local Dependencies:** Zero local packages are required other than Docker.

## References

1. **Huang's LAST Method**:
   Huang M, Frith MC. *Probability-Based Sequence Comparison Finds Pre-Eutherian Nuclear Mitochondrial DNA Segments in Mammalian Genomes.* J Comput Biol. 2026 Feb 2:15578666261416560. doi: [10.1177/15578666261416560](https://doi.org/10.1177/15578666261416560).
   GitHub Repository: [https://github.com/Koumokuyou/NUMTs](https://github.com/Koumokuyou/NUMTs)
