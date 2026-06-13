# Data Dictionary — DSA-NUMT Pipeline

This data dictionary explains the output formats and column definitions across all critical stages of the Diploid Sequence Alignment (DSA) NUMT pipeline:

1. **Stage 1: Classification** (Custom Structural Classification from LAST output)
2. **Stage 3: Haplotype Merge** (Per-individual zygosity determination)
3. **Stage 5: Final Callset** (Cross-sample population frequency and VCF)

---

## Stage 1: Classification & Validation

### `results/01_classification/<sample>.numt_classification.annotated.tsv`

This is the foundational file containing all NUMTs identified in a single haplotype assembly (e.g., `HG002.paternal`), structurally classified using custom logic based on LAST alignments, and biologically validated via BLAST.

| Column | Meaning |
|--------|---------|
| `nu_chr` | The target nuclear chromosome (e.g., `chr1`). |
| `nu_start` | The estimated start position on the nuclear chromosome. |
| `nu_end` | The estimated end position on the nuclear chromosome. |
| `locus_size` | Length of the inserted NUMT on the nuclear reference (`nu_end - nu_start`). |
| `chrM_start` | Start coordinate on the mitochondrial genome. |
| `chrM_end` | End coordinate on the mitochondrial genome. |
| `chrM_span` | Span on the mitochondrial genome. Can be negative if inserted in reverse orientation. |
| `numt_category` | Structural category assigned by Huang's logic. *See NUMT Category Definitions below.* |
| `strand` | Orientation of the insertion (`+` or `-`). |
| `valid_blastn` | Validated by DNA sequence homology (`blastn`) against `chrM` (`YES`/`NO`). |
| `valid_blastx` | Validated by Protein sequence homology (`blastx`) against `mito_proteins` (`YES`/`NO`). |
| `is_complex` | Flag indicating highly complex structural variants (Categories C, D, E). |

**NUMT Category Definitions (Custom Structural Logic):**

| Category | Definition |
|----------|------------|
| `A_single_block` | A clean, contiguous sequence inserted from the mitochondria. |
| `B_divergence_gap` | A single mitochondrial insertion that contains a small gap or divergence, causing LAST to split the alignment slightly. |
| `C_distal_block` | Two structurally distinct blocks originating from distant regions of the mitochondrial genome inserted adjacently. |
| `D_inversion` | Two structurally distinct blocks from the mitochondria, where one is inserted in the opposite orientation. |
| `E_duplication` | The same mitochondrial region inserted multiple times in tandem at the same locus. |

---

## Stage 3: Haplotype Merge (Per-Individual)

### `results/03_per_individual/<donor_id>.merged_haplotypes.tsv`

This file merges the `maternal` and `paternal` calls for a single donor to determine zygosity and structural consensus.

| Column | Meaning |
|--------|---------|
| `hg38_chr` | The target nuclear chromosome (lifted over to GRCh38). |
| `hg38_start` | The GRCh38 start position. |
| `hg38_end` | The GRCh38 end position. |
| `zygosity` | Zygosity of the insertion (`HOM` = Homozygous, `HET_MAT` = Heterozygous Maternal, `HET_PAT` = Heterozygous Paternal). |
| `mat_category` | Structural category observed on the maternal haplotype. |
| `pat_category` | Structural category observed on the paternal haplotype. |
| `mat_chrM_range` | Mitochondrial coordinates observed on the maternal haplotype. |
| `pat_chrM_range` | Mitochondrial coordinates observed on the paternal haplotype. |
| `structural_polymorphism` | Flag indicating if the maternal and paternal alleles are structurally different. *See Structural Polymorphism Definitions below.* |

**Structural Polymorphism Definitions:**

| Flag | Definition |
|------|------------|
| `NO` | Both alleles share the exact same structural category (e.g., both are `A_single_block`). |
| `WITHIN_TIER` | Both alleles share the same fundamental tier (e.g., `A` vs `B`), suggesting a minor divergence rather than a completely different insertion event. |
| `YES` | Alleles represent fundamentally different structural events (e.g., `A_single_block` vs `D_inversion`). |

---

## Stage 5: Final Callset (Cross-Sample Population Catalog)

Outputs are located in `results/05_final_callset/`.

### 1. `HPRC_NUMT_callset.tsv`

The master catalog across all donors in the cohort.

| Column | Meaning |
|--------|---------|
| `hg38_chr` / `start` / `end`| Consensus GRCh38 coordinates spanning all donors. |
| `locus_size` | Length of the consensus window. |
| `allele_count` | Number of haplotypes (alleles) containing this NUMT. |
| `allele_freq` | Allele Frequency (`allele_count / (N_individuals * 2)`). |
| `n_individuals` | Number of donors carrying this NUMT. |
| `n_hom` / `n_het` | Number of donors Homozygous and Heterozygous. |
| `is_ref_numt` | Flag indicating if this insertion is already present in the GRCh38 Reference (`YES`/`NO`). `NO` denotes a Novel/Non-reference NUMT. |
| `max_category` | The most severe structural category observed for this NUMT across the entire population. |
| `category_vote` | A comma-separated count of all categories observed across haplotypes (e.g., `A_single_block:10, B_divergence_gap:2`). |
| `cross_sample_concordance` | Flag indicating structural agreement across the population. (`CONCORDANT` or `DISCORDANT`). |
| `<donor>_zygosity` | Zygosity state for each specific donor. |
| `<donor>_category` | The structural category observed in each specific donor. |

### 2. `HPRC_NUMT_callset.vcf`

A standard, multi-sample VCF representing the population catalog. Suitable for genome browsers (IGV) and downstream variant tools.

**VCF Columns & FILTER:**

| Field / Column | Meaning |
|----------------|---------|
| `<donor_id>` | Sample format (e.g., `HG002`). |
| `FILTER` | Always `PASS` (assuming it survived mapping constraints). |
| `ALT` | `<INS:MT>`, signifying a Nuclear Mitochondrial Insertion. |

**INFO Fields:**

| Field | Meaning |
|-------|---------|
| `END` | End position on the nuclear reference. |
| `SVLEN` | Length of the inserted sequence. |
| `AF` | Allele frequency (weighted by sample fraction). |
| `AC` | Allele count. |
| `NHOM` / `NHET`| Number of donors homozygous/heterozygous. |
| `NUMT_CAT` | Highest NUMT structural category observed (`A`-`E`). |
| `IS_COMPLEX` | Flag if `NUMT_CAT` is C, D, or E. |
| `IS_REF` | Present if the variant is a Reference NUMT (exists in GRCh38). Absent if Novel. |
| `BLAST_VAL` | Flag indicating the sequence passed BLAST homology validation. |
| `STRUCT_POLY`| Structural polymorphism flag across the population (`NO`/`WITHIN_TIER`/`YES`). |
| `CHRM_RANGE` | The mitochondrial coordinates of the inserted sequence. |

**FORMAT Fields (`GT:CAT:MZ:PZ`):**

| Field | Meaning |
|-------|---------|
| `GT` | Genotype (`0/1` = HET, `1/1` = HOM, `0/0` = REF/Not Detected). |
| `CAT` | Structural category for this specific donor (e.g., `A_single_block`). |
| `MZ` / `PZ` | Maternal and Paternal sub-zygosities. |

### 3. `HPRC_NUMT_callset.bed`

Standard BED6+ format file for quick visualization in UCSC Genome Browser. Contains `hg38_chr`, `start`, `end`, `NUMT_ID`, and basic category information.
