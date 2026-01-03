# Usage

`selscape` provides example data to test the workflow, and supports analyzing your own datasets with flexible configuration.

## Quick Start with Example Data

### 1. Download Example Data

`selscape` includes a workflow to download example data from the 1000 Genomes Project:
```bash
# Download example data (YRI and CHS populations, chromosomes 20-21)
snakemake -s examples/get_example_data.smk -c 1
```

This creates the following structure:
```
examples/data/Human/
├── raw/
│   ├── full_chr20.vcf.gz
│   ├── full_chr20.vcf.gz.tbi
│   ├── full_chr21.vcf.gz
│   └── full_chr21.vcf.gz.tbi
├── ancestral_alleles/
│   └── homo_sapiens_ancestor_GRCh38/
├── repeats/
│   ├── hg38.rmsk.autosomes.bed
│   ├── hg38.seg.dups.autosomes.bed
│   └── hg38.simple.repeats.autosomes.bed
├── annotation/
│   ├── Human.gtf.gz
│   └── gene2go.gz
└── metadata/
    └── example_metadata.txt
```

### 2. Run Analysis on Example Data
```bash
# Local execution (single core)
snakemake -c 1 --configfile config/main.yaml

# Local execution (multiple cores)
snakemake -c 4 --configfile config/main.yaml

# HPC cluster (SLURM)
snakemake -c 1 --configfile config/main.yaml --profile config/slurm/ -j 50
```

### 3. Generate HTML Report

Snakemake can automatically generate an interactive HTML report with main results, plots, and tables:
```bash
# Generate report after workflow completion
snakemake -c 1 --configfile config/main.yaml --report

# Customize report location and file name
snakemake -c 1 --configfile config/main.yaml --report /path/to/my_report.html
```

The report for the analysis of the example data can be found [here](report.html).

## Analyzing Your Own Data

### 1. Prepare Your Data

#### VCF Files

Organize your VCF files according to your chosen naming scheme. The configuration will tell Selscape how to find them.

Example structures:
```bash
# Option 1: Chromosome in filename
resources/data/
├── chr1.vcf.gz
├── chr1.vcf.gz.tbi
├── chr2.vcf.gz
└── chr2.vcf.gz.tbi

# Option 2: Prefix + chromosome
resources/data/
├── Sample_chr1.vcf.gz
├── Sample_chr1.vcf.gz.tbi
├── Sample_chr2.vcf.gz
└── Sample_chr2.vcf.gz.tbi
```

**Requirements:**

- All VCF files must be bgzipped (`.vcf.gz`)
- All VCF files must be indexed (`.vcf.gz.tbi`)
- Biallelic SNPs only (can be filtered by selscape)

#### Metadata File

Create a tab-separated file mapping samples to populations:
```
Sample	Population
NA19001	YRI
NA19002	YRI
NA19003	YRI
HG00419	CHS
HG00420	CHS
HG00421	CHS
```

**Requirements:**

- Header line: `Sample` and `Population`
- Tab-separated
- Sample IDs must match VCF file samples

#### Ancestral Alleles (Optional)

If available, ancestral alleles enable polarized analyses (selscan, unfolded BetaScan/dadi):
```bash
ancestral_alleles/
├── anc_chr1.bed.gz
├── anc_chr1.bed.gz.tbi
├── anc_chr2.bed.gz
└── anc_chr2.bed.gz.tbi
```

**Format:** BED file with ancestral base in 4th column:
```
chr1	0	1	A
chr1	1	2	C
chr1	2	3	G
```

### 2. Configure Analysis

Edit `config/main.yaml`:
```yaml
# Species identification
species: "YourSpecies"
tax_id: 9606  # Update for your species
ref_genome: "your_ref"
ploidy: 2

# Populations to analyze
populations:
  - Pop1
  - Pop2

# VCF file configuration - ADJUST TO YOUR NAMING SCHEME
data_folder: "resources/data"
vcf_prefix: "chr"  # or "Sample_chr", "full_chr", etc.
vcf_suffix: ".vcf.gz"

# Sample metadata
metadata: "resources/metadata/samples.txt"

# Chromosomes to analyze
chromosomes:
  - 1
  - 2
  - 3
  # ... add all chromosomes

# Ancestral alleles (optional)
anc_alleles:
  path: "resources/ancestral_alleles"
  prefix: "anc_chr"

# Annotation files
genome_annotation: "resources/annotation/genome.gtf.gz"
gene2go: "resources/annotation/gene2go.gz"

# Repeat regions (optional but recommended)
rmsk: "resources/repeats/repeats.bed"
seg_dup: "resources/repeats/seg_dups.bed"
sim_rep: "resources/repeats/simple_repeats.bed"

# Analysis tool configuration
selscan_config: "config/selscan.yaml"
betascan_config: "config/betascan.yaml"
dadi_config: "config/dadi-cli.yaml"
scikit_allel_config: "config/scikit-allel.yaml"
```

### 3. Run Analysis
```bash
# Preview what will be executed (dry run)
snakemake -np --configfile /path/to/your/main.yaml

# Local execution
snakemake -c 1 --configfile /path/to/your/main.yaml

# HPC cluster (SLURM)
snakemake -c 1 --configfile /path/to/your/main.yaml --profile config/slurm/ --j 50
```

## Running Specific Analyses

You can run only specific methods by requesting specific output files:
```bash
# Only BetaScan with B1
snakemake -c 1 --configfile /path/to/your/main.yaml -R --until plot_gowinda_enrichment_betascan

# Only selscan with single population statistics (requires ancestral alleles)
snakemake -c 1 --configfile /path/to/your/main.yaml -R --until plot_gowinda_enrichment_selscan

# Only scikit-allel with Tajima's D
snakemake -c 1 --configfile /path/to/your/main.yaml -R --until plot_gowinda_enrichment_tajima_d

# Only dadi-cli for DFE inference
snakemake -c 1 --configfile /path/to/your/main.yaml -R --until plot_fitted_dfe
```

## Clean Restart

To remove all results and start fresh:
```bash
# Remove all results (keeps raw data)
rm -rf results/ logs/

# Remove everything including downloaded data
rm -rf resources/ results/ logs/
```

## Workflow Options

| Command | Purpose |
|---------|---------|
| `-c 1` or `--cores 1` | Use 1 core (safest for testing) |
| `-c 4` | Use 4 parallel jobs (local) |
| `-j or --jobs 50` | Submit up to 50 jobs (cluster) |
| `--profile config/slurm/` | Use SLURM cluster configuration |
| `-n` or `--dry-run` | Preview what will be executed |
| `-p` | Print shell commands |
| `--forceall` | Re-run all steps |
| `--unlock` | Unlock working directory after crash |

## Troubleshooting

**Workflow locked:**
```bash
snakemake --unlock --configfile /path/to/your/main.yaml
```

**Re-run specific rule:**
```bash
snakemake --forcerun <rule_name> -c 1 --configfile /path/to/your/main.yaml
```

**Check what will run:**
```bash
snakemake -np --configfile /path/to/your/main.yaml
```

For more information, visit the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html).
