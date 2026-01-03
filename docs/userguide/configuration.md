# Configuration

Selscape uses multiple YAML configuration files to control different aspects of the analysis:

| File | Purpose |
|------|---------|
| `config/main.yaml` | Core settings (species, populations, data paths) |
| `config/selscan.yaml` | selscan parameters |
| `config/betascan.yaml` | BetaScan parameters |
| `config/scikit-allel.yaml` | scikit-allel parameters |
| `config/dadi-cli.yaml` | dadi-cli parameters |
| `config/slurm/config.yaml` | HPC cluster settings (optional) |

## main.yaml

### Example
```yaml
# Species identification
species: "Human"
tax_id: 9606
ref_genome: "hg38"
ploidy: 2

# Populations to analyze
populations:
  - YRI
  - CHS

# VCF file configuration
data_folder: "examples/data/Human/raw"
vcf_prefix: "full_chr"
vcf_suffix: ".vcf.gz"
# Files should be named: {vcf_prefix}{chromosome}{vcf_suffix}
# Example: full_chr21.vcf.gz

metadata: "examples/data/Human/metadata/example_metadata.txt"

# Chromosomes to analyze
chromosomes:
  - 20
  - 21

# Ancestral alleles for polarization (optional but recommended)
anc_alleles:
  path: "examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh38"
  prefix: "homo_sapiens_ancestor_chr"
  # Files: {prefix}{chromosome}.bed.gz

# Annotation files
genome_annotation: "examples/data/Human/annotation/Human.gtf.gz"
gene2go: "examples/data/Human/annotation/gene2go.gz"

# Quality control
hwe_pvalue: 0.001  # Hardy-Weinberg equilibrium threshold

# Repeat regions (optional)
rmsk: "examples/data/Human/repeats/hg38.rmsk.autosomes.bed"
seg_dup: "examples/data/Human/repeats/hg38.seg.dups.autosomes.bed"
sim_rep: "examples/data/Human/repeats/hg38.simple.repeats.autosomes.bed"

# Method configuration files
selscan_config: "config/selscan.yaml"
betascan_config: "config/betascan.yaml"
dadi_config: "config/dadi-cli.yaml"
scikit_allel_config: "config/scikit-allel.yaml"
```

### Parameters

| Parameter | Type | Description | Required |
|-----------|------|-------------|---------|
| species | string | Species name | Y |
| tax_id | integer | NCBI taxonomy ID | Y |
| ref_genome | string | Reference genome build | Y |
| anc_alleles.path | string |Base folder for ancestral allele information |N |
| anc_alleles.prefix | string |Ancestral allele information filename prefix | N |
| genome_annotation | string | Path to genome annotation GTF file | Y |
| gene2go | string | Path to gene2go mapping file | Y |
| hwe_pvalue | float |Hardy-Weinberg equilibrium p-value threshold | Y |
| hwe_pvalue | float |Hardy-Weinberg equilibrium p-value threshold | Y |
| rmsk | string | Path to repeat masker BED file | N |
| seg_dup | string | Path to segmental duplications BED file | N |
| sim_rep | string | Path to simple repeats BED file | N |
| ploidy | integer | Organism ploidy | Y |
| populations | list | List of population IDs to analyze | Y |
| data_folder | string |Base folder for input VCF files | Y |
| vcf_prefix | string | VCF filename prefix |Y |
| vcf_suffix | string | VCF filename suffix |Y |
| metadata | string |Path to sample metadata file |Y |
| chromosomes | list |List of chromosomes to analyze | Y |
| betascan_config | string | Path to BetaScan configuration file | Y |
| selscan_config | string | Path to selscan configuration file | Y |
| dadi_config | string | Path to dadi-cli configuration file | Y |
| scikit_allel_config | string | Path to scikit-allel configuration file | Y |




## betascan.yaml

### Example
```yaml
# Data type
unfolded: True  # Use polarized data (requires ancestral alleles)

# Core allele frequency parameter
core_frq: 0.15

# Allele frequency filters
min_af: 0.05  # Minimum allele frequency
max_af: 0.95  # Maximum allele frequency

# Statistical thresholds
top_proportion: 0.0005  # Top proportion for outlier identification

# Manhattan plot settings
manhattan_plot_width: 640
manhattan_plot_height: 240
manhattan_plot_color1: "#56B4E9"
manhattan_plot_color2: "#F0E442"
```

### Parameters

| Parameter | Type | Description | Required |
|-----------|------|-------------|---------|
| unfolded | boolean | Use polarized/unfolded allele frequency spectrum | Y |
| core_frq | float | Core allele frequency parameter | Y |
| min_af | float | Minimum allele frequency | Y |
| max_af | float | Maximum allele frequency | Y |
| top_proportion | float | Top proportion for outlier identification | Y |


## selscan.yaml

### Example
```yaml
# Data type
unphased: true  # Set to false if you have phased data

# Within-population statistics
wp_stats:
  - ihs   
  - nsl    

# Cross-population statistics  
xp_stats:
  - xpehh  
  - xpnsl

# Minor allele frequency threshold
maf: 0.05

# Statistical thresholds
top_proportion: 0.0005

# Manhattan plot settings
manhattan_plot_width: 640
manhattan_plot_height: 240
manhattan_plot_color1: "#56B4E9"
manhattan_plot_color2: "#F0E442"
```

### Parameters

| Parameter | Type | Description | Required |
|-----------|------|-------------|---------|
| unphased | boolean |Use unphased (true) or phased (false) data | Y |
| wp_stats | list | Within-population statistics to compute | Y |
| xp_stats | list | Cross-population statistics to compute | Y |
| maf | float | Minor allele frequency threshold | Y |
| top_proportion | float | Top proportion for outlier identification | Y |

## scikit-allel.yaml

### Example
```yaml
# Within-population statistics
wp_stats:
  - windowed_tajima_d   
  - moving_tajima_d  

# Windowed approach (base pair windows)
windowed_window_sizes: [100_000]
windowed_step_size_ratios: [1]  # 1 = non-overlapping windows

# Moving approach (SNP windows)
moving_window_sizes: [100]
moving_step_size_ratios: [1]

# Cross-population statistics 
xp_stats:
  - delta_moving_tajima_d

# XP moving approach (SNP windows)
delta_moving_window_sizes: [100]
delta_moving_step_size_ratios: [1]

# Statistical thresholds
top_proportion: 0.05

# Manhattan plot settings
manhattan_plot_width: 640
manhattan_plot_height: 240
manhattan_plot_color1: "#56B4E9"
manhattan_plot_color2: "#F0E442"
```

### Parameters

| Parameter | Type | Description | Required |
|-----------|------|-------------|---------|
| wp_stats | list | Within-population statistics to compute | Y |
| xp_stats | list | Cross-population statistics to compute  | Y |
| wtjd_window_sizes | list | Window sizes in base pairs for windowed tajima's d (wtjd) statistic | Y |
| wtjd_step_size_ratios | list | Step size as fraction of window in base pair for wtjd statistic| Y |
| mtjd_window_sizes | list | Window sizes in SNPs for moving tajima's d (mtjd) statistic| Y |
| mtjd_step_size_ratios | list | Step size as fraction of window in SNPs for mtjd statistic | Y |
| dtjd_moving_window_sizes | list | Window sizes in SNPs for xp moving delta tajima's d (dtjd) statistic | Y |
| dtjd_moving_step_size_ratios | list | Step size as fraction of window in SNPs for cross population dtjd statistic| Y |
| top_proportion | float | Top proportion for outlier identification | Y |

## dadi-cli.yaml

### Example
```yaml
# Data type
unfolded: True

# Demographic models
demog_1d: two_epoch
demog_1d_p0: "0.5 5 0.5"
demog_1d_ub: "10 10 1"
demog_1d_lb: "10e-5 10e-5 0"
 
# DFE models
dfe_1d: lognormal
dfe_1d_p0: "5 5 0.5"
dfe_1d_ub: "100 1000 1"
dfe_1d_lb: "10e-5 10e-5 0"

# Grid sizes
dm_grid_size: "300 400 500"
dfe_grid_size: "300 400 500"

# Gamma points for generating caches
gamma_pts: 100

# Ratio of non-synonymous/synonymous mutation rates
ratio: 2.31
  
# Number of optimization runs
optimizations: 50
  
# Bootstrap parameters
bootstrap_replicates: 100
chunk_size: 1000000

# Frequency spectrum flags
mask_singletons: false
```

### Parameters

| Parameter | Type | Description | Required |
|-----------|------|-------------|---------|
| unfolded | boolean | Use polarized/unfolded allele frequency spectrum | Y |
| demog_1d | string | 1D demographic model name | Y |
| demog_1d_p0 | string | Initial parameters for 1D demographic model | Y |
| demog_1d_ub | string | Upper bounds for 1D demographic model parameters | Y |
| demog_1d_lb | string | Lower bounds for 1D demographic model parameters | Y |
| dfe_1d | string | 1D DFE model name | Y |
| dfe_1d_p0 | string | Initial parameters for 1D DFE model | Y |
| dfe_1d_ub | string | Upper bounds for 1D DFE model parameters  | Y |
| dfe_1d_lb | string | Lower bounds for 1D DFE model parameters | Y |
| dm_grid_size | string | Grid sizes for demographic inference | Y |
| dfe_grid_size | string | Grid for DFE inference | Y |
| gamma_pts | integer | Number of gamma points for cache generation | Y |
| ratio | float | Ratio of non-synonymous to synonymous mutation rates | Y |
| optimizations | integer | Number of optimization runs | Y |
| bootstrap_replicates | integer | Number of bootstrap replicates | Y |
| chunk_size | integer | Chunk size for bootstrapping | Y |
| mask_singletons | boolean | Whether to mask singleton variants | Y |


## Important Notes

1. **Ancestral Alleles**: If `anc_alleles` is not provided, polarization-dependent analyses (selscan, betascan with unfolded=True, dadi with unfolded=True) will not run.

2. **VCF File Naming**: VCF files must follow the pattern: `{data_folder}/{vcf_prefix}{chromosome}{vcf_suffix}`

    - Example with `vcf_prefix: "full_chr"`: `examples/data/Human/raw/full_chr21.vcf.gz`
    - Example with `vcf_prefix: "chr"`: `resources/data/chr21.vcf.gz`
    - Example with `vcf_prefix: ""`: `resources/data/21.vcf.gz`

3. **Metadata Format**: Tab-separated file with columns: `Sample` and `Population`

4. **Repeat Masking**: Optional but recommended for balancing selection analysis. Set to empty string or null to disable.

5. **Step Size Ratios**: Values from 0 to 1 represent the fraction of window size to step.

    - Example: window size = 100, step size ratio = 0.1 → step size = 10

Please check [Snakemake Configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuration) documentation for additional information.
