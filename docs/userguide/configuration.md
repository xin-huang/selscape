# Configuration

## Essential Setup (main.yaml)

Most users only need to modify these core settings:

## Species and Populations (Human Example)
```yaml
# Main configuration for species and data paths

testing_mode:  true  # Set to false for production analysis

# Species identification
species: "Human"

# Populations to analyze
populations:
  - YRI
  - CHS

# Data location and file naming
data:
  # Start with a small subset for testing 
  chromosomes:
     - 21
     - 22
  
  # Tell selscape how your VCF files are named
  vcf_files:
    base_path: "resources/data"
    file_prefix: ""
    file_suffix: ".vcf.gz"
    chr_name: "{ppl}.chr{i}"  # {ppl} = population, {i} = chromosome

  # Path to your sample metadata file
  metadata: "resources/data/1KG_metadata.txt"

# Reference genome
reference:
  species_code: "hg38"
  human_code: "hg38"
```

## Method Parameters

**Betascan (Balancing Selection)**
```yaml

core_frequencies:
  - 0.15  # Most sensitive for detecting balancing selection

frequency_filters:
  min_af: 0.05  # Remove rare variants
  max_af: 0.95  # Remove nearly fixed variants

folding:
  use_folded: true  # Use if you don't have ancestral states
  ```
**Selscan (Positive Selection)**
```yaml
data_type:
  unphased: true  # Set to false if you have phased data

statistics:
  within_population: [ihs, nsl]
  cross_population: [xpehh, xpnsl]

frequency_filters:
  maf_thresholds: [0.05]  # Minor allele frequency cutoff
  ```

## Advanced Configuration
**main.yaml Parameters**

|       Parameter         |Type                        |Description                     | Default        | 
|----------------|-------------------------------|-----------------------------|-----------------------------| 
|testing_mode|boolean|Use example data for testing|false|
|species|string         |Species identifier for outputs           | "Human"| 
|populations         |list           |Population code to analyze          | ["YRI", "CHS"]| 
|data.chromosomes        |list|Chromosomes to include | [21, 22]  |
|data.vcf_files.base_path|string         |VCF file directory          | "resources/data"| 
|data.vcf_files.chr_name        |string         |Chromosome naming pattern          | "{ppl}.chr{i}"| 
|reference.species_code      |string|Reference genome build | "hg38"  |
|data.metadata|string|Sample metadata file path|"resources/data/1KG_metadata.txt"|
|data.vcf_files.file_prefix|string|VCF filename prefix|""|
|data.vcf_files.file_suffix|string|VCF filename suffix|".vcf.gz"|
|data_sources|dict|Data source configurations|see config|
|reference.annotation_reference|string|Reference genome for annotation databases|"hg38"|

**betascan.yaml Parameters**
|       Parameter         |Type                        |Description                     | Default        | 
|----------------|-------------------------------|-----------------------------|-----------------------------| 
|core_frequencies|list           |Core allele frequencies for B1           | [0.15]| 
|frequency_filters.min_af         |float          |Minimum allele frequency          |0.05 | 
|frequency_filters.max_af          |float|Maximum allele frequency|0.95 | 
|folding.use_folded          |boolean|Use folded SFS|true | 
|thresholds.top_percent         |float|Top % for candidates|0.0005 | 
|quality_control.hwe_pvalue|float|Hardy-Weinberg p-value threshold|0.001|
|quality_control.mask_repeats|boolean|Exclude repetitive regions|true|

**selscan.yaml Parameters**
|       Parameter         |Type                        |Description                     | Default        | 
|----------------|-------------------------------|-----------------------------|-----------------------------| 
|data_type.unphased|boolean            |Handle unphased data           | true| 
|statistics.within_population         |list    |      Within-pop statistics |[ihs, nsl] | 
|statistics.cross_population          |list|Cross-pop statistics|[xpehh, xpnsl] | 
|frequency_filters.maf_thresholds          |list|MAF cutoffs|[0.05] |

**dadi-cli.yaml Parameters**
|       Parameter         |Type                        |Description                     | Default        | 
|----------------|-------------------------------|-----------------------------|-----------------------------| 
|flags.ratio|float           |Transition/transversion ratio          | 2.31| 
|flags.optimizations        |integer|Number of optimization runs           | 50| 
|grid_size.inference         |string|Grid size for inference| "300 400 500"| 
|gamma_pts       |integer|Gamma integration points| 2000| 
|grid_size.cache|string|Grid size for caching|"800 1000 1200"|
|flags.bootstrap_replicates|integer|Bootstrap replicates|100|
|flags.chunk_size|integer|Chunk size for bootstrap|1000000|

**annovar.yaml Parameters**
|Parameter|Type|Description|Default|
|----------------|-------------------------------|-----------------------------|-----------------------------| 
|species.{species}.genome_build|string|Reference genome build|"hg38"|
|settings.threads|integer|Annotation threads|8|



Please check [Snakemake Configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#configuration) documentation for additional information.
