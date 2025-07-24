# Usage

Selscape provides two modes: **testing mode** for learning the workflow with example data, and **main analysis workflow** for analyzing your own datasets.

## Testing Mode

**1. Configure for Testing**

Ensure `testing_mode: true` in your main configuration:
```bash
# config/main.yaml 
testing_mode:  true
```
**2. Run Test Analysis**

The repository includes small example datasets (YRI and CHS populations, chromosomes 21-22) in `examples/data/Human/`, so you can start testing immediately:

```bash
# Local execution 
snakemake --cores 1

# HPC execution 
snakemake --profile config/slurm/ --jobs 10
```

## Main Analysis Workflow
**1. Prepare your Data**
Place your raw VCF files in the correct directory structure:
```
resources/data/{species}/raw/
├── chr1.vcf.gz
├── chr1.vcf.gz.tbi 
├── chr2.vcf.gz 
├── chr2.vcf.gz.tbi
└── ...
```

**2. Configure for Main Workflow**
```bash
# config/main.yaml
testing_mode: false

species: "YourSpecies"

populations:
  - PopA
  - PopB

data:
  # Update paths to match your data structure
  vcf_files:
    base_path: "resources/data"
    file_prefix: ""
    file_suffix: ".vcf.gz"
    chr_name: "chr{i}"  # Adjust to match your naming

  metadata: "resources/data/YourSpecies_metadata.txt"
  
  chromosomes:
    - 1
    - 2
    # ... add all chromosomes
```


**3. Run Main  Analysis**

Run Analysis on HPC:
```bash
snakemake --profile config/slurm/ -j 30
```

Run Analysis locally:
```bash
snakemake -j 5
```
## Clean Restart
To start fresh (removes all previous results):
```bash
rm -rf resources/ results/ logs/
```

## Basic Options

|Command                 |Purpose                                |
|------------------------|---------------------------------------|
|-j 30                   |Use 30 parallel jobs                   |
|--profile config/slurm/ |Submit jobs to SLURM cluster           |
|-n                      |Dry run (preview what will be executed)|


For further information, please visit the [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html#) documentation website.
