# Selscape

  `selscape` is a flexible workflow framework designed to study patterns of natural selection across different species. It integrates different detection methods such as betascan to detect long-term balancing selection, selscan for positive selection and dadi-cli for distribution of fitness effects (DFE). In addition, to understand biological functions of selected candidates, annotation with ANNOVAR and functional enrichment with Gowinda are performed.


## Key Features

`selscape` is designed for flexibility and supports:

| Feature | Description |
| - | - |
| **Multi-species support** | Configurable for any species through YAML files |
| **Selection methods** | betascan, selscan, dadi-cli integration |
| **Flexible data input** | Customizable VCF file naming and chromosome schemes |
| **Functional analysis** | Gene annotation and GO enrichment |
| **Scalable execution** | Snakemake workflow with cluster support |

## Requirements
`selscape` works on Linux operating systems with the following dependencies:
  - bcftools=1.21
  - bedtools=2.31.1
  - bioconductor-biomart=2.62.0
  - cyvcf2=0.31.0
  - dadi=2.3.0
  - dadi-cli=0.9.1
  - matplotlib=3.9.4
  - nlopt=2.7.1
  - numpy=1.26.4
  - pandas=2.2.3
  - pip=25.0.1
  - plink=1.90b6.21
  - pyliftover=0.4.1
  - pysam=0.23.0
  - python=3.11.13
  - r-base=4.4.3
  - r-ggplot2=3.5.2
  - r-qqman=0.1.9
  - r-readr=2.1.5
  - scipy=1.13.1
  - snakemake=7.32.4


## Installation

Users can install `selscape` using the following steps:

1. **Install Mambaforge** (if not already installed):
   [Mambaforge installation guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

2. **Clone the repository**:
```bash
git clone https://github.com/xin-huang/selscape.git
cd selscape
```
3. **Create and activate the environment**:
```bash
mamba env create -f workflow/envs/env.yaml
mmaba activate selscape
```

Users should manually download [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and place it in `resources/tools`. An institutional email adress is required for registration.

## Help

For detailed configuration and usage instructions, see the documentation.
If you need further help, such as reporting a bug or suggesting a feature, please open an issue.
