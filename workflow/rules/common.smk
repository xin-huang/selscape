# Copyright 2025 Xin Huang and Simon Chen
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html

import numpy as np
import yaml
import os
import sys
from itertools import combinations
from snakemake.utils import validate

# CONFIGURATION LOADING

main_config = config

# Load configs

with open(main_config["selscan_config"], "r") as f:
    selscan_config = yaml.safe_load(f)

with open(main_config["betascan_config"], "r") as f:
    betascan_config = yaml.safe_load(f)

with open(main_config["dadi_config"], "r") as f:
    dadi_config = yaml.safe_load(f)

with open(main_config["scikit_allel_config"], "r") as f:
    scikit_allel_config = yaml.safe_load(f)

# Config Validation
validate(main_config, schema="../schemas/config.schema.yaml")
validate(selscan_config, schema="../schemas/selscan.schema.yaml")
validate(betascan_config, schema="../schemas/betascan.schema.yaml")
validate(dadi_config, schema="../schemas/dadi-cli.schema.yaml")
validate(scikit_allel_config, schema="../schemas/scikit-allel.schema.yaml")


# Load and validate all dataset configs
dataset_configs = {}
for cfg_path in main_config["datasets"]:
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    validate(cfg, schema="../schemas/dataset.schema.yaml")
    dataset_configs[cfg["dataset"]] = cfg


# Flat lists of (dataset, species, population) tuples — used for 1pop rules
DATASET_1POP = [
    (ds, cfg["species"], pop, cfg["ref_genome"])
    for ds, cfg in dataset_configs.items()
    for pop in cfg["populations"]
]

# Flat lists of (dataset, species, pair) tuples — used for 2pop/xp rules
DATASET_2POP = [
    (ds, cfg["species"], "_".join(pair), cfg["ref_genome"])
    for ds, cfg in dataset_configs.items()
    for pair in combinations(cfg["populations"], 2)
]

# Ancestral allele gating policy: all datasets must provide anc_alleles, otherwise disable anc-only globally.
datasets_without_anc = [
    ds for ds, cfg in dataset_configs.items()
    if not cfg.get("anc_alleles")
]
if datasets_without_anc:
    print(
        "Missing anc_alleles for dataset(s): "
        + ", ".join(datasets_without_anc)
        + ". Running workflow in no-ancestral-alleles mode (anc-only rules disabled for all datasets).",
        file=sys.stderr,
    )
    print(f"Circos plots for positive selection will not be generated for: {', '.join(datasets_without_anc)}. To generate them, add anc_alleles to the dataset config.", file=sys.stderr)
    datasets_with_anc = []
else:
    datasets_with_anc = list(dataset_configs.keys())


DATASET_1POP_ANC = [(ds, sp, pop, rg) for ds, sp, pop, rg in DATASET_1POP if ds in datasets_with_anc]
DATASET_2POP_ANC = [(ds, sp, pair, rg) for ds, sp, pair, rg in DATASET_2POP if ds in datasets_with_anc]


all_chromosomes = sorted(set(
    c for cfg in dataset_configs.values()
    for c in cfg["chromosomes"]
), key=str)


def expand_1pop(pattern, anc_only=False, **kw):
    source = DATASET_1POP_ANC if anc_only else DATASET_1POP
    return [
        f
        for ds, sp, pop, rg in source
        for f in expand(pattern, dataset=ds, species=sp, ppl=pop, ref_genome=rg, **kw)
    ]


def expand_2pop(pattern, anc_only=False, **kw):
    source = DATASET_2POP_ANC if anc_only else DATASET_2POP
    return [
        f
        for ds, sp, pair, rg in source
        for f in expand(pattern, dataset=ds, species=sp, pair=pair, ref_genome=rg, **kw)
    ]


SELSCAN_KW = dict(
    method=selscan_config["wp_stats"],
    maf=selscan_config["maf"],
    cutoff=selscan_config["top_proportion"],
)

SELSCAN_XP_KW = dict(
    method=selscan_config["xp_stats"],
    maf=selscan_config["maf"],
    cutoff=selscan_config["top_proportion"],
)

BETASCAN_KW = dict(
    core_frq=betascan_config["core_frq"],
    cutoff=betascan_config["top_proportion"],
)

DADI_1D_KW = dict(
    demog=dadi_config["demog_1d"],
)

TAJIMAD_MOVING_KW = dict(
    method="moving_tajima_d",
    window=scikit_allel_config["mtjd_window_sizes"],
    step=scikit_allel_config["mtjd_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

TAJIMAD_WINDOWED_KW = dict(
    method="windowed_tajima_d",
    window=scikit_allel_config["wtjd_window_sizes"],
    step=scikit_allel_config["wtjd_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

DELTA_TAJIMAD_KW = dict(
    method="delta_tajima_d",
    window=scikit_allel_config["dtjd_window_sizes"],
    step=scikit_allel_config["dtjd_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

WP_SELSCAN_STATS = selscan_config.get("wp_stats") or []
_MTJD_ON = bool(TAJIMAD_MOVING_KW["window"] and TAJIMAD_MOVING_KW["step"])
_WTJD_ON = bool(TAJIMAD_WINDOWED_KW["window"] and TAJIMAD_WINDOWED_KW["step"])

_POS_CIRCOS_ALL_TRACKS = [
    {"stat": "ihs", "enabled": "ihs" in WP_SELSCAN_STATS,
     "name": "iHS", "file": "ihs_scores", "score_col": "normalized_ihs", "color": "#1f77b4"},
    {"stat": "nsl", "enabled": "nsl" in WP_SELSCAN_STATS,
     "name": "nSL", "file": "nsl_scores", "score_col": "normalized_nsl", "color": "#ff7f0e"},
    {"stat": "mtjd", "enabled": _MTJD_ON,
     "name": "mtjd", "file": "mtjd_scores", "score_col": "tajima_d", "color": "#2ca02c"},
    {"stat": "wtjd", "enabled": _WTJD_ON,
     "name": "wtjd", "file": "wtjd_scores", "score_col": "tajima_d", "color": "#d62728"},
]
_POS_CIRCOS_SLOTS = [[65, 75], [50, 60], [35, 45], [20, 30]]

POS_CIRCOS_TRACKS = [
    {k: v for k, v in t.items() if k not in ("stat", "enabled")} | {"r_range": slot}
    for t, slot in zip([t for t in _POS_CIRCOS_ALL_TRACKS if t["enabled"]], _POS_CIRCOS_SLOTS)
]


_BAL_CIRCOS_ALL_TRACKS = [
    {"stat": "b1", "enabled": True,
     "name": "B1", "file": "b1_scores", "score_col": "B1", "color": "#1f77b4"},
    {"stat": "mtjd", "enabled": _MTJD_ON,
     "name": "mtjd", "file": "mtjd_bal_scores", "score_col": "tajima_d", "color": "#2ca02c"},
    {"stat": "wtjd", "enabled": _WTJD_ON,
     "name": "wtjd", "file": "wtjd_bal_scores", "score_col": "tajima_d", "color": "#d62728"},
]
_BAL_CIRCOS_SLOTS = [[60, 75], [40, 55], [20, 35]]

BAL_CIRCOS_TRACKS = [
    {k: v for k, v in t.items() if k not in ("stat", "enabled")} | {"r_range": slot}
    for t, slot in zip([t for t in _BAL_CIRCOS_ALL_TRACKS if t["enabled"]], _BAL_CIRCOS_SLOTS)
]


def get_pos_circos_inputs(wc):
    base = "results/positive_selection"
    enabled = {t["file"] for t in POS_CIRCOS_TRACKS}
    inputs = {}
    if "ihs_scores" in enabled:
        inputs["ihs_scores"] = f"{base}/selscan/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/ihs_{SELSCAN_KW['maf']}/{wc.ppl}.normalized.ihs.scores"
    if "nsl_scores" in enabled:
        inputs["nsl_scores"] = f"{base}/selscan/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/nsl_{SELSCAN_KW['maf']}/{wc.ppl}.normalized.nsl.scores"
    if "mtjd_scores" in enabled:
        inputs["mtjd_scores"] = f"{base}/scikit-allel/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/{TAJIMAD_MOVING_KW['method']}/{TAJIMAD_MOVING_KW['window'][0]}_{TAJIMAD_MOVING_KW['step'][0]}/{wc.ppl}.{TAJIMAD_MOVING_KW['method']}.scores"
    if "wtjd_scores" in enabled:
        inputs["wtjd_scores"] = f"{base}/scikit-allel/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/{TAJIMAD_WINDOWED_KW['method']}/{TAJIMAD_WINDOWED_KW['window'][0]}_{TAJIMAD_WINDOWED_KW['step'][0]}/{wc.ppl}.{TAJIMAD_WINDOWED_KW['method']}.scores"
    return inputs


def get_bal_circos_inputs(wc):
    base = "results/balancing_selection"
    enabled = {t["file"] for t in BAL_CIRCOS_TRACKS}
    inputs = {}
    if "b1_scores" in enabled:
        inputs["b1_scores"] = f"{base}/betascan/{wc.species}/{wc.dataset}/{wc.ppl}/m_{BETASCAN_KW['core_frq']}/{wc.ppl}.{get_ref_genome(wc)}.m_{BETASCAN_KW['core_frq']}.b1.scores"
    if "mtjd_bal_scores" in enabled:
        inputs["mtjd_bal_scores"] = f"{base}/scikit-allel/{wc.species}/{wc.dataset}/{TAJIMAD_MOVING_KW['method']}/{wc.ppl}/{TAJIMAD_MOVING_KW['window'][0]}_{TAJIMAD_MOVING_KW['step'][0]}/{wc.ppl}.{TAJIMAD_MOVING_KW['method']}.merged.scores"
    if "wtjd_bal_scores" in enabled:
        inputs["wtjd_bal_scores"] = f"{base}/scikit-allel/{wc.species}/{wc.dataset}/{TAJIMAD_WINDOWED_KW['method']}/{wc.ppl}/{TAJIMAD_WINDOWED_KW['window'][0]}_{TAJIMAD_WINDOWED_KW['step'][0]}/{wc.ppl}.{TAJIMAD_WINDOWED_KW['method']}.merged.scores"
    return inputs


XP_SELSCAN_STATS = selscan_config.get("xp_stats") or []
XP_ALLEL_STATS = scikit_allel_config.get("xp_stats") or []

_XP_CIRCOS_ALL_TRACKS = [
    {"stat": "xpehh", "enabled": "xpehh" in XP_SELSCAN_STATS,
     "name": "XP-EHH", "file": "xpehh_scores", "score_col": "normalized_xpehh", "color": "#1f77b4"},
    {"stat": "xpnsl", "enabled": "xpnsl" in XP_SELSCAN_STATS,
     "name": "XP-nSL", "file": "xpnsl_scores", "score_col": "normalized_xpnsl", "color": "#ff7f0e"},
    {"stat": "dtjd", "enabled": bool(XP_ALLEL_STATS and DELTA_TAJIMAD_KW["window"] and DELTA_TAJIMAD_KW["step"]),
     "name": "dtjd", "file": "dtjd_scores", "score_col": "delta_tajima_d", "color": "#2ca02c"},
]
_XP_CIRCOS_SLOTS = [[60, 75], [40, 55], [20, 35]]

XP_CIRCOS_TRACKS = [
    {k: v for k, v in t.items() if k not in ("stat", "enabled")} | {"r_range": slot}
    for t, slot in zip([t for t in _XP_CIRCOS_ALL_TRACKS if t["enabled"]], _XP_CIRCOS_SLOTS)
]


def get_xp_circos_inputs(wc):
    base = "results/positive_selection"
    enabled = {t["file"] for t in XP_CIRCOS_TRACKS}
    inputs = {}
    if "xpehh_scores" in enabled:
        inputs["xpehh_scores"] = f"{base}/selscan/{wc.species}/{wc.dataset}/2pop/{wc.pair}/xpehh_{SELSCAN_XP_KW['maf']}/{wc.pair}.normalized.xpehh.scores"
    if "xpnsl_scores" in enabled:
        inputs["xpnsl_scores"] = f"{base}/selscan/{wc.species}/{wc.dataset}/2pop/{wc.pair}/xpnsl_{SELSCAN_XP_KW['maf']}/{wc.pair}.normalized.xpnsl.scores"
    if "dtjd_scores" in enabled:
        inputs["dtjd_scores"] = f"{base}/scikit-allel/{wc.species}/{wc.dataset}/2pop/{wc.pair}/{DELTA_TAJIMAD_KW['method']}/{DELTA_TAJIMAD_KW['window'][0]}_{DELTA_TAJIMAD_KW['step'][0]}/{wc.pair}.{DELTA_TAJIMAD_KW['method']}.merged.scores"
    return inputs


selscan_method_names = {
    "ihs": "iHS",
    "nsl": "nSL",
    "xpehh": "XP-EHH",
    "xpnsl": "XP-nSL",
}

# HELPER FUNCTIONS


def get_dataset_cfg(wildcards):
    """Return the config dict for the given dataset wildcard."""
    return dataset_configs[wildcards.dataset]


def get_anc_allele_bed(wildcards):
    """Get ancestral allele bed files."""
    cfg = get_dataset_cfg(wildcards)
    anc = cfg["anc_alleles"]
    chr_prefix = anc.get("chr_prefix", "")
    return f"{anc['path']}/{anc['prefix']}.{chr_prefix}{wildcards.i}.bed.gz"


def get_vcf_input_path(wildcards):
    """Get VCF input file path."""
    cfg = get_dataset_cfg(wildcards)
    return f"{cfg['data_folder']}/{cfg['vcf_prefix']}{wildcards.i}{cfg['vcf_suffix']}"


def get_metadata(wildcards):
    """Get sample metadata file path for the given dataset."""
    return get_dataset_cfg(wildcards)["metadata"]


def get_genome_annotation(wildcards):
    """Get genome annotation GTF path for the given dataset."""
    return get_dataset_cfg(wildcards)["genome_annotation"]


def get_gene2go(wildcards):
    """Get gene2go mapping file path for the given dataset."""
    return get_dataset_cfg(wildcards)["gene2go"]


def get_rmsk(wildcards):
    """Get repeatmasker BED path for the given dataset (empty string if null)."""
    return get_dataset_cfg(wildcards)["rmsk"] or ""


def get_seg_dup(wildcards):
    """Get segmental duplication BED path for the given dataset (empty string if null)."""
    return get_dataset_cfg(wildcards)["seg_dup"] or ""


def get_sim_rep(wildcards):
    """Get simple repeats BED path for the given dataset (empty string if null)."""
    return get_dataset_cfg(wildcards)["sim_rep"] or ""


def get_hwe_pvalue(wildcards):
    """Get HWE p-value threshold for the given dataset."""
    return get_dataset_cfg(wildcards)["hwe_pvalue"]


def get_ploidy(wildcards):
    """Get ploidy for the given dataset."""
    return get_dataset_cfg(wildcards)["ploidy"]


def get_tax_id(wildcards):
    """Get NCBI taxonomy ID for the given dataset."""
    return get_dataset_cfg(wildcards)["tax_id"]


def get_chromosomes(wildcards):
    """Get chromosome list for the given dataset."""
    return get_dataset_cfg(wildcards)["chromosomes"]


def _top_pct(wildcards) -> str:
    """Return 'Top X%' string based on wildcards.cutoff (float or str)."""
    return f"Top {float(wildcards.cutoff)*100:g}%"


def _vs_pair(wildcards) -> str:
    """Return 'A vs B' string from wildcards.pair formatted as 'A_B'."""
    return " vs ".join(wildcards.pair.split("_"))


def format_method_name(method):
    """Format scikit-allel method name for display in titles."""
    formatted = method.replace("tajima_d", "Tajima's D").replace("_", " ")
    return formatted.title().replace("Tajima'S D", "Tajima's D")


def get_phasing_flag(wildcards):
    """Return --unphased flag if dataset is unphased or polyploid."""
    ploidy = get_dataset_cfg(wildcards)["ploidy"]
    if selscan_config["unphased"] or ploidy > 2:
        return "--unphased"
    return ""


def get_ref_genome(wildcards):
    """Get reference genome build for the given dataset."""
    return get_dataset_cfg(wildcards)["ref_genome"]


def selscan_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for within-population selscan Manhattan plot."""
    return {
        "Dataset": wildcards.dataset,
        "Population": wildcards.ppl,
        "Minor Allele Frequency": wildcards.maf,
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }

def selscan_xp_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for cross-population selscan Manhattan plot."""
    return {
        "Dataset": wildcards.dataset,
        "Populations": _vs_pair(wildcards),
        "Minor Allele Frequency": wildcards.maf,
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }

def betascan_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for betascan plots (Manhattan or Enrichment), includes core frequency."""
    return {
        "Dataset": wildcards.dataset,
        "Population": wildcards.ppl,
        "Core Frequency": str(wildcards.core_frq),
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }

def tajima_d_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for Tajima's D plots (both windowed and moving)."""
    method_name = (
        "Moving Tajima's D"
        if wildcards.method.startswith("moving")
        else "Windowed Tajima's D"
    )
    window_unit = " SNPs" if wildcards.method.startswith("moving") else " bp"
    step_size = int(float(wildcards.step) * int(wildcards.window))
    return {
        "Dataset": wildcards.dataset,
        "Population": wildcards.ppl,
        "Window": f"{wildcards.window}{window_unit}",
        "Step": f"{step_size}{window_unit}",
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }

def delta_tajima_d_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for delta Tajima's D plots (cross-population)."""
    step_size = int(float(wildcards.step) * int(wildcards.window))
    return {
        "Dataset": wildcards.dataset,
        "Populations": _vs_pair(wildcards),
        "Window": f"{wildcards.window} SNPs",
        "Step": f"{step_size} SNPs",
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }

def get_betascan_vcf_dir(wildcards):
    """Return polarized_data or processed_data depending on anc_alleles config."""
    if cfg.get("anc_alleles") and betascan_config["unfolded"]:
        return "polarized_data"
    return "processed_data"


def get_folding_flag(wildcards):
    """Return -fold flag if dataset has no ancestral alleles or betascan is folded."""
    if wildcards.dataset in datasets_with_anc and betascan_config["unfolded"]:
        return "-fold"
    return ""


def get_dadi_vcf_dir(wildcards):
    """Return polarized_data or processed_data depending on dadi unfolded config."""
    if wildcards.dataset in datasets_with_anc and dadi_config["unfolded"]:
        return "polarized_data"
    return "processed_data"


def get_polarization_flag(wildcards):
    """Return --polarized flag if dataset has ancestral alleles and dadi is unfolded."""
    if wildcards.dataset in datasets_with_anc and dadi_config["unfolded"]:
        return "--polarized"
    return ""

def fitted_1pop_dm_labels(wildcards, type: str = "Model Fit") -> dict[str, str]:
    """Labels for 1-population demographic model fit plots."""
    return {
        "Dataset": wildcards.dataset,
        "Population": wildcards.ppl,
        "Demographic Model": wildcards.demog,
        "Type": type,
    }

def fitted_dfe_labels(wildcards, type: str = "Model Fit") -> dict[str, str]:
    """Labels for 1-population DFE model fit plots."""
    return {
        "Dataset": wildcards.dataset,
        "Population": wildcards.ppl,
        "Demographic Model": wildcards.demog,
        "DFE Model": wildcards.dfe,
        "Type": type,
    }

def add_selscan_title(wildcards, input):
    """Generate title for selscan plots and tables."""
    cutoff_pct = float(wildcards.cutoff) * 100
    pop_id = wildcards.get("ppl") or wildcards.get("pair")

    if hasattr(input, "scores"):
        return " ".join([
            f"{pop_id}",
            f"(MAF={wildcards.maf},",
            f"Top {cutoff_pct:.2f}%)",
        ])

    return " ".join([
        f"{pop_id}",
        selscan_method_names[wildcards.method],
        f"(MAF={wildcards.maf},",
        f"Top {cutoff_pct:.2f}%)",
    ])


def add_scikit_allel_title(wildcards, input=None):
    """Generate title for scikit-allel plots and tables."""
    window = int(wildcards.window)
    step = int(float(wildcards.step) * window)
    cutoff_pct = float(wildcards.cutoff) * 100
    pop_id = wildcards.get("ppl") or wildcards.get("pair")

    method = wildcards.method
    if method == "delta_tajima_d":
        unit = "SNPs"
        method_name = "Delta Tajima's D"
    else:
        unit = "SNPs" if method == "moving_tajima_d" else "bp"
        method_name = format_method_name(method)

    if input and hasattr(input, "scores"):
        return " ".join([
            f"{pop_id}",
            f"(Window size={window} {unit},",
            f"Step size={step} {unit},",
            f"Top {cutoff_pct:.2f}%)",
        ])

    return " ".join([
        f"{pop_id}",
        method_name,
        f"(Window size={window} {unit},",
        f"Step size={step} {unit},",
        f"Top {cutoff_pct:.2f}%)",
    ])


def add_betascan_title(wildcards, input):
    """Generate title for BetaScan plots and tables."""
    cutoff_pct = float(wildcards.cutoff) * 100

    if hasattr(input, "scores"):
        return " ".join([
            f"{wildcards.ppl}",
            f"(Core Freq={wildcards.core_frq},",
            f"Top {cutoff_pct:.2f}%)",
        ])

    return " ".join([
        f"{wildcards.ppl}",
        "B1",
        f"(Core Freq={wildcards.core_frq},",
        f"Top {cutoff_pct:.2f}%)",
    ])


def add_dm_title(wildcards, input):
    """Add title for dadi-cli demographic model plots tables."""
    demog_fmt = wildcards.demog.replace('_', ' ').title()

    if hasattr(input, "tsv"):
        return f"{wildcards.ppl} {demog_fmt} Demographic Model (Top 10 Bestfits, {dadi_config['optimizations']} optimizations)"

    return f"{wildcards.ppl} {demog_fmt} Demographic Model Fit"


def add_dfe_title(wildcards, input):
    """Add title for dadi-cli DFE plots and tables."""
    demog_fmt = wildcards.demog.replace('_', ' ').title()
    dfe_fmt = wildcards.dfe.replace('_', ' ').title()
    base = f"{wildcards.ppl} {dfe_fmt} DFE ({demog_fmt})"

    if hasattr(input, "dfe_popt"):
        return f"{base} Mutation Proportions"

    if hasattr(input, "tsv"):
        if "godambe" in str(input.tsv):
            return f"{base} Estimated 95% Uncertainties ({dadi_config['bootstrap_replicates']} bootstrap replicates, chunk size={dadi_config['chunk_size']} bp)"
        return f"{base} (Top 10 Bestfits, {dadi_config['optimizations']} optimizations)"

    return f"{base} Model Fit"


def get_nomisid_flag(wildcards):
    """Return --nomisid if dataset has no ancestral alleles or dadi is folded."""
    if wildcards.dataset not in datasets_with_anc or not dadi_config["unfolded"]:
        return "--nomisid"
    return ""


def get_dadi_param(key, wildcards):
    """Return dadi parameter string, stripping misid value for folded datasets."""
    value = dadi_config[key]
    if wildcards.dataset not in datasets_with_anc or not dadi_config["unfolded"]:
        value = " ".join(value.split()[:-1])
    return value

def get_dfe_bestfit_files(wildcards):
    """Bestfit files for all populations in wildcards.dataset."""
    return [
        f
        for ds, sp, pop, rg in DATASET_1POP
        if ds == wildcards.dataset
        for f in expand(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}.InferDFE.bestfits",
            dataset=ds, species=sp, ppl=pop, ref_genome=rg,
            **DADI_1D_KW, dfe=dadi_config["dfe_1d"],
        )
    ]


def get_dfe_ci_files(wildcards):
    """Godambe CI files for all populations in wildcards.dataset."""
    return [
        f
        for ds, sp, pop, rg in DATASET_1POP
        if ds == wildcards.dataset
        for f in expand(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}.{ref_genome}.{demog}.{dfe}.godambe.ci",
            dataset=ds, species=sp, ppl=pop, ref_genome=rg,
            **DADI_1D_KW, dfe=dadi_config["dfe_1d"],
        )
    ]


def get_dfe_populations(wildcards):
    """Population list for wildcards.dataset, in DATASET_1POP order."""
    return [pop for ds, sp, pop, rg in DATASET_1POP if ds == wildcards.dataset]


def get_dfe_datasets(wildcards):
    """Dataset column values, same length as get_dfe_populations."""
    return [wildcards.dataset] * len(get_dfe_populations(wildcards))


def expand_1pop_circos(pattern, anc_only=False):
    source = DATASET_1POP_ANC if anc_only else DATASET_1POP
    return [
        f
        for ds, sp, pop, rg in source
        if ds in datasets_with_circos
        for f in expand(pattern, dataset=ds, species=sp, ppl=pop, ref_genome=rg)
    ]

def expand_2pop_circos(pattern, anc_only=False):
    source = DATASET_2POP_ANC if anc_only else DATASET_2POP
    return [
        f
        for ds, sp, pair, rg in source
        if ds in datasets_with_circos
        for f in expand(pattern, dataset=ds, species=sp, pair=pair, ref_genome=rg)
    ]


datasets_with_circos = [
    ds for ds, cfg in dataset_configs.items()
    if cfg.get("chr_bed") and cfg.get("cytoband")
]

def circos_labels(wildcards, type: str = "Circos Plot") -> dict[str, str]:
    """Labels for circos plots (positive and balancing selection)."""
    return {
        "Dataset": wildcards.dataset,
        "Population": wildcards.ppl,
        "Type": type,
    }

def get_chr_bed(wildcards):
    """Get chromosome sizes BED path for the given dataset (empty string if null)."""
    return get_dataset_cfg(wildcards).get("chr_bed") or ""


def get_cytoband(wildcards):
    """Get cytoband annotation path for the given dataset (empty string if null)."""
    return get_dataset_cfg(wildcards).get("cytoband") or ""


wildcard_constraints:
    cutoff=r"[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?",
    maf=r"[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?",
    window=r"[0-9_]+",
    step=r"[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?",
 
 
def expand_dataset(pattern, anc_only=False, **kw):
    source = DATASET_1POP_ANC if anc_only else DATASET_1POP
    datasets = {(ds, sp, rg) for ds, sp, _pop, rg in source}
    return [
        f
        for ds, sp, rg in sorted(datasets)
        for f in expand(pattern, dataset=ds, species=sp, ref_genome=rg, **kw)
    ]


def expand_2pop_focal(pattern, anc_only=False, **kw):
    """Like expand_2pop, but one entry per pair and focal population."""
    source = DATASET_2POP_ANC if anc_only else DATASET_2POP
    return [
        f
        for ds, sp, pair, rg in source
        for pop in pair.split("_")
        for f in expand(pattern, dataset=ds, species=sp, pair=pair, focal=pop, ref_genome=rg, **kw)
    ]


def get_xp_focal_directions(wildcards):
    """(pair, focal, reference) for both directions of every population pair in a dataset."""
    populations = dataset_configs[wildcards.dataset]["populations"]
    return [
        ("_".join((pop1, pop2)), focal, reference)
        for pop1, pop2 in combinations(populations, 2)
        for focal, reference in ((pop1, pop2), (pop2, pop1))
    ]
