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


SELSCAN_KW = dict(
    species=main_config["species"],
    ppl=main_config["populations"],
    method=selscan_config["wp_stats"],
    maf=selscan_config["maf"],
    cutoff=selscan_config["top_proportion"],
)

SELSCAN_XP_KW = dict(
    species=main_config["species"],
    pair=["_".join(pair) for pair in combinations(main_config["populations"], 2)],
    method=selscan_config["xp_stats"],
    maf=selscan_config["maf"],
    cutoff=selscan_config["top_proportion"],
)

BETASCAN_KW = dict(
    species=main_config["species"],
    ppl=main_config["populations"],
    ref_genome=main_config["ref_genome"],
    core_frq=betascan_config["core_frq"],
    cutoff=betascan_config["top_proportion"],
)

DADI_1D_KW = dict(
    species=main_config["species"],
    ppl=main_config["populations"],
    ref_genome=main_config["ref_genome"],
    demog=dadi_config["demog_1d"],
)

TAJIMAD_MOVING_KW = dict(
    species=main_config["species"],
    ppl=main_config["populations"],
    method="moving_tajima_d",
    window=scikit_allel_config["mtjd_window_sizes"],
    step=scikit_allel_config["mtjd_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

TAJIMAD_WINDOWED_KW = dict(
    species=main_config["species"],
    ppl=main_config["populations"],
    method="windowed_tajima_d",
    window=scikit_allel_config["wtjd_window_sizes"],
    step=scikit_allel_config["wtjd_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

DELTA_TAJIMAD_KW = dict(
    species=main_config["species"],
    pair=["_".join(pair) for pair in combinations(main_config["populations"], 2)],
    method="delta_tajima_d",
    window=scikit_allel_config["dtjd_window_sizes"],
    step=scikit_allel_config["dtjd_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

selscan_method_names = {
    "ihs": "iHS",
    "nsl": "nSL",
    "xpehh": "XP-EHH",
    "xpnsl": "XP-nSL",
}

# HELPER FUNCTIONS


def get_anc_allele_bed(wildcards):
    """Get ancestral allele bed files."""
    return f"{main_config['anc_alleles']['path']}/{main_config['anc_alleles']['prefix']}{wildcards.i}.bed.gz"


def get_vcf_input_path(wildcards):
    """Get VCF input file path."""
    return f"{main_config['data_folder']}/{main_config['vcf_prefix']}{wildcards.i}{main_config['vcf_suffix']}"


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

def selscan_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for within-population selscan Manhattan plot."""
    return {
        "Population": wildcards.ppl,
        "Minor Allele Frequency": wildcards.maf,
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }


def selscan_xp_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for cross-population selscan Manhattan plot."""
    return {
        "Populations": _vs_pair(wildcards),
        "Minor Allele Frequency": wildcards.maf,
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }


def betascan_labels(wildcards, type: str = "Manhattan Plot") -> dict[str, str]:
    """Labels for betascan plots (Manhattan or Enrichment), includes core frequency."""
    return {
        "Population": wildcards.ppl,
        "Core Frequency": str(wildcards.core_frq),
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }


def tajima_d_labels(wildcards, type: str = "Plot") -> dict[str, str]:
    """Labels for Tajima's D plots (both windowed and moving)."""
    method_name = (
        "Moving Tajima's D"
        if wildcards.method.startswith("moving")
        else "Windowed Tajima's D"
    )
    window_unit = " SNPs" if wildcards.method.startswith("moving") else " bp"
    step_size = int(float(wildcards.step) * int(wildcards.window))
    return {
        "Population": wildcards.ppl,
        "Window": f"{wildcards.window}{window_unit}",
        "Step": f"{step_size}{window_unit}",
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }


def delta_tajima_d_labels(wildcards, type: str = "Plot") -> dict[str, str]:
    """Labels for delta Tajima's D plots (cross-population)."""
    step_size = int(float(wildcards.step) * int(wildcards.window))
    return {
        "Populations": _vs_pair(wildcards),
        "Window": f"{wildcards.window} SNPs",
        "Step": f"{step_size} SNPs",
        "Threshold": _top_pct(wildcards),
        "Type": type,
    }


def fitted_1pop_dm_labels(wildcards, type: str = "Model Fit") -> dict[str, str]:
    """Labels for 1-population demographic model fit plots."""
    return {
        "Population": wildcards.ppl,
        "Demographic Model": wildcards.demog,
        "Type": type,
    }


def fitted_dfe_labels(wildcards, type: str = "Model Fit") -> dict[str, str]:
    """Labels for 1-population DFE model fit plots."""
    return {
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
    """ Generate title for BetaScan plots and tables."""
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
