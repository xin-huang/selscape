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
    window=scikit_allel_config["moving_window_sizes"],
    step=scikit_allel_config["moving_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

TAJIMAD_WINDOWED_KW = dict(
    species=main_config["species"],
    ppl=main_config["populations"],
    method="windowed_tajima_d",
    window=scikit_allel_config["windowed_window_sizes"],
    step=scikit_allel_config["windowed_step_size_ratios"],
    cutoff=scikit_allel_config["top_proportion"],
)

DELTA_TAJIMAD_KW = dict(
    species=main_config["species"],
    pair=["_".join(pair) for pair in combinations(main_config["populations"], 2)],
    method="delta_moving_tajima_d",
    window=scikit_allel_config["delta_moving_window_sizes"],
    step=scikit_allel_config["delta_moving_step_size_ratios"],
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
