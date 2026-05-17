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


rule merge_dfe_confidence_intervals:
    input:
        bestfit_files=expand_1pop(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}.InferDFE.bestfits",
            **DADI_1D_KW, dfe=dadi_config["dfe_1d"],
        ),
        ci_files=expand_1pop(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}.{ref_genome}.{demog}.{dfe}.godambe.ci",
            **DADI_1D_KW, dfe=dadi_config["dfe_1d"],
        ),
    output:
        merged="results/plots/dfe/combined.dfe_params.tsv",
    params:
        populations=expand_1pop("{ppl}"),
        datasets=expand_1pop("{dataset}"),
    log:
        "logs/deleterious_dfe/merge_dfe_confidence_intervals.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/merge_dfe_ci.py"


rule plot_dfe_confidence_intervals:
    input:
        data=rules.merge_dfe_confidence_intervals.output.merged,
    output:
        plot=report(
            "results/plots/dfe/combined.dfe_params.svg",
            category="Distribution of Fitness Effects",
            subcategory="DFE Parameters",
            labels={"Type": "DFE Confidence Intervals"},
        ),
    params:
        populations=expand_1pop("{ppl}"),
        population_groups=lambda _: main_config.get("population_groups", {}),
        mu_ylim=None,
        sigma_ylim=None,
    log:
        "logs/deleterious_dfe/plot_dfe_confidence_intervals.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_dfe_params.py"
