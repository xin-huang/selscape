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
        bestfit_files=get_dfe_bestfit_files,
        ci_files=get_dfe_ci_files,
    output:
        merged="results/plots/dfe/{species}/{dataset}/{dataset}.dfe_params.tsv",
    params:
        populations=get_dfe_populations,
        datasets=get_dfe_datasets,
    log:
        "logs/deleterious_dfe/merge_dfe_confidence_intervals.{species}.{dataset}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/merge_dfe_ci.py"


rule plot_dfe_confidence_intervals:
    input:
        data=rules.merge_dfe_confidence_intervals.output.merged,
    output:
        plot=report(
            "results/plots/dfe/{species}/{dataset}/{species}.{dataset}.dfe_params.svg",
            category="Distribution of Fitness Effects",
            subcategory="DFE Parameters",
            labels={"Dataset": "{dataset}", "Type": "DFE Confidence Intervals"},
        ),
    params:
        populations=get_dfe_populations,
        population_groups=lambda _: main_config.get("population_groups", {}),
        mu_ylim=None,
        sigma_ylim=None,
    log:
        "logs/deleterious_dfe/plot_dfe_confidence_intervals.{species}.{dataset}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_dfe_params.py"


rule make_positive_selection_circos:
    input:
        ihs_scores=lambda wc: f"results/positive_selection/selscan/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/ihs_{SELSCAN_KW['maf']}/{wc.ppl}.normalized.ihs.scores",
        nsl_scores=lambda wc: f"results/positive_selection/selscan/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/nsl_{SELSCAN_KW['maf']}/{wc.ppl}.normalized.nsl.scores",
        mtjd_scores=lambda wc: f"results/positive_selection/scikit-allel/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/moving_tajima_d/{TAJIMAD_MOVING_KW['window'][0]}_{TAJIMAD_MOVING_KW['step'][0]}/{wc.ppl}.moving_tajima_d.scores",
        wtjd_scores=lambda wc: f"results/positive_selection/scikit-allel/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/windowed_tajima_d/{TAJIMAD_WINDOWED_KW['window'][0]}_{TAJIMAD_WINDOWED_KW['step'][0]}/{wc.ppl}.windowed_tajima_d.scores",
        chr_bed=get_chr_bed,
        cytoband=get_cytoband,
    output:
        plot=report(
            "results/plots/circos/{species}/{dataset}/{ppl}/{ppl}_positive_selection_circos_scores.png",
            category="Circos Plots",
            subcategory="Positive Selection",
            labels=lambda wildcards: {
                "Population": wildcards.ppl,
                "Type": "Circos Plot",
            },
        ),
    params:
        population="{ppl}",
        ref_genome=get_ref_genome,
        tracks=[
            {
                "name": "iHS",
                "file": "ihs_scores",
                "score_col": "normalized_ihs",
                "r_range": [65, 75],
                "color": "#1f77b4",
            },
            {
                "name": "nSL",
                "file": "nsl_scores",
                "score_col": "normalized_nsl",
                "r_range": [50, 60],
                "color": "#ff7f0e",
            },
            {
                "name": "mtjd",
                "file": "mtjd_scores",
                "score_col": "tajima_d",
                "r_range": [35, 45],
                "color": "#2ca02c",
            },
            {
                "name": "wtjd",
                "file": "wtjd_scores",
                "score_col": "tajima_d",
                "r_range": [20, 30],
                "color": "#d62728",
            },
        ],
    resources:
        mem_mb=32000,
    log:
        "logs/circos/make_positive_selection_circos.{species}.{dataset}.{ppl}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_circos_scores.py"


rule make_balancing_selection_circos:
    input:
        b1_scores=lambda wc: f"results/balancing_selection/betascan/{wc.species}/{wc.dataset}/{wc.ppl}/m_{BETASCAN_KW['core_frq']}/{wc.ppl}.{get_ref_genome(wc)}.m_{BETASCAN_KW['core_frq']}.b1.scores",
        mtjd_bal_scores=lambda wc: f"results/balancing_selection/scikit-allel/{wc.species}/{wc.dataset}/moving_tajima_d/{wc.ppl}/{TAJIMAD_MOVING_KW['window'][0]}_{TAJIMAD_MOVING_KW['step'][0]}/{wc.ppl}.moving_tajima_d.merged.scores",
        wtjd_bal_scores=lambda wc: f"results/balancing_selection/scikit-allel/{wc.species}/{wc.dataset}/windowed_tajima_d/{wc.ppl}/{TAJIMAD_WINDOWED_KW['window'][0]}_{TAJIMAD_WINDOWED_KW['step'][0]}/{wc.ppl}.windowed_tajima_d.merged.scores",
        chr_bed=get_chr_bed,
        cytoband=get_cytoband,
    output:
        plot=report(
            "results/plots/circos/{species}/{dataset}/{ppl}/{ppl}_balancing_selection_circos_scores.png",
            category="Circos Plots",
            subcategory="Balancing Selection",
            labels=lambda wildcards: {
                "Population": wildcards.ppl,
                "Type": "Circos Plot",
            },
        ),
    params:
        population="{ppl}",
        ref_genome=get_ref_genome,
        tracks=[
            {
                "name": "B1",
                "file": "b1_scores",
                "score_col": "B1",
                "r_range": [60, 75],
                "color": "#1f77b4",
            },
            {
                "name": "mtjd",
                "file": "mtjd_bal_scores",
                "score_col": "tajima_d",
                "r_range": [40, 55],
                "color": "#2ca02c",
            },
            {
                "name": "wtjd",
                "file": "wtjd_bal_scores",
                "score_col": "tajima_d",
                "r_range": [20, 35],
                "color": "#d62728",
            },
        ],
    resources:
        mem_mb=32000,
    log:
        "logs/circos/make_balancing_selection_circos.{species}.{dataset}.{ppl}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_circos_scores.py"
