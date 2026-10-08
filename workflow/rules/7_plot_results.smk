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
            category="Positive Selection",
            subcategory="Circos Plots",
            labels=circos_labels,
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
            category="Balancing Selection",
            subcategory="Circos Plots",
            labels=circos_labels,
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


rule plot_selscan_upset:
    input:
        genes=lambda wc: expand(
            "results/positive_selection/selscan/{species}/{dataset}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_upset/{species}/{dataset}/positive_selection/{dataset}.{method}.maf_{maf}.top_{cutoff}.upset.svg",
            category="Within-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Positive",
                "Threshold": _top_pct(wildcards),
                "Type": "UpSet (shared outliers between groups)",
            },
        ),
        table="results/plots/wp_upset/{species}/{dataset}/positive_selection/{dataset}.{method}.maf_{maf}.top_{cutoff}.upset.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        max_intersections=15,
        title=lambda wc: (
            f"{wc.dataset} {selscan_method_names[wc.method]} "
            f"(MAF={wc.maf}, Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_selscan_upset.{species}.{dataset}.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_upset.py"


rule plot_tajima_d_upset:
    input:
        genes=lambda wc: expand(
            "results/positive_selection/scikit-allel/{species}/{dataset}/1pop/{ppl}/{method}/{window}_{step}/{ppl}.{method}.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_upset/{species}/{dataset}/positive_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.upset.svg",
            category="Within-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Positive",
                "Window": f"{wildcards.window} {'SNPs' if wildcards.method == 'moving_tajima_d' else 'bp'}",
                "Threshold": _top_pct(wildcards),
                "Type": "UpSet (shared outliers between groups)",
            },
        ),
        table="results/plots/wp_upset/{species}/{dataset}/positive_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.upset.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        max_intersections=15,
        title=lambda wc: (
            f"{wc.dataset} {format_method_name(wc.method)}, positive selection "
            f"(Window size={wc.window} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
            f"Step size={int(float(wc.step) * int(wc.window))} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
            f"Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_tajima_d_upset.{species}.{dataset}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_upset.py"


rule plot_betascan_upset:
    input:
        genes=lambda wc: expand(
            "results/balancing_selection/betascan/{species}/{dataset}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            ref_genome=get_ref_genome(wc),
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_upset/{species}/{dataset}/balancing_selection/{dataset}.m_{core_frq}.b1.top_{cutoff}.upset.svg",
            category="Within-Population Overview",
            subcategory="B1",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Balancing",
                "Core Frequency": wildcards.core_frq,
                "Threshold": _top_pct(wildcards),
                "Type": "UpSet (shared outliers between groups)",
            },
        ),
        table="results/plots/wp_upset/{species}/{dataset}/balancing_selection/{dataset}.m_{core_frq}.b1.top_{cutoff}.upset.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        max_intersections=15,
        title=lambda wc: (
            f"{wc.dataset} B1 "
            f"(Core Freq={wc.core_frq}, Top {float(wc.cutoff) * 100:.2f}%)"
rule make_xp_selection_circos:
    input:
        unpack(get_xp_circos_inputs),
        chr_bed=get_chr_bed,
        cytoband=get_cytoband,
    output:
        plot=report(
            "results/plots/circos/{species}/{dataset}/{pair}/{pair}_xp_selection_circos_scores.png",
            category="Positive Selection",
            subcategory="Circos Plots",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Population": wildcards.pair,
                "Type": "Circos Plot",
            },
        ),
    params:
        population="{pair}",
        ref_genome=get_ref_genome,
        tracks=XP_CIRCOS_TRACKS,
    resources:
        mem_mb=32000,
    log:
        "logs/circos/make_xp_selection_circos.{species}.{dataset}.{pair}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_circos_scores.py"


rule plot_selscan_xp_upset:
    input:
        genes=lambda wc: [
            f
            for pair, focal, _reference in get_xp_focal_directions(wc)
            for f in expand(
                "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.focal_{focal}.outlier.genes",
                pair=pair, focal=focal, allow_missing=True,
            )
        ],
    output:
        plot=report(
            "results/plots/xp_upset/{species}/{dataset}/{dataset}.{method}.maf_{maf}.top_{cutoff}.upset.svg",
            category="Cross-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Threshold": _top_pct(wildcards),
                "Type": "UpSet (shared candidates between groups)",
            },
        ),
        table="results/plots/xp_upset/{species}/{dataset}/{dataset}.{method}.maf_{maf}.top_{cutoff}.upset.tsv",
    params:
        focal_pairs=lambda wc: [[focal, reference] for _pair, focal, reference in get_xp_focal_directions(wc)],
        population_groups=lambda _: main_config.get("population_groups", {}),
        max_intersections=15,
        title=lambda wc: (
            f"{wc.dataset} {selscan_method_names[wc.method]} "
            f"(MAF={wc.maf}, Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_betascan_upset.{species}.{dataset}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_upset.py"


rule plot_tajima_d_balancing_upset:
    input:
        genes=lambda wc: expand(
            "results/balancing_selection/scikit-allel/{species}/{dataset}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_upset/{species}/{dataset}/balancing_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.upset.svg",
            category="Within-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Balancing",
                "Window": f"{wildcards.window} {'SNPs' if wildcards.method == 'moving_tajima_d' else 'bp'}",
                "Threshold": _top_pct(wildcards),
                "Type": "UpSet (shared outliers between groups)",
            },
        ),
        table="results/plots/wp_upset/{species}/{dataset}/balancing_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.upset.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        max_intersections=15,
        title=lambda wc: (
            f"{wc.dataset} {format_method_name(wc.method)}, balancing selection "
            f"(Window size={wc.window} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
            f"Step size={int(float(wc.step) * int(wc.window))} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
        "logs/plots/plot_selscan_xp_upset.{species}.{dataset}.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_xp_upset.py"


rule plot_delta_tajima_d_upset:
    input:
        genes=lambda wc: [
            f
            for pair, focal, _reference in get_xp_focal_directions(wc)
            for f in expand(
                "results/positive_selection/scikit-allel/{species}/{dataset}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.focal_{focal}.outlier.genes",
                pair=pair, focal=focal, allow_missing=True,
            )
        ],
    output:
        plot=report(
            "results/plots/xp_upset/{species}/{dataset}/{dataset}.{method}.{window}_{step}.top_{cutoff}.upset.svg",
            category="Cross-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Window": f"{wildcards.window} SNPs",
                "Threshold": _top_pct(wildcards),
                "Type": "UpSet (shared candidates between groups)",
            },
        ),
        table="results/plots/xp_upset/{species}/{dataset}/{dataset}.{method}.{window}_{step}.top_{cutoff}.upset.tsv",
    params:
        focal_pairs=lambda wc: [[focal, reference] for _pair, focal, reference in get_xp_focal_directions(wc)],
        population_groups=lambda _: main_config.get("population_groups", {}),
        max_intersections=15,
        title=lambda wc: (
            f"{wc.dataset} Delta Tajima's D "
            f"(Window size={wc.window} SNPs, "
            f"Step size={int(float(wc.step) * int(wc.window))} SNPs, "
            f"Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_tajima_d_balancing_upset.{species}.{dataset}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_upset.py"


rule plot_selscan_wp_matrix:
    input:
        genes=lambda wc: expand(
            "results/positive_selection/selscan/{species}/{dataset}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_matrix/{species}/{dataset}/positive_selection/{dataset}.{method}.maf_{maf}.top_{cutoff}.matrix.svg",
            category="Within-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Positive",
                "Threshold": _top_pct(wildcards),
                "Type": "Matrix (shared outlier genes)",
            },
        ),
        table="results/plots/wp_matrix/{species}/{dataset}/positive_selection/{dataset}.{method}.maf_{maf}.top_{cutoff}.matrix.tsv",
    params:
        "logs/plots/plot_delta_tajima_d_upset.{species}.{dataset}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_xp_upset.py"


rule plot_selscan_xp_matrix:
    input:
        genes=lambda wc: [
            f
            for pair, focal, _reference in get_xp_focal_directions(wc)
            for f in expand(
                "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.focal_{focal}.outlier.genes",
                pair=pair, focal=focal, allow_missing=True,
            )
        ],
    output:
        plot=report(
            "results/plots/xp_matrix/{species}/{dataset}/{dataset}.{method}.maf_{maf}.top_{cutoff}.matrix.svg",
            category="Cross-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Threshold": _top_pct(wildcards),
                "Type": "Matrix (candidates per direction)",
            },
        ),
        table="results/plots/xp_matrix/{species}/{dataset}/{dataset}.{method}.maf_{maf}.top_{cutoff}.matrix.tsv",
    params:
        focal_pairs=lambda wc: [[focal, reference] for _pair, focal, reference in get_xp_focal_directions(wc)],
        population_groups=lambda _: main_config.get("population_groups", {}),
        title=lambda wc: (
            f"{wc.dataset} {selscan_method_names[wc.method]} "
            f"(MAF={wc.maf}, Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_selscan_wp_matrix.{species}.{dataset}.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_wp_matrix.py"


rule plot_tajima_d_wp_matrix:
    input:
        genes=lambda wc: expand(
            "results/positive_selection/scikit-allel/{species}/{dataset}/1pop/{ppl}/{method}/{window}_{step}/{ppl}.{method}.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_matrix/{species}/{dataset}/positive_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.matrix.svg",
            category="Within-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Positive",
                "Window": f"{wildcards.window} {'SNPs' if wildcards.method == 'moving_tajima_d' else 'bp'}",
                "Threshold": _top_pct(wildcards),
                "Type": "Matrix (shared outlier genes)",
            },
        ),
        table="results/plots/wp_matrix/{species}/{dataset}/positive_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.matrix.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        title=lambda wc: (
            f"{wc.dataset} {format_method_name(wc.method)}, positive selection "
            f"(Window size={wc.window} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
            f"Step size={int(float(wc.step) * int(wc.window))} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
            f"Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_tajima_d_wp_matrix.{species}.{dataset}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_wp_matrix.py"


rule plot_betascan_wp_matrix:
    input:
        genes=lambda wc: expand(
            "results/balancing_selection/betascan/{species}/{dataset}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            ref_genome=get_ref_genome(wc),
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_matrix/{species}/{dataset}/balancing_selection/{dataset}.m_{core_frq}.b1.top_{cutoff}.matrix.svg",
            category="Within-Population Overview",
            subcategory="B1",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Balancing",
                "Core Frequency": wildcards.core_frq,
                "Threshold": _top_pct(wildcards),
                "Type": "Matrix (shared outlier genes)",
            },
        ),
        table="results/plots/wp_matrix/{species}/{dataset}/balancing_selection/{dataset}.m_{core_frq}.b1.top_{cutoff}.matrix.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        title=lambda wc: (
            f"{wc.dataset} B1 "
            f"(Core Freq={wc.core_frq}, Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_betascan_wp_matrix.{species}.{dataset}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_wp_matrix.py"


rule plot_tajima_d_balancing_wp_matrix:
    input:
        genes=lambda wc: expand(
            "results/balancing_selection/scikit-allel/{species}/{dataset}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.outlier.genes",
            ppl=get_dataset_cfg(wc)["populations"],
            allow_missing=True,
        ),
    output:
        plot=report(
            "results/plots/wp_matrix/{species}/{dataset}/balancing_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.matrix.svg",
            category="Within-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Selection": "Balancing",
                "Window": f"{wildcards.window} {'SNPs' if wildcards.method == 'moving_tajima_d' else 'bp'}",
                "Threshold": _top_pct(wildcards),
                "Type": "Matrix (shared outlier genes)",
            },
        ),
        table="results/plots/wp_matrix/{species}/{dataset}/balancing_selection/{dataset}.{method}.{window}_{step}.top_{cutoff}.matrix.tsv",
    params:
        population_groups=lambda _: main_config.get("population_groups", {}),
        title=lambda wc: (
            f"{wc.dataset} {format_method_name(wc.method)}, balancing selection "
            f"(Window size={wc.window} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
            f"Step size={int(float(wc.step) * int(wc.window))} {'SNPs' if wc.method == 'moving_tajima_d' else 'bp'}, "
        "logs/plots/plot_selscan_xp_matrix.{species}.{dataset}.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_xp_matrix.py"


rule plot_delta_tajima_d_matrix:
    input:
        genes=lambda wc: [
            f
            for pair, focal, _reference in get_xp_focal_directions(wc)
            for f in expand(
                "results/positive_selection/scikit-allel/{species}/{dataset}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.focal_{focal}.outlier.genes",
                pair=pair, focal=focal, allow_missing=True,
            )
        ],
    output:
        plot=report(
            "results/plots/xp_matrix/{species}/{dataset}/{dataset}.{method}.{window}_{step}.top_{cutoff}.matrix.svg",
            category="Cross-Population Overview",
            subcategory="{method}",
            labels=lambda wildcards: {
                "Dataset": wildcards.dataset,
                "Window": f"{wildcards.window} SNPs",
                "Threshold": _top_pct(wildcards),
                "Type": "Matrix (candidates per direction)",
            },
        ),
        table="results/plots/xp_matrix/{species}/{dataset}/{dataset}.{method}.{window}_{step}.top_{cutoff}.matrix.tsv",
    params:
        focal_pairs=lambda wc: [[focal, reference] for _pair, focal, reference in get_xp_focal_directions(wc)],
        population_groups=lambda _: main_config.get("population_groups", {}),
        title=lambda wc: (
            f"{wc.dataset} Delta Tajima's D "
            f"(Window size={wc.window} SNPs, "
            f"Step size={int(float(wc.step) * int(wc.window))} SNPs, "
            f"Top {float(wc.cutoff) * 100:.2f}%)"
        ),
    resources:
        mem_mb=8000,
    log:
        "logs/plots/plot_tajima_d_balancing_wp_matrix.{species}.{dataset}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_wp_matrix.py"
        "logs/plots/plot_delta_tajima_d_matrix.{species}.{dataset}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_xp_matrix.py"
