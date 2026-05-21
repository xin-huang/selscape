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


ruleorder: plot_fitted_dfe > plot_fitted_1pop_dm

mask_singletons_flag = "--mask-singletons" if dadi_config["mask_singletons"] else ""


rule create_pop_info:
    input:
        metadata=get_metadata,
    output:
        pop_info="results/dadi/{species}/{dataset}/dfe/{ppl}/pop.list",
    log:
        "logs/deleterious_dfe/create_pop_info.{species}.{dataset}.{ppl}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        grep -w {wildcards.ppl} {input.metadata} > {output.pop_info} 2> {log}
        """


rule generate_1pop_fs:
    input:
        syn_vcf=lambda wc: f"results/{get_dadi_vcf_dir(wc)}/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/{wc.ppl}.biallelic.syn.snps.{wc.ref_genome}.vcf.gz",
        nonsyn_vcf=lambda wc: f"results/{get_dadi_vcf_dir(wc)}/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/{wc.ppl}.biallelic.nonsyn.snps.{wc.ref_genome}.vcf.gz",
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        syn_fs="results/dadi/{species}/{dataset}/dfe/{ppl}/fs/{ppl}.syn.{ref_genome}.dadi.fs",
        nonsyn_fs="results/dadi/{species}/{dataset}/dfe/{ppl}/fs/{ppl}.nonsyn.{ref_genome}.dadi.fs",
    params:
        polarized_flag=get_polarization_flag,
        mask_singletons_flag=f"{mask_singletons_flag}",
    log:
        "logs/deleterious_dfe/generate_1pop_fs.{species}.{dataset}.{ppl}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-ids {wildcards.ppl} --pop-info {input.pop_info} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --output {output.syn_fs} {params.polarized_flag} {params.mask_singletons_flag} 2> {log}
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-ids {wildcards.ppl} --pop-info {input.pop_info} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) --output {output.nonsyn_fs} {params.polarized_flag} {params.mask_singletons_flag} 2>> {log}
        """


rule infer_1pop_dm_warm_up:
    input:
        fs=rules.generate_1pop_fs.output.syn_fs,
    output:
        p0="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDM/{ppl}.{ref_genome}.{demog}.InferDM.opts.0",
    params:
        output_prefix="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDM/{ppl}.{ref_genome}.{demog}",
        p0=lambda wc: get_dadi_param("demog_1d_p0", wc),
        ubounds=lambda wc: get_dadi_param("demog_1d_ub", wc),
        lbounds=lambda wc: get_dadi_param("demog_1d_lb", wc),
        grid_size=dadi_config["dm_grid_size"],
        optimizations=dadi_config["optimizations"],
        nomisid_flag=get_nomisid_flag,
    resources:
        time=720,
        cpus=16,
    log:
        "logs/deleterious_dfe/infer_1pop_dm_warm_up.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --optimizations {params.optimizations} {params.nomisid_flag} 2> {log}
        """


rule infer_1pop_dm_fine_tune:
    input:
        fs=rules.generate_1pop_fs.output.syn_fs,
        p0=rules.infer_1pop_dm_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDM/{ppl}.{ref_genome}.{demog}.InferDM.bestfits",
    params:
        output_prefix="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDM/{ppl}.{ref_genome}.{demog}",
        ubounds=lambda wc: get_dadi_param("demog_1d_ub", wc),
        lbounds=lambda wc: get_dadi_param("demog_1d_lb", wc),
        grid_size=dadi_config["dm_grid_size"],
        optimizations=dadi_config["optimizations"],
        nomisid_flag=get_nomisid_flag,
    resources:
        time=720,
        cpus=16,
    log:
        "logs/deleterious_dfe/infer_1pop_dm_fine_tune.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus}  --force-convergence {params.optimizations} {params.nomisid_flag} 2> {log}
        """


rule get_1pop_dm_top_10_bestfits:
    input:
        bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
    output:
        bestfit="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDM/{ppl}.{ref_genome}.{demog}.InferDM.top10.bestfits.tsv",
    log:
        "logs/deleterious_dfe/get_1pop_dm_top_10_bestfits.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        awk '
            /^# Converged results/ {{on=1; next}}
            /^# Top 100 results/   {{exit}}
            on                     {{print; if(++n==10) exit}}
        ' {input.bestfit} > {output.bestfit} 2> {log}
        """


rule convert_1pop_dm_top_10_bestfits_html:
    input:
        tsv=rules.get_1pop_dm_top_10_bestfits.output.bestfit,
    output:
        html=report(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/html/{ppl}.{ref_genome}.{demog}.InferDM.top10.bestfits.html",
            category="Distribution of Fitness Effects",
            subcategory="Single Population",
            labels=lambda wildcards: fitted_1pop_dm_labels(
                wildcards, type="Bestfit Table"
            ),
        ),
    params:
        title=add_dm_title,
    log:
        "logs/deleterious_dfe/convert_1pop_dm_top_10_bestfits_html.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule generate_1d_cache:
    input:
        bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        cache="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.spectra.bpkl",
    params:
        grid_size=dadi_config["dfe_grid_size"],
        gamma_pts=dadi_config["gamma_pts"],
    resources:
        time=720,
        cpus=16,
    log:
        "logs/deleterious_dfe/generate_1d_cache.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli GenerateCache --model {wildcards.demog}_sel --demo-popt {input.bestfit} --sample-size $(awk 'END {{print NR * 2}}' {input.pop_info}) --grids {params.grid_size} --gamma-pts {params.gamma_pts} --output {output.cache} --cpus {resources.cpus} --cache-type cache1d 2> {log}
        """


rule infer_dfe_warm_up:
    input:
        fs=rules.generate_1pop_fs.output.nonsyn_fs,
        cache=rules.generate_1d_cache.output.cache,
        dm_bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
    output:
        p0="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}.InferDFE.opts.0",
    params:
        output_prefix="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}",
        p0=lambda wc: get_dadi_param("dfe_1d_p0", wc),
        ubounds=lambda wc: get_dadi_param("dfe_1d_ub", wc),
        lbounds=lambda wc: get_dadi_param("dfe_1d_lb", wc),
        ratio=dadi_config["ratio"],
        optimizations=dadi_config["optimizations"],
        nomisid_flag=get_nomisid_flag,
    resources:
        time=720,
        cpus=16,
    log:
        "logs/deleterious_dfe/infer_dfe_warm_up.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli InferDFE --fs {input.fs} --cache1d {input.cache} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf1d {wildcards.dfe} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --ratio {params.ratio} --cpus {resources.cpus} --optimizations {params.optimizations} {params.nomisid_flag} 2> {log}
        """


rule infer_dfe_fine_tune:
    input:
        fs=rules.generate_1pop_fs.output.nonsyn_fs,
        cache=rules.generate_1d_cache.output.cache,
        dm_bestfit=rules.infer_1pop_dm_fine_tune.output.bestfit,
        p0=rules.infer_dfe_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}.InferDFE.bestfits",
    params:
        output_prefix="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}",
        ubounds=lambda wc: get_dadi_param("dfe_1d_ub", wc),
        lbounds=lambda wc: get_dadi_param("dfe_1d_lb", wc),
        ratio=dadi_config["ratio"],
        optimizations=dadi_config["optimizations"],
        nomisid_flag=get_nomisid_flag,
    resources:
        time=720,
        cpus=16,
    log:
        "logs/deleterious_dfe/infer_dfe_fine_tune.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli InferDFE --fs {input.fs} --cache1d {input.cache} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} --pdf1d {wildcards.dfe} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --ratio {params.ratio} --cpus {resources.cpus} --force-convergence {params.optimizations} {params.nomisid_flag} 2> {log}
        """


rule get_1pop_dfe_top_10_bestfits:
    input:
        bestfit=rules.infer_dfe_fine_tune.output.bestfit,
    output:
        bestfit="results/dadi/{species}/{dataset}/dfe/{ppl}/InferDFE/{ppl}.{ref_genome}.{demog}.{dfe}.InferDFE.top10.bestfits.tsv",
    log:
        "logs/deleterious_dfe/get_1pop_dfe_top_10_bestfits.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        awk '
            /^# Converged results/ {{on=1; next}}
            /^# Top 100 results/   {{exit}}
            on                     {{print; if(++n==10) exit}}
        ' {input.bestfit} > {output.bestfit} 2> {log}
        """


rule convert_1pop_dfe_top_10_bestfits_html:
    input:
        tsv=rules.get_1pop_dfe_top_10_bestfits.output.bestfit,
    output:
        html=report(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/html/{ppl}.{ref_genome}.{demog}.{dfe}.InferDFE.top10.bestfits.html",
            category="Distribution of Fitness Effects",
            subcategory="Single Population",
            labels=lambda wildcards: fitted_dfe_labels(
                wildcards, type="Bestfit Table"
            ),
        ),
    params:
        title=add_dfe_title,
    log:
        "logs/deleterious_dfe/convert_1pop_dfe_top_10_bestfits_html.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule dfe_godambe_ci:
    input:
        syn_vcf=lambda wc: f"results/{get_dadi_vcf_dir(wc)}/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/{wc.ppl}.biallelic.syn.snps.{wc.ref_genome}.vcf.gz",
        nonsyn_vcf=lambda wc: f"results/{get_dadi_vcf_dir(wc)}/{wc.species}/{wc.dataset}/1pop/{wc.ppl}/{wc.ppl}.biallelic.nonsyn.snps.{wc.ref_genome}.vcf.gz",
        pop_info=rules.create_pop_info.output.pop_info,
        nonsyn_fs=rules.generate_1pop_fs.output.nonsyn_fs,
        cache=rules.generate_1d_cache.output.cache,
        dfe_bestfit=rules.infer_dfe_fine_tune.output.bestfit,
    output:
        dfe_godambe_ci="results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}.{ref_genome}.{demog}.{dfe}.godambe.ci",
    params:
        syn_dir="results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_syn",
        nonsyn_dir="results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_non",
        syn_output_prefix="results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_syn/{ppl}.{ref_genome}.syn",
        nonsyn_output_prefix="results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}_bootstrapping_non/{ppl}.{ref_genome}.nonsyn",
        bootstrap_reps=dadi_config["bootstrap_replicates"],
        chunk_size=dadi_config["chunk_size"],
        polarized_flag=get_polarization_flag,
        mask_singletons_flag=f"{mask_singletons_flag}",
    log:
        "logs/deleterious_dfe/dfe_godambe_ci.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( [ -d {params.syn_dir} ] || mkdir -p {params.syn_dir} ) 2> {log}
        ( [ -d {params.nonsyn_dir} ] || mkdir -p {params.nonsyn_dir} ) 2>> {log}
        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-info {input.pop_info} --pop-ids {wildcards.ppl} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) {params.polarized_flag} {params.mask_singletons_flag} --bootstrap {params.bootstrap_reps} --chunk-size {params.chunk_size} --output {params.syn_output_prefix} 2>> {log}
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-info {input.pop_info} --pop-ids {wildcards.ppl} --projections $(awk 'END {{print NR * 2}}' {input.pop_info}) {params.polarized_flag} {params.mask_singletons_flag} --bootstrap {params.bootstrap_reps} --chunk-size {params.chunk_size} --output {params.nonsyn_output_prefix} 2>> {log}
        dadi-cli StatDFE --fs {input.nonsyn_fs} --dfe-popt {input.dfe_bestfit} --cache1d {input.cache} --pdf1d {wildcards.dfe} --bootstrapping-nonsynonymous-dir {params.nonsyn_dir} --bootstrapping-synonymous-dir {params.syn_dir} --output {output.dfe_godambe_ci}  2>> {log}
        """


rule parse_dfe_godambe_ci_table:
    input:
        bestfit=rules.infer_dfe_fine_tune.output.bestfit,
        ci=rules.dfe_godambe_ci.output.dfe_godambe_ci,
    output:
        tsv="results/dadi/{species}/{dataset}/dfe/{ppl}/StatDFE/{ppl}.{ref_genome}.{demog}.{dfe}.godambe.ci.tsv",
    log:
        "logs/deleterious_dfe/parse_dfe_godambe_ci_table.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        awk -F "\\t" -v OFS="\\t" '/^# Log\(likelihood\)/ {{
            printf "step_size\\t"
            for (i=2; i<NF; i++) printf "%s%s", $i, (i<NF-1?OFS:ORS)
            exit
        }}' {input.bestfit} > {output.tsv} 2> {log}

        sed -n '/^Estimated 95% uncerts/ {{
            :join
            /]/!{{N; b join}}
            s/\\n/ /g
            p
        }}' {input.ci} |\
        awk '{{
            match($0, /step size ([0-9.]+)/, m)
            step = m[1]
            match($0, /\[([0-9.eE+\- ]+)\]/, a)
            n = split(a[1], arr)
            printf "%s", step
            for (i = 1; i <= n; i++) printf "\t%s", arr[i]
            print ""
        }}' >> {output.tsv} 2>> {log}
        """


rule dfe_godambe_ci_table_html:
    input:
        tsv=rules.parse_dfe_godambe_ci_table.output.tsv,
    output:
        html=report(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/html/{ppl}.{ref_genome}.{demog}.{dfe}.godambe.ci.html",
            category="Distribution of Fitness Effects",
            subcategory="Single Population",
            labels=lambda wildcards: fitted_dfe_labels(
                wildcards, type="Estimated 95% Uncerts"
            ),
        ),
    params:
        title=add_dfe_title,
    log:
        "logs/reports/dfe_ci_table_html.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule plot_fitted_1pop_dm:
    input:
        fs=rules.generate_1pop_fs.output.syn_fs,
        dm_popt=rules.infer_1pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        fs_plot="results/dadi/{species}/{dataset}/dfe/{ppl}/plots/{ppl}.{ref_genome}.{demog}.fitted.png",
    params:
        ploidy=get_ploidy,
    log:
        "logs/deleterious_dfe/plot_fitted_dm.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli Plot --fs {input.fs} --demo-popt {input.dm_popt} --model {wildcards.demog} --projections $(awk 'END {{print NR * {params.ploidy}}}' {input.pop_info}) --output {output.fs_plot} 2> {log}
        """


rule plot_fitted_dfe:
    input:
        fs=rules.generate_1pop_fs.output.nonsyn_fs,
        dfe_popt=rules.infer_dfe_fine_tune.output.bestfit,
        cache=rules.generate_1d_cache.output.cache,
        pop_info=rules.create_pop_info.output.pop_info,
    output:
        fs_plot="results/dadi/{species}/{dataset}/dfe/{ppl}/plots/{ppl}.{ref_genome}.{demog}.{dfe}.fitted.png",
    params:
        ploidy=get_ploidy,
    log:
        "logs/deleterious_dfe/plot_fitted_dfe.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        dadi-cli Plot --fs {input.fs} --cache1d {input.cache} --dfe-popt {input.dfe_popt} --pdf1d {wildcards.dfe} --projections $(awk 'END {{print NR * {params.ploidy}}}' {input.pop_info}) --output {output.fs_plot} 2> {log}
        """


rule plot_mutation_proportions:
    input:
        dfe_popt=rules.infer_dfe_fine_tune.output.bestfit,
    output:
        plot=report(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/plots/{ppl}.{ref_genome}.{demog}.{dfe}.fitted.mut.prop.png",
            category="Distribution of Fitness Effects",
            subcategory="Single Population",
            labels=lambda wildcards: fitted_dfe_labels(
                wildcards, type="Proportion Plot"
            ),
        ),
    params:
        title=add_dfe_title,
    log:
        "logs/deleterious_dfe/plot_mutation_proportion.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_mutation_prop.py"


rule wrap_fitted_1pop_dm_html:
    input:
        plot=rules.plot_fitted_1pop_dm.output.fs_plot,
    output:
        html=report(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/html/{ppl}.{ref_genome}.{demog}.dm.fitted.html",
            category="Distribution of Fitness Effects",
            subcategory="Single Population",
            labels=lambda wildcards: fitted_1pop_dm_labels(wildcards, type="Model Fit Plot"),
        ),
    params:
        title=add_dm_title,
    log:
        "logs/deleterious_dfe/wrap_fitted_1pop_dm_html.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot2html.py"


rule wrap_fitted_dfe_html:
    input:
        plot=rules.plot_fitted_dfe.output.fs_plot,
    output:
        html=report(
            "results/dadi/{species}/{dataset}/dfe/{ppl}/html/{ppl}.{ref_genome}.{demog}.{dfe}.dfe.fitted.html",
            category="Distribution of Fitness Effects",
            subcategory="Single Population",
            labels=lambda wildcards: fitted_dfe_labels(wildcards, type="Model Fit Plot"),
        ),
    params:
        title=add_dfe_title,
    log:
        "logs/deleterious_dfe/wrap_fitted_dfe_html.{species}.{dataset}.{ppl}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot2html.py"
