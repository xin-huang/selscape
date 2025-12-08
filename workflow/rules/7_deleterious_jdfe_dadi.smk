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


ruleorder: plot_fitted_jdfe > plot_fitted_2pop_dm


vcf_dir = "processed_data"
polarization_flag = ""
if "anc_alleles" in main_config and main_config["anc_alleles"] and dadi_config["unfolded"]:
    vcf_dir = "polarized_data"
    polarization_flag = "--polarized"
mask_singletons_flag = "--mask-singletons" if dadi_config["mask_singletons"] else ""


rule generate_2pop_fs:
    input:
        syn_vcf=f"results/{vcf_dir}/{{species}}/2pop/{{pair}}/{{pair}}.biallelic.syn.snps.{{ref_genome}}.vcf.gz",
        nonsyn_vcf=f"results/{vcf_dir}/{{species}}/2pop/{{pair}}/{{pair}}.biallelic.nonsyn.snps.{{ref_genome}}.vcf.gz",
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        syn_fs="results/dadi/{species}/jdfe/{pair}/fs/{pair}.syn.{ref_genome}.dadi.fs",
        nonsyn_fs="results/dadi/{species}/jdfe/{pair}/fs/{pair}.nonsyn.{ref_genome}.dadi.fs",
    params:
        ploidy=main_config["ploidy"],
        polarized_flag=f'{polarization_flag}',
        mask_singletons_flag=f'{mask_singletons_flag}',
    log:
        "logs/deleterious_jdfe/generate_2pop_fs.{species}.{pair}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"

        dadi-cli GenerateFs \
            --vcf {input.syn_vcf} \
            --pop-ids $pop1 $pop2 \
            --pop-info {input.pop_info} \
            --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') \
            --output {output.syn_fs} \
            {params.polarized_flag} {params.mask_singletons_flag} 2> {log}

        dadi-cli GenerateFs \
            --vcf {input.nonsyn_vcf} \
            --pop-ids $pop1 $pop2 \
            --pop-info {input.pop_info} \
            --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') \
            --output {output.nonsyn_fs} \
            {params.polarized_flag} {params.mask_singletons_flag} 2>> {log}
        """


rule infer_2pop_dm_warm_up:
    input:
        fs=rules.generate_2pop_fs.output.syn_fs,
    output:
        p0="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{ref_genome}.{demog}.InferDM.opts.0",
    params:
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{ref_genome}.{demog}",
        p0=dadi_config["demog_2d_p0"],
        ubounds=dadi_config["demog_2d_ub"],
        lbounds=dadi_config["demog_2d_lb"],
        grid_size=dadi_config["dm_grid_size"],
        optimizations=dadi_config["optimizations"],
    resources:
        time=1440,
        cpus=16,
    log:
        "logs/deleterious_jdfe/infer_2pop_dm_warm_up.{species}.{pair}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} \
            --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --optimizations {params.optimizations} 2> {log}
        """


rule infer_2pop_dm_fine_tune:
    input:
        fs=rules.generate_2pop_fs.output.syn_fs,
        p0=rules.infer_2pop_dm_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{ref_genome}.{demog}.InferDM.bestfits",
    params:
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{ref_genome}.{demog}",
        ubounds=dadi_config["demog_2d_ub"],
        lbounds=dadi_config["demog_2d_lb"],
        grid_size=dadi_config["dm_grid_size"],
        optimizations=dadi_config["optimizations"],
    resources:
        time=1440,
        cpus=16,
    log:
        "logs/deleterious_jdfe/infer_2pop_dm_fine_tune.{species}.{pair}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        dadi-cli InferDM --fs {input.fs} --model {wildcards.demog} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} \
            --grids {params.grid_size} --output-prefix {params.output_prefix} --cpus {resources.cpus} --force-convergence {params.optimizations} 2> {log}
        """


rule get_2pop_dm_top_10_bestfits:
    input:
        bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
    output:
        bestfit="results/dadi/{species}/jdfe/{pair}/InferDM/{pair}.{ref_genome}.{demog}.InferDM.top10.bestfits.tsv",
    log:
        "logs/deleterious_jdfe/get_2pop_dm_top_10_bestfits.{species}.{pair}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        awk '
            /^# Converged results/ {{on=1; next}}
            /^# Top 100 results/   {{exit}}
            on                     {{print; if(++n==10) exit}}
        ' {input.bestfit} > {output.bestfit} 2> {log}
        """


rule convert_2pop_dm_top_10_bestfits_html:
    input:
        tsv=rules.get_2pop_dm_top_10_bestfits.output.bestfit,
    output:
        html=report(
            "results/dadi/{species}/jdfe/{pair}/html/{pair}.{ref_genome}.{demog}.InferDM.top10.bestfits.html",
            category="Distribution of Fitness Effects",
            subcategory="Two Populations",
            labels=lambda wildcards: fitted_2pop_dm_labels(wildcards, type="Bestfit Table"),
        ),
    log:
        "logs/deleterious_jdfe/convert_2pop_dm_top_10_bestfits_html.{species}.{pair}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml",
    script:
        "../scripts/tsv2html.R"


rule generate_2d_cache:
    input:
        bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        cache="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.spectra.bpkl",
    params:
        grid_size=dadi_config["jdfe_grid_size"],
        gamma_pts=dadi_config["gamma_pts"],
    resources:
        time=1440,
        cpus=16,
        mem_gb=500,
    log:
        "logs/deleterious_jdfe/generate_2d_cache.{species}.{pair}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"

        dadi-cli GenerateCache --model {wildcards.demog}_sel --demo-popt {input.bestfit} \
            --sample-size $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * 2}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * 2}}') \
            --grids {params.grid_size} --gamma-pts {params.gamma_pts} --output {output.cache} --cpus {resources.cpus} --cache-type cache2d 2> {log}
        """


rule get_dfe_bestfits:
    input:
        cache=rules.generate_2d_cache.output.cache,
    output:
        constants="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.{dfe}.constants",
    params:
        demog_1d=dadi_config["demog_1d"],
    log:
        "logs/deleterious_jdfe/get_dfe_bestfits.{species}.{pair}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"

        awk '
        FNR==NR {{
          if ($0 ~ /Converged/) {{ getline; getline; a1=$2; a2=$3 }}
          next
        }}
        {{
          if ($0 ~ /Converged/) {{ getline; getline; b1=$2; b2=$3 }}
        }}
        END {{
          printf "%.2f %.2f %.2f %.2f\\n", a1, b1, a2, b2
        }}
        ' "results/dadi/{wildcards.species}/dfe/$pop1/InferDFE/$pop1.{wildcards.ref_genome}.{params.demog_1d}.{wildcards.dfe}.InferDFE.bestfits" \
          "results/dadi/{wildcards.species}/dfe/$pop2/InferDFE/$pop2.{wildcards.ref_genome}.{params.demog_1d}.{wildcards.dfe}.InferDFE.bestfits" \
          > {output.constants} 2> {log}
        """


rule infer_jdfe_warm_up:
    input:
        fs=rules.generate_2pop_fs.output.nonsyn_fs,
        cache2d=rules.generate_2d_cache.output.cache,
        dm_bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
        constants=rules.get_dfe_bestfits.output.constants,
    output:
        p0="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}.InferDFE.opts.0",
    params:
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}",
        p0=dadi_config["dfe_2d_p0"],
        ubounds=dadi_config["dfe_2d_ub"],
        lbounds=dadi_config["dfe_2d_lb"],
        ratio=dadi_config["ratio"],
        optimizations=dadi_config["optimizations"],
    resources:
        time=1440,
        cpus=16,
        mem_gb=500,
    log:
        "logs/deleterious_jdfe/infer_jdfe_warm_up.{species}.{pair}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        constants=$(cat {input.constants})
        dadi-cli InferDFE --fs {input.fs} --cache2d {input.cache2d} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} \
            --pdf2d biv_{wildcards.dfe} --p0 {params.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --constants $constants -1 -1 \
            --ratio {params.ratio} --cpus {resources.cpus} --optimizations {params.optimizations} 2> {log}
        """


rule infer_jdfe_fine_tune:
    input:
        fs=rules.generate_2pop_fs.output.nonsyn_fs,
        cache2d=rules.generate_2d_cache.output.cache,
        dm_bestfit=rules.infer_2pop_dm_fine_tune.output.bestfit,
        constants=rules.get_dfe_bestfits.output.constants,
        p0=rules.infer_jdfe_warm_up.output.p0,
    output:
        bestfit="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}.InferDFE.bestfits",
    params:
        output_prefix="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}",
        ubounds=dadi_config["dfe_2d_ub"],
        lbounds=dadi_config["dfe_2d_lb"],
        ratio=dadi_config["ratio"],
        optimizations=dadi_config["optimizations"],
    resources:
        time=1440,
        cpus=16,
        mem_gb=500,
    log:
        "logs/deleterious_jdfe/infer_jdfe_fine_tune.{species}.{pair}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        constants=$(cat {input.constants})
        dadi-cli InferDFE --fs {input.fs} --cache2d {input.cache2d} --demo-popt {input.dm_bestfit} --output-prefix {params.output_prefix} \
            --pdf2d biv_{wildcards.dfe} --bestfit-p0-file {input.p0} --ubounds {params.ubounds} --lbounds {params.lbounds} --constants $constants -1 -1 \
            --ratio {params.ratio} --cpus {resources.cpus} --force-convergence {params.optimizations} 2> {log}
        """


rule get_2pop_jdfe_top_10_bestfits:
    input:
        bestfit=rules.infer_jdfe_fine_tune.output.bestfit,
    output:
        bestfit="results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}.InferDFE.top10.bestfits.tsv",
    log:
        "logs/deleterious_dfe/get_2pop_jdfe_top_10_bestfits.{species}.{pair}.{ref_genome}.{demog}.biv_{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        awk '
            /^# Converged results/ {{on=1; next}}
            /^# Top 100 results/   {{exit}}
            on                     {{print; if(++n==10) exit}}
        ' {input.bestfit} > {output.bestfit} 2> {log}
        """


rule convert_2pop_jdfe_top_10_bestfits_html:
    input:
        tsv=rules.get_2pop_jdfe_top_10_bestfits.output.bestfit,
    output:
        html=report(
            "results/dadi/{species}/jdfe/{pair}/html/{pair}.{ref_genome}.{demog}.biv_{dfe}.InferDFE.top10.bestfits.html",
            category="Distribution of Fitness Effects",
            subcategory="Two Populations",
            labels=lambda wildcards: fitted_jdfe_labels(wildcards, type="Bestfit Table"),
        ),
    log:
        "logs/deleterious_dfe/convert_2pop_jdfe_top_10_bestfits_html.{species}.{pair}.{ref_genome}.{demog}.biv_{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    script:
        "../scripts/tsv2html.R"


rule jdfe_godambe_ci:
    input:
        syn_vcf=f"results/{vcf_dir}/{{species}}/2pop/{{pair}}/{{pair}}.biallelic.syn.snps.{{ref_genome}}.vcf.gz",
        nonsyn_vcf=f"results/{vcf_dir}/{{species}}/2pop/{{pair}}/{{pair}}.biallelic.nonsyn.snps.{{ref_genome}}.vcf.gz",
        pop_info=rules.create_pair_info.output.pair_info,
        nonsyn_fs=rules.generate_2pop_fs.output.nonsyn_fs,
        cache2d=rules.generate_2d_cache.output.cache,
        dfe_bestfit=rules.infer_jdfe_fine_tune.output.bestfit,
    output:
        jdfe_godambe_ci="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}.godambe.ci",
    params:
        syn_dir="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_syn",
        nonsyn_dir="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_non",
        syn_output_prefix="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_syn/{pair}.{ref_genome}.syn",
        nonsyn_output_prefix="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}_bootstrapping_non/{pair}.{ref_genome}.nonsyn",
        ploidy=main_config["ploidy"],
        bootstrap_reps=dadi_config["bootstrap_replicates"],
        chunk_size=dadi_config["chunk_size"],
        polarized_flag=f"{polarization_flag}",
        mask_singletons_flag=f"{mask_singletons_flag}",
    log:
        "logs/deleterious_jdfe/jdfe_godambe_ci.{species}.{pair}.{ref_genome}.{demog}.{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"
        
        ( [ -d {params.syn_dir} ] || mkdir -p {params.syn_dir} ) 2> {log}
        ( [ -d {params.nonsyn_dir} ] || mkdir -p {params.nonsyn_dir} ) 2>> {log}
        
        dadi-cli GenerateFs --vcf {input.syn_vcf} --pop-info {input.pop_info} --pop-ids $pop1 $pop2 \
            --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') \
            {params.polarized_flag} {params.mask_singletons_flag} --bootstrap {params.bootstrap_reps} --chunk-size {params.chunk_size} \
            --output {params.syn_output_prefix} 2>> {log}
            
        dadi-cli GenerateFs --vcf {input.nonsyn_vcf} --pop-info {input.pop_info} --pop-ids $pop1 $pop2 \
            --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') \
            {params.polarized_flag} {params.mask_singletons_flag} --bootstrap {params.bootstrap_reps} --chunk-size {params.chunk_size} \
            --output {params.nonsyn_output_prefix} 2>> {log}
            
        dadi-cli StatDFE --fs {input.nonsyn_fs} --dfe-popt {input.dfe_bestfit} --cache2d {input.cache2d} \
            --pdf2d biv_{wildcards.dfe} --bootstrapping-nonsynonymous-dir {params.nonsyn_dir} --bootstrapping-synonymous-dir {params.syn_dir} \
            --output {output.jdfe_godambe_ci} 2>> {log}
        """


rule parse_jdfe_godambe_ci_table:
    input:
        bestfit=rules.infer_jdfe_fine_tune.output.bestfit,
        ci=rules.jdfe_godambe_ci.output.jdfe_godambe_ci,
    output:
        tsv="results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}.{ref_genome}.{demog}.biv_{dfe}.godambe.ci.tsv",
    log:
        "logs/deleterious_jdfe/parse_jdfe_godambe_ci_table.{species}.{pair}.{ref_genome}.{demog}.biv_{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
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


rule jdfe_ci_table_html:
    input:
        tsv=rules.parse_jdfe_godambe_ci_table.output.tsv,
    output:
        html=report(
            "results/dadi/{species}/jdfe/{pair}/html/{pair}.{ref_genome}.{demog}.biv_{dfe}.godambe.ci.html",
            category="Distribution of Fitness Effects",
            subcategory="Two Populations",
            labels=lambda wildcards: fitted_jdfe_labels(wildcards, type="Estimated 95% Uncerts"),
        ),
    log:
        "logs/reports/jdfe_ci_table_html.{species}.{pair}.{ref_genome}.{demog}.biv_{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    script:
        "../scripts/tsv2html.R"


rule plot_fitted_2pop_dm:
    input:
        fs=rules.generate_2pop_fs.output.syn_fs,
        dm_popt=rules.infer_2pop_dm_fine_tune.output.bestfit,
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        fs_plot=report(
            "results/dadi/{species}/jdfe/{pair}/plots/{pair}.{ref_genome}.{demog}.fitted.png",
            category="Distribution of Fitness Effects",
            subcategory="Two Populations",
            labels=fitted_2pop_dm_labels,
        ),
    params:
        ploidy=main_config["ploidy"],
    log:
        "logs/deleterious_jdfe/plot_fitted_2pop_dm.{species}.{pair}.{ref_genome}.{demog}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"
        
        dadi-cli Plot --fs {input.fs} --demo-popt {input.dm_popt} --model {wildcards.demog} \
            --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') \
            --output {output.fs_plot} 2> {log}
        """


rule plot_fitted_jdfe:
    input:
        fs=rules.generate_2pop_fs.output.nonsyn_fs,
        jdfe_popt=rules.infer_jdfe_fine_tune.output.bestfit,
        cache2d=rules.generate_2d_cache.output.cache,
        pop_info=rules.create_pair_info.output.pair_info,
    output:
        fs_plot=report(
            "results/dadi/{species}/jdfe/{pair}/plots/{pair}.{ref_genome}.{demog}.biv_{dfe}.fitted.png",
            category="Distribution of Fitness Effects",
            subcategory="Two Populations",
            labels=fitted_jdfe_labels,
        ),
    params:
        ploidy=main_config["ploidy"],
    log:
        "logs/deleterious_jdfe/plot_fitted_jdfe.{species}.{pair}.{ref_genome}.{demog}.biv_{dfe}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"
        
        dadi-cli Plot --fs {input.fs} --cache2d {input.cache2d} --dfe-popt {input.jdfe_popt} \
            --pdf2d biv_{wildcards.dfe} \
            --projections $(grep -w $pop1 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') $(grep -w $pop2 {input.pop_info} | awk 'END {{print NR * {params.ploidy}}}') \
            --output {output.fs_plot} 2> {log}
        """
