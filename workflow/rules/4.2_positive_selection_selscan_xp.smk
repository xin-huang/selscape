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


rule extract_pair_snps:
    input:
        vcf=rules.polarize_2pop.output.vcf,
        pair_info=rules.create_pair_info.output.pair_info,
    output:
        vcf1=temp(
            "results/polarized_data/{species}/{dataset}/2pop/{pair}/pop1.chr{i}.biallelic.snps.vcf.gz"
        ),
        vcf2=temp(
            "results/polarized_data/{species}/{dataset}/2pop/{pair}/pop2.chr{i}.biallelic.snps.vcf.gz"
        ),
        idx1=temp(
            "results/polarized_data/{species}/{dataset}/2pop/{pair}/pop1.chr{i}.biallelic.snps.vcf.gz.tbi"
        ),
        idx2=temp(
            "results/polarized_data/{species}/{dataset}/2pop/{pair}/pop2.chr{i}.biallelic.snps.vcf.gz.tbi"
        ),
        map=temp(
            "results/polarized_data/{species}/{dataset}/2pop/{pair}/chr{i}.biallelic.snps.map"
        ),
    log:
        "logs/positive_selection/extract_pair_snps.{species}.{dataset}.{pair}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        pair="{wildcards.pair}"
        pop1="${{pair%%_*}}"
        pop2="${{pair##*_}}"
        ( bcftools view {input.vcf} -S <(grep -w $pop1 {input.pair_info} | awk '{{print $1}}') --force-samples | bgzip -c > {output.vcf1} ) 2> {log}
        ( bcftools view {input.vcf} -S <(grep -w $pop2 {input.pair_info} | awk '{{print $1}}') --force-samples | bgzip -c > {output.vcf2} ) 2>> {log}
        ( bcftools query -f "%CHROM\\t%CHROM:%POS:%REF:%ALT\\t%POS\\t%POS\\n" {input.vcf} > {output.map} ) 2>> {log}
        tabix -p vcf {output.vcf1} 2>> {log}
        tabix -p vcf {output.vcf2} 2>> {log}
        """


rule estimate_selscan_xp_scores:
    input:
        selscan=rules.download_selscan.output.selscan,
        vcf1=rules.extract_pair_snps.output.vcf1,
        vcf2=rules.extract_pair_snps.output.vcf2,
        map=rules.extract_pair_snps.output.map,
    output:
        out=temp(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.chr{i}.{method}.out"
        ),
        formatted_out=temp(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.chr{i}.{method}.formatted.out"
        ),
        log=temp(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.chr{i}.{method}.log"
        ),
    params:
        output_prefix="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.chr{i}",
        phasing_flag=get_phasing_flag,
    resources:
        cpus=8,
    log:
        "logs/positive_selection/estimate_selscan_xp_scores.{species}.{dataset}.{pair}.{method}.{maf}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        {input.selscan} --vcf {input.vcf1} --vcf-ref {input.vcf2} --map {input.map} --{wildcards.method} --out {params.output_prefix} --threads {resources.cpus} --maf {wildcards.maf} {params.phasing_flag} 2> {log}
        head -1 {output.out} > {output.formatted_out} 2>> {log}
        ( sed '1d' {output.out} | awk -v chr={wildcards.i} 'BEGIN{{OFS="\\t"}}{{print chr,$2,$3,$4,$5,$6,$7,$8}}' >> {output.formatted_out} ) 2>> {log}
        """


rule normalize_selscan_xp_scores:
    input:
        norm=rules.download_selscan.output.norm,
        scores=lambda wc: expand(
            rules.estimate_selscan_xp_scores.output.formatted_out,
            i=get_chromosomes(wc),
            allow_missing=True,
        ),
    output:
        done=temp(touch(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.{method}.normalized.done"
        )),
        log=temp(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.{method}.log"
        ),
    log:
        "logs/positive_selection/normalize_selscan_xp_scores.{species}.{dataset}.{pair}.{method}.{maf}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        {input.norm} --files {input.scores} --log {output.log} --{wildcards.method} 2> {log}
        """


rule merge_selscan_xp_scores:
    input:
        done=rules.normalize_selscan_xp_scores.output.done,
    params:
        scores=lambda wc: expand(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.chr{i}.{method}.formatted.out.norm",
            i=get_chromosomes(wc),
            species=wc.species,
            dataset=wc.dataset,
            pair=wc.pair,
            method=wc.method,
            maf=wc.maf,
            allow_missing=True,
        ),
    output:
        merged_scores="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.scores",
    log:
        "logs/positive_selection/merge_selscan_xp_scores.{species}.{dataset}.{pair}.{method}.{maf}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( cat {params.scores} | grep -v id | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$9}}' | sed '1iSNP\\tCHR\\tBP\\tnormalized_{wildcards.method}' > {output.merged_scores} ) 2> {log}
        """


rule plot_selscan_xp:
    input:
        scores=rules.merge_selscan_xp_scores.output.merged_scores,
    output:
        outliers="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.outliers.scores",
        plot=report(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.scores.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=selscan_xp_labels,
        ),
    params:
        title=add_selscan_title,
        score_column="normalized_{method}",
        use_absolute="TRUE",
        cutoff="{cutoff}",
        width=selscan_config["manhattan_plot_width"],
        height=selscan_config["manhattan_plot_height"],
        color1=selscan_config["manhattan_plot_color1"],
        color2=selscan_config["manhattan_plot_color2"],
    resources:
        mem_mb=16000,
    log:
        "logs/positive_selection/plot_selscan_xp.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/manhattan.R"


rule annotate_selscan_xp_outliers:
    input:
        outliers="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.outliers.scores",
        annotation=lambda wc: expand(
            rules.annotate_biallelic_snps.output.txt,
            i=get_chromosomes(wc),
            ref_genome=get_ref_genome(wc),
            allow_missing=True,
        ),
    output:
        annotated_outliers="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.annotated.outliers",
    resources:
        mem_mb=32000,
    log:
        "logs/positive_selection/annotate_selscan_xp_outliers.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/functional_enrichment/get_annotated_outliers.py"


rule get_selscan_xp_outlier_genes:
    input:
        selscan_xp_outliers=rules.annotate_selscan_xp_outliers.output.annotated_outliers,
    output:
        selscan_xp_genes="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.outlier.genes",
    log:
        "logs/positive_selection/get_selscan_xp_outlier_genes.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.selscan_xp_outliers} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.selscan_xp_genes} ) 2> {log} || true
        sed -i '1iGene' {output.selscan_xp_genes} 2>> {log}
        """


rule selscan_xp_outlier_genes_table_html:
    input:
        tsv="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.outlier.genes",
    output:
        html=report(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.outlier.genes.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: selscan_xp_labels(wildcards, type="Gene List"),
        ),
    params:
        title=add_selscan_title,
    log:
        "logs/positive_selection/selscan_xp_outlier_genes_table_html.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule enrichment_selscan_xp_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        outliers=rules.annotate_selscan_xp_outliers.output.annotated_outliers,
        total=lambda wc: expand(
            rules.polarize_2pop.output.vcf,
            i=get_chromosomes(wc),
            allow_missing=True,
        ),
    output:
        outlier_snps="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.outlier.snps.tsv",
        total_snps="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.total.snps.tsv",
        enrichment="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_mb=32000,
        cpus=8,
    log:
        "logs/positive_selection/enrichment_selscan_xp_gowinda.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        sed '1d' {input.outliers} | awk '{{print "chr"$1"\\t"$2}}' > {output.outlier_snps} 2> {log}
        for i in {input.total}; do
            bcftools query -f "%CHROM\\t%POS\\n" $i
        done > {output.total_snps} 2>> {log}
        java -Xmx{resources.mem_mb}m -jar {input.gowinda} \
            --snp-file {output.total_snps} \
            --candidate-snp-file {output.outlier_snps} \
            --gene-set-file {input.go2gene} \
            --annotation-file {input.gtf} \
            --simulations 1000000 \
            --min-significance 1 \
            --gene-definition gene \
            --threads {resources.cpus} \
            --output-file {output.enrichment} \
            --mode gene \
            --min-genes 1 >> {log} 2>&1 || true
        sed -i '1iGO_ID\\tavg_genes_sim\\tgenes_found\\tp_value\\tp_adjusted\\tgenes_uniq\\tgenes_max\\tgenes_total\\tdescription\\tgene_list' {output.enrichment}
        """


rule selscan_xp_enrichment_results_table_html:
    input:
        tsv="results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.tsv",
    output:
        html=report(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: selscan_xp_labels(
                wildcards, type="Enrichment Table"),
        ),
    params:
        title=add_selscan_title,
    log:
        "logs/positive_selection/selscan_xp_enrichment_results_table_html.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule plot_gowinda_enrichment_selscan_xp:
    input:
        enrichment=rules.enrichment_selscan_xp_gowinda.output.enrichment,
    output:
        plot=report(
            "results/positive_selection/selscan/{species}/{dataset}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: selscan_xp_labels(wildcards, type="Enrichment Plot"),
        ),
    params:
        title=add_selscan_title,
    resources:
        mem_mb=8000,
    log:
        "logs/positive_selection/plot_gowinda_enrichment_selscan_xp.{species}.{dataset}.{pair}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_gowinda_enrichment.py"
