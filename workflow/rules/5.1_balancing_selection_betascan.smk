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



rule get_allele_counts:
    input:
        vcf=lambda wc: f"results/{get_betascan_vcf_dir(wc)}/{wc.dataset}/{wc.species}/1pop/{wc.ppl}/{wc.ppl}.chr{wc.i}.biallelic.snps.repeats.removed.vcf.gz",
    output:
        ac=temp(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/ac/{ppl}.{ref_genome}.chr{i}.ac"
        ),
    params:
        ploidy=get_ploidy,
        min_af=betascan_config["min_af"],
        max_af=betascan_config["max_af"],
    log:
        "logs/balancing_selection/get_allele_counts.{dataset}.{species}.{ppl}.{ref_genome}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        bcftools view -i '(INFO/AC>{params.ploidy}*N_SAMPLES*{params.min_af}) && (INFO/AC<{params.ploidy}*N_SAMPLES*{params.max_af})' {input.vcf} \
            | bcftools query -f '%POS\t%INFO/AC\t%INFO/AN\n' > {output.ac} 2>> {log}
        """


rule estimate_b1_scores:
    input:
        ac=rules.get_allele_counts.output.ac,
        betascan=rules.download_betascan.output.betascan,
    output:
        scores=temp(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.chr{i}.b1.scores"
        ),
    params:
        folding_flag=get_folding_flag,
    log:
        "logs/balancing_selection/estimate_b1_scores.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( python {input.betascan} -i {input.ac} -m {wildcards.core_frq} {params.folding_flag} | \
        grep -v Position | \
        awk -v chr={wildcards.i} '{{print chr"\\t"$0}}' > {output.scores} ) 2> {log}
        """


rule merge_b1_scores:
    input:
        scores=lambda wc: expand(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.chr{i}.b1.scores",
            dataset=wc.dataset, species=wc.species, ppl=wc.ppl,
            core_frq=wc.core_frq, ref_genome=wc.ref_genome,
            i=get_chromosomes(wc),
        ),
    output:
        merged_scores="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.scores",
    log:
        "logs/balancing_selection/merge_b1_scores.{dataset}.{species}.{ppl}.m_{core_frq}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( cat {input.scores} | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$3}}' | sed '1iSNP\\tCHR\\tBP\\tB1' > {output.merged_scores} ) 2> {log}
        """


rule annotate_betascan_outliers:
    input:
        outliers="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outliers.scores",
        annotation=lambda wc: expand(
            "results/annotated_data/{dataset}/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
            dataset=wc.dataset, species=wc.species,
            i=get_chromosomes(wc),
            ref_genome=get_ref_genome(wc),
        ),
    output:
        annotated_outliers="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.annotated.outliers",
    resources:
        mem_gb=32,
    log:
        "logs/balancing_selection/annotate_betascan_outliers.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/functional_enrichment/get_annotated_outliers.py"


rule get_betascan_outlier_genes:
    input:
        betascan_outliers=rules.annotate_betascan_outliers.output.annotated_outliers,
    output:
        betascan_genes="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outlier.genes",
    log:
        "logs/balancing_selection/get_betascan_outlier_genes.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.betascan_outliers} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.betascan_genes} ) 2> {log} || true
        sed -i '1iGene' {output.betascan_genes} 2>> {log}
        """


rule enrichment_betascan_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        outliers=rules.annotate_betascan_outliers.output.annotated_outliers,
        total=lambda wc: expand(
            f"results/{get_betascan_vcf_dir(wc)}/{{dataset}}/{{species}}/1pop/{{ppl}}/{{ppl}}.chr{{i}}.biallelic.snps.vcf.gz",
            dataset=wc.dataset, species=wc.species, ppl=wc.ppl,
            i=get_chromosomes(wc),
        ),
    output:
        outlier_snps="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outlier.snps.tsv",
        total_snps="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.total.snps.tsv",
        enrichment="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/balancing_selection/enrichment_betascan_gowinda.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.outliers} | awk '{{print "chr"$1"\\t"$2}}' > {output.outlier_snps} ) 2> {log}

        for i in {input.total}; do
            bcftools query -f "%CHROM\\t%POS\\n" $i
        done > {output.total_snps} 2>> {log}

        java -Xmx{resources.mem_gb}g -jar {input.gowinda} \
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

        sed -i '1iGO_ID\\tavg_genes_sim\\tgenes_found\\tp_value\\tp_adjusted\\tgenes_uniq\\tgenes_max\\tgenes_total\\tdescription\\tgene_list' {output.enrichment} 2>> {log}
        """


rule betascan_outlier_genes_table_html:
    input:
        tsv=rules.get_betascan_outlier_genes.output.betascan_genes,
    output:
        html=report(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outlier.genes.html",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(wildcards, type="Gene List"),
        ),
    params:
        title=add_betascan_title,
    log:
        "logs/balancing_selection/betascan_outlier_genes_table_html.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule betascan_enrichment_results_table_html:
    input:
        tsv=rules.enrichment_betascan_gowinda.output.enrichment,
    output:
        html=report(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.enrichment.html",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(
                wildcards, type="Enrichment Table"
            ),
        ),
    params:
        title=add_betascan_title,
    log:
        "logs/balancing_selection/betascan_enrichment_results_table_html.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule plot_betascan:
    input:
        scores=rules.merge_b1_scores.output.merged_scores,
    output:
        outliers="results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.outliers.scores",
        plot=report(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.scores.png",
            category="Balancing Selection",
            subcategory="B1",
            labels=betascan_labels,
        ),
    params:
        title=add_betascan_title,
        score_column="B1",
        use_absolute="FALSE",
        cutoff="{cutoff}",
        width=betascan_config["manhattan_plot_width"],
        height=betascan_config["manhattan_plot_height"],
        color1=betascan_config["manhattan_plot_color1"],
        color2=betascan_config["manhattan_plot_color2"],
    resources:
        mem_gb=16,
    log:
        "logs/balancing_selection/plot_betascan.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/manhattan.R"


rule plot_gowinda_enrichment_betascan:
    input:
        enrichment=rules.enrichment_betascan_gowinda.output.enrichment,
    output:
        plot=report(
            "results/balancing_selection/betascan/{dataset}/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.enrichment.png",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(wildcards, type="Enrichment Plot"),
        ),
    params:
        title=add_betascan_title,
    resources:
        mem_gb=8,
    log:
        "logs/balancing_selection/plot_gowinda_enrichment_betascan.{dataset}.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_gowinda_enrichment.py"
