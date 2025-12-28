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


vcf_dir = "processed_data"
folding_flag = ""
if (
    "anc_alleles" in main_config
    and main_config["anc_alleles"]
    and betascan_config["unfolded"]
):
    vcf_dir = "polarized_data"
    folding_flag = "-fold"


rule get_allele_counts:
    input:
        vcf=f"results/{vcf_dir}/{{species}}/1pop/{{ppl}}/{{ppl}}.chr{{i}}.biallelic.snps.repeats.removed.vcf.gz",
    output:
        ac=temp(
            "results/balancing_selection/betascan/{species}/{ppl}/ac/{ppl}.{ref_genome}.chr{i}.ac"
        ),
    params:
        ploidy=main_config["ploidy"],
        min_af=betascan_config["min_af"],
        max_af=betascan_config["max_af"],
    log:
        "logs/balancing_selection/get_allele_counts.{species}.{ppl}.{ref_genome}.chr{i}.log",
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
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.chr{i}.b1.scores"
        ),
    params:
        folding_flag=f"{folding_flag}",
    log:
        "logs/balancing_selection/estimate_b1_scores.{species}.{ppl}.{ref_genome}.m_{core_frq}.chr{i}.log",
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
        scores=expand(
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.chr{i}.b1.scores",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        merged_scores="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.scores",
    log:
        "logs/balancing_selection/merge_b1_scores.{species}.{ppl}.m_{core_frq}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( cat {input.scores} | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$3}}' | sed '1iSNP\\tCHR\\tBP\\tB1' > {output.merged_scores} ) 2> {log}
        """


rule annotate_betascan_candidates:
    input:
        outliers="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.candidates.scores",
        annotation=expand(
            "results/annotated_data/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
            i=main_config["chromosomes"],
            ref_genome=main_config["ref_genome"],
            allow_missing=True,
        ),
    output:
        annotated_candidates="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.annotated.candidates",
    resources:
        mem_gb=32,
    log:
        "logs/balancing_selection/annotate_betascan_candidates.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/get_annotated_candidates.py"


rule get_betascan_candidate_genes:
    input:
        betascan_candidates=rules.annotate_betascan_candidates.output.annotated_candidates,
    output:
        betascan_genes="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.candidate.genes",
    log:
        "logs/balancing_selection/get_betascan_candidate_genes.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.betascan_candidates} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.betascan_genes} ) 2> {log} || true
        sed -i '1iGene' {output.betascan_genes} 2>> {log}
        """


rule enrichment_betascan_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        candidates=rules.annotate_betascan_candidates.output.annotated_candidates,
        total=expand(
            f"results/{vcf_dir}/{{species}}/1pop/{{ppl}}/{{ppl}}.chr{{i}}.biallelic.snps.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        candidate_snps="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.candidate.snps.tsv",
        total_snps="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.total.snps.tsv",
        enrichment="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/balancing_selection/enrichment_betascan_gowinda.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.candidates} | awk '{{print "chr"$1"\\t"$2}}' > {output.candidate_snps} ) 2> {log}

        for i in {input.total}; do
            bcftools query -f "%CHROM\\t%POS\\n" $i
        done > {output.total_snps} 2>> {log}

        java -Xmx{resources.mem_gb}g -jar {input.gowinda} \
            --snp-file {output.total_snps} \
            --candidate-snp-file {output.candidate_snps} \
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


rule betascan_candidate_genes_table_html:
    input:
        tsv=rules.get_betascan_candidate_genes.output.betascan_genes,
    output:
        html=report(
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.candidate.genes.html",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(wildcards, type="Gene List"),
        ),
    params:
        title=lambda w: f"{w.ppl} B1 (Core Freq={w.core_frq}, Top {float(w.cutoff)*100:.2f}%) CANDIDATE GENES",
    log:
        "logs/balancing_selection/betascan_candidate_genes_table_html.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule betascan_enrichment_results_table_html:
    input:
        tsv=rules.enrichment_betascan_gowinda.output.enrichment,
    output:
        html=report(
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.enrichment.html",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(
                wildcards, type="Enrichment Table"
            ),
        ),
    params:
        title=lambda w: f"{w.ppl} B1 (Core Freq={w.core_frq}, Top {float(w.cutoff)*100:.2f}%) ENRICHMENT",
    log:
        "logs/balancing_selection/betascan_enrichment_results_table_html.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule plot_betascan:
    input:
        scores=rules.merge_b1_scores.output.merged_scores,
    output:
        candidates="results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.candidates.scores",
        plot=report(
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.scores.png",
            category="Balancing Selection",
            subcategory="B1",
            labels=betascan_labels,
        ),
    params:
        title=lambda w: f"{w.ppl} (Core Freq={w.core_frq}, Top {float(w.cutoff)*100:.2f}%)",
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
        "logs/balancing_selection/plot_betascan.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/manhattan.R"


rule plot_gowinda_enrichment_betascan:
    input:
        enrichment=rules.enrichment_betascan_gowinda.output.enrichment,
    output:
        count_plot=report(
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.enrichment.png",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(
                wildcards, type="Enrichment Plot"
            ),
        ),
        qscore_plot=report(
            "results/balancing_selection/betascan/{species}/{ppl}/m_{core_frq}/{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.gowinda.qscore.png",
            category="Balancing Selection",
            subcategory="B1",
            labels=lambda wildcards: betascan_labels(wildcards, type="Q-Score Plot"),
        ),
    params:
        title=lambda w: f"{w.ppl} B1 ENRICHMENT (Core Freq={w.core_frq}, Top {float(w.cutoff)*100:.2f}%)",
    resources:
        mem_gb=8,
    log:
        "logs/balancing_selection/plot_gowinda_enrichment_betascan.{species}.{ppl}.{ref_genome}.m_{core_frq}.b1.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/plot_gowinda_enrichment.py"
