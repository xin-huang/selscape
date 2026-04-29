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


rule calc_delta_tajima_d:
    input:
        vcf=rules.polarize_2pop.output.vcf,
        pair_info=rules.create_pair_info.output.pair_info,
    output:
        scores=temp(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/chr{i}.{method}.scores.txt"
        ),
    params:
        window_size="{window}",
        step_size_ratio="{step}",
    log:
        "logs/positive_selection/calc_delta_tajima_d.{dataset}.{species}.{pair}.{method}.{window}_{step}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/selection_analysis/calc_delta_tajima_d.py"


rule format_delta_tajima_d:
    input:
        scores=rules.calc_delta_tajima_d.output.scores,
    output:
        formatted=temp(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/chr{i}.{method}.formatted.txt"
        ),
    params:
        chrom="{i}",
    log:
        "logs/positive_selection/format_delta_tajima_d.{dataset}.{species}.{pair}.{method}.{window}_{step}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        awk -v chr="{params.chrom}" 'BEGIN{{OFS="\\t"}}
             NR==1{{print "SNP", "CHR", "BP", "delta_tajima_d", "window_start", "window_end", "n_snps"}}
             NR>1 {{print chr":"$1, chr, $1, $4, $1, $2, $3}}' \
        {input.scores} > {output.formatted} 2> {log}
        """


rule merge_delta_tajima_d_scores:
    input:
        scores=lambda wc: expand(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/chr{i}.{method}.formatted.txt",
            dataset=wc.dataset, species=wc.species, pair=wc.pair,
            method=wc.method, window=wc.window, step=wc.step,
            i=get_chromosomes(wc),
        ),
    output:
        merged_scores="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.merged.scores",
    log:
        "logs/positive_selection/merge_delta_tajima_d_scores.{dataset}.{species}.{pair}.{method}.{window}_{step}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        head -1 {input.scores[0]} > {output.merged_scores} 2> {log}
        for file in {input.scores}; do
            tail -n +2 $file >> {output.merged_scores}
        done 2>> {log}
        """


rule plot_delta_tajima_d:
    input:
        scores=rules.merge_delta_tajima_d_scores.output.merged_scores,
    output:
        plot=report(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.scores.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(
                wildcards, type="Manhattan Plot"
            ),
        ),
        outliers="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.outliers.scores",
    params:
        title=add_scikit_allel_title,
        score_column="delta_tajima_d",
        cutoff="{cutoff}",
        use_absolute="TRUE",
        width=scikit_allel_config["manhattan_plot_width"],
        height=scikit_allel_config["manhattan_plot_height"],
        color1=scikit_allel_config["manhattan_plot_color1"],
        color2=scikit_allel_config["manhattan_plot_color2"],
    log:
        "logs/positive_selection/plot_delta_tajima_d.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/manhattan.R"


rule extract_delta_tajima_d_outlier_variants:
    input:
        scores=rules.plot_delta_tajima_d.output.outliers,
        vcfs=lambda wc: expand(
            "results/processed_data/{dataset}/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz",
            dataset=wc.dataset, species=wc.species, pair=wc.pair,
            i=get_chromosomes(wc),
        ),
    output:
        regions=temp("results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.outliers.bed"),
        variants=temp("results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.outliers.variants"),
    log:
        "logs/positive_selection/extract_delta_tajima_d_outlier_variants.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        r"""
        ( sed '1d' {input.scores} | awk '{{print "chr"$2"\t"$5"\t"$6}}' > {output.regions} ) 2> {log}
    
        for i in {input.vcfs}; do
            bcftools view -H -R {output.regions} $i | awk '{{print $1"\t"$2}}'
        done | sort -u | sed '1iCHR\tBP' > {output.variants} 2>> {log}
        """


rule annotate_delta_tajima_d_outliers:
    input:
        outliers=rules.extract_delta_tajima_d_outlier_variants.output.variants,
        annotation=lambda wc: expand(
            "results/annotated_data/{dataset}/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
            dataset=wc.dataset, species=wc.species,
            i=get_chromosomes(wc),
            ref_genome=get_ref_genome(wc),
        ),
    output:
        annotated_outliers="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.annotated.outliers",
    resources:
        mem_gb=32,
    log:
        "logs/positive_selection/annotate_delta_tajima_d_outliers.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/functional_enrichment/get_annotated_outliers.py"


rule get_delta_tajima_d_outlier_genes:
    input:
        delta_outliers=rules.annotate_delta_tajima_d_outliers.output.annotated_outliers,
    output:
        delta_genes="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.outlier.genes",
    log:
        "logs/positive_selection/get_delta_tajima_d_outlier_genes.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.delta_outliers} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.delta_genes} ) 2> {log} || true
        sed -i '1iGene' {output.delta_genes} 2>> {log}
        """


rule delta_tajima_d_outlier_genes_table_html:
    input:
        tsv=rules.get_delta_tajima_d_outlier_genes.output.delta_genes,
    output:
        html=report(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.outlier.genes.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(
                wildcards, type="Gene List"
            ),
        ),
    params:
        title=add_scikit_allel_title,
    log:
        "logs/positive_selection/delta_tajima_d_outlier_genes_table_html.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log", 
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule enrichment_delta_tajima_d_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        outliers=rules.annotate_delta_tajima_d_outliers.output.annotated_outliers,
        total=lambda wc: expand(
            "results/processed_data/{dataset}/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz",
            dataset=wc.dataset, species=wc.species, pair=wc.pair,
            i=get_chromosomes(wc),
        ),
    output:
        outlier_snps="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.outlier.snps.tsv",
        total_snps="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.total.snps.tsv",
        enrichment="results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/positive_selection/enrichment_delta_tajima_d_gowinda.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        r"""
        sed '1d' {input.outliers} | awk '{{print "chr"$1"\t"$2}}' > {output.outlier_snps} 2> {log}

        for i in {input.total}; do
            bcftools query -f "%CHROM\t%POS\n" $i
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

        sed -i '1iGO_ID\tavg_genes_sim\tgenes_found\tp_value\tp_adjusted\tgenes_uniq\tgenes_max\tgenes_total\tdescription\tgene_list' {output.enrichment} 2>> {log}
        """


rule delta_tajima_d_enrichment_results_table_html:
    input:
        tsv=rules.enrichment_delta_tajima_d_gowinda.output.enrichment,
    output:
        html=report(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.gowinda.enrichment.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(
                wildcards, type="Enrichment Table"
            ),
        ),
    params:
        title=add_scikit_allel_title,
    log:
        "logs/positive_selection/delta_tajima_d_enrichment_results_table_html.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/tsv2html.R"


rule plot_gowinda_enrichment_delta_tajima_d:
    input:
        enrichment=rules.enrichment_delta_tajima_d_gowinda.output.enrichment,
    output:
        plot=report(
            "results/positive_selection/scikit-allel/{dataset}/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.gowinda.enrichment.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(wildcards, type="Enrichment Plot"),
        ),
    params:
        title=add_scikit_allel_title,
    resources:
        mem_gb=8,
    log:
        "logs/positive_selection/plot_gowinda_enrichment_delta_tajima_d.{dataset}.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/visualization/plot_gowinda_enrichment.py"
