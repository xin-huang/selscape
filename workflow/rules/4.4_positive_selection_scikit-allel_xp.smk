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


def _add_delta_tajima_d_title(wildcards, input):
    window = int(wildcards.window)
    step = int(float(wildcards.step) * window)
    cutoff_pct = float(wildcards.cutoff) * 100

    if hasattr(input, "scores"):
        return " ".join([
            f"{wildcards.pair}",
            f"(Window={window} SNPs,",
            f"Step={step} SNPs,",
            f"Top {cutoff_pct:.2f}%)",
        ])

    return " ".join([
        f"{wildcards.pair}",
        f"Delta Tajima's D",
        f"(Window={window} SNPs,",
        f"Step={step} SNPs,",
        f"Top {cutoff_pct:.2f}%)",
    ])


rule calc_delta_tajima_d:
    input:
        vcf=rules.polarize_2pop.output.vcf,
        pair_info=rules.create_pair_info.output.pair_info,
    output:
        scores=temp(
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/chr{i}.{method}.scores.txt"
        ),
    params:
        window_size="{window}",
        step_size_ratio="{step}",
    log:
        "logs/positive_selection/calc_delta_tajima_d.{species}.{pair}.{method}.{window}_{step}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/calc_delta_tajima_d.py"


rule format_delta_tajima_d:
    input:
        scores=rules.calc_delta_tajima_d.output.scores,
    output:
        formatted=temp(
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/chr{i}.{method}.formatted.txt"
        ),
    params:
        chrom="{i}",
    log:
        "logs/positive_selection/format_delta_tajima_d.{species}.{pair}.{method}.{window}_{step}.chr{i}.log",
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
        scores=expand(
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/chr{i}.{method}.formatted.txt",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        merged_scores="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.merged.scores",
    log:
        "logs/positive_selection/merge_delta_tajima_d_scores.{species}.{pair}.{method}.{window}_{step}.log",
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
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.scores.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(
                wildcards, type="Manhattan Plot"
            ),
        ),
        candidates="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.candidates.scores",
    params:
        title=_add_delta_tajima_d_title,
        score_column="delta_tajima_d",
        cutoff="{cutoff}",
        use_absolute="TRUE",
        width=scikit_allel_config["manhattan_plot_width"],
        height=scikit_allel_config["manhattan_plot_height"],
        color1=scikit_allel_config["manhattan_plot_color1"],
        color2=scikit_allel_config["manhattan_plot_color2"],
    log:
        "logs/positive_selection/plot_delta_tajima_d.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/manhattan.R"


rule extract_delta_tajima_d_candidate_variants:
    input:
        scores=rules.plot_delta_tajima_d.output.candidates,
        vcfs=expand(
            "results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz",
            species=main_config["species"],
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        regions=temp("results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.candidates.bed"),
        variants=temp("results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.candidates.variants"),
    log:
        "logs/positive_selection/extract_delta_tajima_d_candidate_variants.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        r"""
        ( sed '1d' {input.scores} | awk '{{print "chr"$2"\t"$5"\t"$6}}' > {output.regions} ) 2> {log}
    
        for i in {input.vcfs}; do
            bcftools view -H -R {output.regions} $i | awk '{{print $1"\t"$2}}'
        done | sort -u | sed '1iCHR\tBP' > {output.variants} 2>> {log}
        """


rule annotate_delta_tajima_d_candidates:
    input:
        outliers=rules.extract_delta_tajima_d_candidate_variants.output.variants,
        annotation=expand(
            "results/annotated_data/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
            species=main_config["species"],
            i=main_config["chromosomes"],
            ref_genome=main_config["ref_genome"],
            allow_missing=True,
        ),
    output:
        annotated_candidates="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.annotated.candidates",
    resources:
        mem_gb=32,
    log:
        "logs/positive_selection/annotate_delta_tajima_d_candidates.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/get_annotated_candidates.py"


rule get_delta_tajima_d_candidate_genes:
    input:
        delta_candidates=rules.annotate_delta_tajima_d_candidates.output.annotated_candidates,
    output:
        delta_genes="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.candidate.genes",
    log:
        "logs/positive_selection/get_delta_tajima_d_candidate_genes.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.delta_candidates} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.delta_genes} ) 2> {log} || true
        sed -i '1iGene' {output.delta_genes} 2>> {log}
        """


rule delta_tajima_d_candidate_genes_table_html:
    input:
        tsv=rules.get_delta_tajima_d_candidate_genes.output.delta_genes,
    output:
        html=report(
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.candidate.genes.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(
                wildcards, type="Gene List"
            ),
        ),
    params:
        title=_add_delta_tajima_d_title,
    log:
        "logs/positive_selection/delta_tajima_d_candidate_genes_table_html.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log", 
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule enrichment_delta_tajima_d_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        candidates=rules.annotate_delta_tajima_d_candidates.output.annotated_candidates,
        total=expand(
            "results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz",
            species=main_config["species"],
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        candidate_snps="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.candidate.snps.tsv",
        total_snps="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.total.snps.tsv",
        enrichment="results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/positive_selection/enrichment_delta_tajima_d_gowinda.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        r"""
        sed '1d' {input.candidates} | awk '{{print "chr"$1"\t"$2}}' > {output.candidate_snps} 2> {log}

        for i in {input.total}; do
            bcftools query -f "%CHROM\t%POS\n" $i
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

        sed -i '1iGO_ID\tavg_genes_sim\tgenes_found\tp_value\tp_adjusted\tgenes_uniq\tgenes_max\tgenes_total\tdescription\tgene_list' {output.enrichment} 2>> {log}
        """


rule delta_tajima_d_enrichment_results_table_html:
    input:
        tsv=rules.enrichment_delta_tajima_d_gowinda.output.enrichment,
    output:
        html=report(
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.gowinda.enrichment.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(
                wildcards, type="Enrichment Table"
            ),
        ),
    params:
        title=_add_delta_tajima_d_title,
    log:
        "logs/positive_selection/delta_tajima_d_enrichment_results_table_html.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule plot_gowinda_enrichment_delta_tajima_d:
    input:
        enrichment=rules.enrichment_delta_tajima_d_gowinda.output.enrichment,
    output:
        plot=report(
            "results/positive_selection/scikit-allel/{species}/2pop/{pair}/{method}/{window}_{step}/{pair}.{method}.top_{cutoff}.gowinda.enrichment.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: delta_tajima_d_labels(wildcards, type="Enrichment Plot"),
        ),
    params:
        title=_add_delta_tajima_d_title,
    resources:
        mem_gb=8,
    log:
        "logs/positive_selection/plot_gowinda_enrichment_delta_tajima_d.{species}.{pair}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/plot_gowinda_enrichment.py"
