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


rule calc_tajima_d_balancing:
    input:
        vcf=rules.remove_repeats.output.vcf,
    output:
        scores=temp(
            "results/balancing_selection/scikit-allel/{species}/1pop/{ppl}/{method}/{window}_{step}/chr{i}.{method}.scores.txt"
        ),
    params:
        window_size="{window}",
        step_size_ratio="{step}",
    log:
        "logs/balancing_selection/calc_tajima_d_balancing.{species}.{ppl}.{method}.{window}_{step}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/calc_tajima_d.py"


rule format_tajima_d_balancing:
    input:
        scores=rules.calc_tajima_d_balancing.output.scores,
    output:
        formatted=temp(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/chr{i}.{method}.tajima_d.formatted.txt"
        ),
    params:
        chrom="{i}",
        method="{method}",
    log:
        "logs/balancing_selection/format_tajima_d_balancing.{species}.{ppl}.{method}.{window}_{step}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        awk -v chr="{params.chrom}" 'BEGIN{{OFS="\\t"}}
            NR==1{{print "SNP", "CHR", "BP", "tajima_d", "window_start", "window_end", "n_snps"}}
            NR>1 && $4>0 {{print chr":"$1, chr, $1, $4, $1, $2, $3}}' \
            {input.scores} > {output.formatted} 2> {log}
        """


rule merge_tajima_d_balancing_scores:
    input:
        scores=expand(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/chr{i}.{method}.tajima_d.formatted.txt",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        merged_scores="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.merged.scores",
    log:
        "logs/balancing_selection/merge_tajima_d_balancing_scores.{species}.{ppl}.{method}.{window}_{step}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        head -1 {input.scores[0]} > {output.merged_scores} 2> {log}
        for file in {input.scores}; do
            tail -n +2 $file >> {output.merged_scores}
        done 2>> {log}
        """


rule plot_tajima_d_balancing:
    input:
        scores=rules.merge_tajima_d_balancing_scores.output.merged_scores,
    output:
        plot=report(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.scores.png",
            category="Balancing Selection",
            subcategory="{method}",
            labels=lambda wildcards: tajima_d_labels(wildcards, type="Manhattan Plot"),
        ),
        candidates="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.candidates.scores",
    params:
        title=lambda w: f"{w.ppl} (Window={w.window}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Step={int(float(w.step) * int(w.window))}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Top {float(w.cutoff)*100:.2f}%)",
        score_column="tajima_d",
        cutoff="{cutoff}",
        use_absolute="FALSE",
        ylab="Tajima's D",
        width=scikit_allel_config["manhattan_plot_width"],
        height=scikit_allel_config["manhattan_plot_height"],
        color1=scikit_allel_config["manhattan_plot_color1"],
        color2=scikit_allel_config["manhattan_plot_color2"],
    log:
        "logs/balancing_selection/plot_tajima_d_balancing.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/manhattan.R"


rule extract_tajima_d_balancing_candidate_variants:
    input:
        scores=rules.plot_tajima_d_balancing.output.candidates,
        vcfs=expand(
            "results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        regions=temp("results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.candidates.bed"),
        variants=temp("results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.candidates.variants"),
    log:
        "logs/balancing_selection/extract_tajima_d_balancing_candidate_variants.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        r"""
        ( sed '1d' {input.scores} | awk '{{print "chr"$2"\t"$5"\t"$6}}' > {output.regions} ) 2> {log}
    
        for i in {input.vcfs}; do
            bcftools view -H -R {output.regions} $i | awk '{{print $1"\t"$2}}'
        done | sort -u | sed '1iCHR\tBP' > {output.variants} 2>> {log}
        """


rule annotate_tajima_d_balancing_candidates:
    input:
        outliers=rules.extract_tajima_d_balancing_candidate_variants.output.variants,
        annotation=expand(
            "results/annotated_data/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
            i=main_config["chromosomes"],
            ref_genome=main_config["ref_genome"],
            allow_missing=True,
        ),
    output:
        annotated_candidates="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.annotated.candidates",
    resources:
        mem_gb=32,
    log:
        "logs/balancing_selection/annotate_tajima_d_balancing_candidates.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/get_annotated_candidates.py"


rule get_tajima_d_balancing_candidate_genes:
    input:
        tajima_d_candidates=rules.annotate_tajima_d_balancing_candidates.output.annotated_candidates,
    output:
        tajima_d_genes="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.candidate.genes",
    log:
        "logs/balancing_selection/get_tajima_d_balancing_candidate_genes.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.tajima_d_candidates} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.tajima_d_genes} ) 2> {log} || true
        sed -i '1iGene' {output.tajima_d_genes} 2>> {log}
        """


rule tajima_d_balancing_candidate_genes_table_html:
    input:
        tsv=rules.get_tajima_d_balancing_candidate_genes.output.tajima_d_genes,
    output:
        html=report(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.candidate.genes.html",
            category="Balancing Selection",
            subcategory="{method}",
            labels=lambda wildcards: tajima_d_labels(wildcards, type="Gene List"),
        ),
    params:
        title=lambda w: f"{w.ppl} {w.method.upper().replace('_', ' ')} (Window={w.window}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Step={int(float(w.step) * int(w.window))}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Top {float(w.cutoff)*100:.2f}%) CANDIDATE GENES ",
    log:
        "logs/balancing_selection/tajima_d_balancing_candidate_genes_table_html.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule enrichment_tajima_d_balancing_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        candidates=rules.annotate_tajima_d_balancing_candidates.output.annotated_candidates,
        total=expand(
            "results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        candidate_snps="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.candidate.snps.tsv",
        total_snps="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.total.snps.tsv",
        enrichment="results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/balancing_selection/enrichment_tajima_d_balancing_gowinda.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
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

        sed -i '1iGO_ID\tavg_genes_sim\tgenes_found\tp_value\tp_adjusted\tgenes_uniq\tgenes_max\tgenes_total\tdescription\tgene_list' {output.enrichment}
        """


rule tajima_d_balancing_enrichment_results_table_html:
    input:
        tsv=rules.enrichment_tajima_d_balancing_gowinda.output.enrichment,
    output:
        html=report(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.gowinda.enrichment.html",
            category="Balancing Selection",
            subcategory="{method}",
            labels=lambda wildcards: tajima_d_labels(
                wildcards, type="Enrichment Table"
            ),
        ),
    params:
        title=lambda w: f"{w.ppl} {w.method.upper().replace('_', ' ')} (Window={w.window}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Step={int(float(w.step) * int(w.window))}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Top {float(w.cutoff)*100:.2f}%) ENRICHMENT",
    log:
        "logs/balancing_selection/tajima_d_balancing_enrichment_results_table_html.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule plot_gowinda_enrichment_tajima_d_balancing:
    input:
        enrichment=rules.enrichment_tajima_d_balancing_gowinda.output.enrichment,
    output:
        count_plot=report(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.gowinda.enrichment.png",
            category="Balancing Selection",
            subcategory="{method}",
            labels=lambda wildcards: tajima_d_labels(
                wildcards, type="Enrichment Plot"
            ),
        ),
        qscore_plot=report(
            "results/balancing_selection/scikit-allel/{species}/{method}/{ppl}/{window}_{step}/{ppl}.{method}.top_{cutoff}.gowinda.qscore.png",
            category="Balancing Selection",
            subcategory="{method}",
            labels=lambda wildcards: tajima_d_labels(wildcards, type="Q-Score Plot"),
        ),
    params:
        title=lambda w: f"{w.ppl} {w.method.upper().replace('_', ' ')} ENRICHMENT (Window={w.window}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Step={int(float(w.step) * int(w.window))}{' SNPs' if w.method == 'moving_tajima_d' else ' bp'}, Top {float(w.cutoff)*100:.2f}%)",
    resources:
        mem_gb=8,
    log:
        "logs/balancing_selection/plot_gowinda_enrichment_tajima_d_balancing.{species}.{ppl}.{method}.{window}_{step}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/plot_gowinda_enrichment.py"
