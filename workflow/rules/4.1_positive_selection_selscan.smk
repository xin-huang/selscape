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


phasing_flag = ""

if selscan_config["unphased"] or main_config["ploidy"] > 2:
    phasing_flag = "--unphased"


rule extract_snp_pos:
    input:
        vcf=rules.polarize_1pop.output.vcf,
    output:
        map=temp(
            "results/polarized_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps.map"
        ),
    log:
        "logs/positive_selection/extract_snp_pos.{species}.{ppl}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        bcftools query -f "%CHROM\\t%CHROM:%POS:%REF:%ALT\\t%POS\\t%POS\\n" {input.vcf} > {output.map} 2> {log}
        """


rule estimate_selscan_scores:
    input:
        selscan=rules.download_selscan.output.selscan,
        vcf=rules.polarize_1pop.output.vcf,
        map=rules.extract_snp_pos.output.map,
    output:
        out=temp(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.chr{i}.{method}.out"
        ),
        formatted_out=temp(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.chr{i}.{method}.formatted.out"
        ),
        log=temp(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.chr{i}.{method}.log"
        ),
    params:
        output_prefix="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.chr{i}",
        phasing_flag=f"{phasing_flag}",
    resources:
        cpus=8,
    log:
        "logs/positive_selection/estimate_selscan_scores.{species}.{ppl}.{method}.{maf}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        {input.selscan} --vcf {input.vcf} --map {input.map} --{wildcards.method} --out {params.output_prefix} --threads {resources.cpus} --maf {wildcards.maf} {params.phasing_flag} 2> {log}
        awk -v chr={wildcards.i} 'BEGIN{{OFS="\\t"}}{{print chr,$2,$3,$4,$5,$6}}' {output.out} > {output.formatted_out} 2>> {log}
        """


rule normalize_selscan_scores:
    input:
        norm=rules.download_selscan.output.norm,
        scores=expand(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.chr{i}.{method}.formatted.out",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        scores=temp(
            expand(
                "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.chr{i}.{method}.formatted.out.100bins.norm",
                i=main_config["chromosomes"],
                allow_missing=True,
            )
        ),
        log=temp(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.{method}.log"
        ),
    log:
        "logs/positive_selection/normalize_selscan_scores.{species}.{ppl}.{method}.{maf}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        {input.norm} --files {input.scores} --log {output.log} --{wildcards.method} 2> {log}
        """


rule merge_selscan_scores:
    input:
        scores=rules.normalize_selscan_scores.output.scores,
    output:
        merged_scores="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.scores",
    log:
        "logs/positive_selection/merge_selscan_scores.{species}.{ppl}.{method}.{maf}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        cat {input.scores} | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$7}}' | sed '1iSNP\\tCHR\\tBP\\tnormalized_{wildcards.method}' > {output.merged_scores} 2> {log}
        """


rule plot_selscan:
    input:
        scores=rules.merge_selscan_scores.output.merged_scores,
    output:
        candidates="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.candidates.scores",
        plot=report(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.scores.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=selscan_labels,
        ),
    params:
        title=lambda w: f"{w.ppl} (MAF={w.maf}, Top {float(w.cutoff)*100:.2f}%)",
        score_column="normalized_{method}",
        use_absolute="TRUE",
        cutoff="{cutoff}",
        width=selscan_config["manhattan_plot_width"],
        height=selscan_config["manhattan_plot_height"],
        color1=selscan_config["manhattan_plot_color1"],
        color2=selscan_config["manhattan_plot_color2"],
    resources:
        mem_gb=16,
    log:
        "logs/positive_selection/plot_selscan.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/manhattan.R"


rule annotate_selscan_candidates:
    input:
        outliers="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.candidates.scores",
        annotation=expand(
            "results/annotated_data/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
            i=main_config["chromosomes"],
            ref_genome=main_config["ref_genome"],
            allow_missing=True,
        ),
    output:
        annotated_candidates="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.annotated.candidates",
    resources:
        mem_gb=32,
    log:
        "logs/positive_selection/annotate_selscan_candidates.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/get_annotated_candidates.py"


rule get_selscan_candidate_genes:
    input:
        selscan_candidates=rules.annotate_selscan_candidates.output.annotated_candidates,
    output:
        selscan_genes="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.candidate.genes",
    log:
        "logs/positive_selection/get_selscan_candidate_genes.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( sed '1d' {input.selscan_candidates} | grep -v ";" | awk '{{print $7}}' | sort | uniq > {output.selscan_genes} ) 2> {log} || true
        sed -i '1iGene' {output.selscan_genes} 2>> {log}
        """


rule selscan_candidate_genes_table_html:
    input:
        tsv="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.candidate.genes",
    output:
        html=report(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.candidate.genes.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: selscan_labels(wildcards, type="Gene List"),
        ),
    params:
        title=lambda w: f"{w.ppl} {w.method.upper()} (MAF={w.maf}, Top {float(w.cutoff)*100:.2f}%) CANDIDATE GENES",
    log:
        "logs/positive_selection/selscan_candidate_genes_table_html.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule enrichment_selscan_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        candidates=rules.annotate_selscan_candidates.output.annotated_candidates,
        total=expand(
            "results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        candidate_snps="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.candidate.snps.tsv",
        total_snps="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.total.snps.tsv",
        enrichment="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/positive_selection/enrichment_selscan_gowinda.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        sed '1d' {input.candidates} | awk '{{print "chr"$1"\\t"$2}}' > {output.candidate_snps} 2> {log}

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

        sed -i '1iGO_ID\\tavg_genes_sim\\tgenes_found\\tp_value\\tp_adjusted\\tgenes_uniq\\tgenes_max\\tgenes_total\\tdescription\\tgene_list' {output.enrichment}
        """


rule selscan_enrichment_results_table_html:
    input:
        tsv="results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.tsv",
    output:
        html=report(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.html",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: selscan_labels(
                wildcards, type="Enrichment Table"),
            ),
    params:
        title=lambda w: f"{w.ppl} {w.method.upper()} (MAF={w.maf}, Top {float(w.cutoff)*100:.2f}%) ENRICHMENT",
    log:
        "logs/positive_selection/selscan_enrichment_results_table_html.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule plot_gowinda_enrichment_selscan:
    input:
        enrichment=rules.enrichment_selscan_gowinda.output.enrichment,
    output:
        plot=report(
            "results/positive_selection/selscan/{species}/1pop/{ppl}/{method}_{maf}/{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.gowinda.enrichment.png",
            category="Positive Selection",
            subcategory="{method}",
            labels=lambda wildcards: selscan_labels(wildcards, type="Enrichment Plot"),
            ),
    params:
        title=lambda w: f"{w.ppl} {w.method.upper()} ENRICHMENT (MAF={w.maf}, Top {float(w.cutoff)*100:.2f}%)",
    resources:
        mem_gb=8,
    log:
        "logs/positive_selection/plot_gowinda_enrichment_selscan.{species}.{ppl}.normalized.{method}.maf_{maf}.top_{cutoff}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/plot_gowinda_enrichment.py"
