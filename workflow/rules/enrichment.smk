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


rule run_gowinda:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=expand(rules.convert_ncbi_go.output.go2gene, species=main_config["species"]),
        gtf=expand(rules.convert_ncbi_gtf.output.gtf, species=main_config["species"]),
        outliers="{prefix}.annotated.outliers",
    output:
        outlier_snps="{prefix}.outlier.snps.tsv",
        total_snps="{prefix}.total.snps.tsv",
        enrichment="{prefix}.gowinda.enrichment.tsv",
    resources:
        mem_gb=32,
        cpus=8,
    log:
        "logs/enrichment/{prefix}.gowinda.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        sed '1d' {input.outliers} | awk '{{print "chr"$1"\\t"$2}}' > {output.outlier_snps} 2> {log}
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


rule gowinda_html:
    input:
        tsv="{prefix}.gowinda.enrichment.tsv",
    output:
        html="{prefix}.gowinda.enrichment.html",
    params:
        title="",
    log:
        "logs/enrichment/{prefix}.html.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/tsv2html.R"


rule gowinda_plot:
    input:
        enrichment="{prefix}.gowinda.enrichment.tsv",
    output:
        plot="{prefix}.gowinda.enrichment.png",
    resources:
        mem_gb=8,
    log:
        "logs/enrichment/{prefix}.plot.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/plot_gowinda_enrichment.py"


ruleorder: enrichment_selscan_gowinda > enrichment_selscan_xp_gowinda > enrichment_tajima_d_gowinda > enrichment_delta_tajima_d_gowinda > enrichment_betascan_gowinda > enrichment_tajima_d_balancing_gowinda > run_gowinda

ruleorder: selscan_enrichment_results_table_html > selscan_xp_enrichment_results_table_html > tajima_d_enrichment_results_table_html > delta_tajima_d_enrichment_results_table_html > betascan_enrichment_results_table_html > tajima_d_balancing_enrichment_results_table_html > gowinda_html

ruleorder: plot_gowinda_enrichment_selscan > plot_gowinda_enrichment_selscan_xp > plot_gowinda_enrichment_tajima_d > plot_gowinda_enrichment_delta_tajima_d > plot_gowinda_enrichment_betascan > plot_gowinda_enrichment_tajima_d_balancing > gowinda_plot	
