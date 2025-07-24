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


rule extract_biallelic_snps:
    input:
        vcf=get_vcf_input_path,
    output:
        vcf="results/processed_data/{species}/all/chr{i}.biallelic.snps.vcf.gz",
        idx="results/processed_data/{species}/all/chr{i}.biallelic.snps.vcf.gz.tbi",
    log:
        "logs/preprocess/extract_biallelic_snps.{species}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule create_pair_info:
    input:
        metadata=main_config["metadata"],
    output:
        pair_info="results/samples/{species}/{pair}/{pair}.list",
    log:
        "logs/preprocess/create_pair_info.{species}.{pair}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)
        awk -v pop=$pop1 '$2 == pop {{print $1"\\t"$2}}' {input.metadata} > {output.pair_info} 2> {log}
        awk -v pop=$pop2 '$2 == pop {{print $1"\\t"$2}}' {input.metadata} >> {output.pair_info} 2>> {log}
        """


rule annotate_biallelic_snps:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        ref_gene=rules.download_annovar_db.output.ref_gene,
    output:
        avinput="results/annotated_data/{species}/all/chr{i}.biallelic.snps.{ref_genome}.avinput",
        txt="results/annotated_data/{species}/all/chr{i}.biallelic.snps.{ref_genome}_multianno.txt",
    resources:
        cpus=8,
        mem_gb=32,
    params:
        output_prefix="results/annotated_data/{species}/all/chr{i}.biallelic.snps",
        db_dir="resources/tools/annovar/{ref_genome}_db",
    log:
        "logs/preprocess/annotate_biallelic_snps.{species}.chr{i}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\n" {input.vcf} > {output.avinput} 2> {log}
        resources/tools/annovar/table_annovar.pl {output.avinput} {params.db_dir} \
            -buildver {wildcards.ref_genome} \
            -out {params.output_prefix} \
            -protocol refGene \
            -operation g \
            -nastring . \
            -polish \
            -remove \
            -thread {resources.cpus} >> {log} 2>&1
        """


rule extract_pop_data:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        metadata=main_config["metadata"],
    output:
        vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz",
        idx="results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz.tbi",
    log:
        "logs/preprocess/extract_pop_data.{species}.{ppl}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools view {input.vcf} -S <(sed '1d' {input.metadata} | grep -w {wildcards.ppl} | awk '{{print $1}}') --force-samples | \
        bcftools view -i "AC>0 && AC<AN" | \
        bcftools annotate -x ^FORMAT/GT | \
        bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule extract_pair_data:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        samples=rules.create_pair_info.output.pair_info,
    output:
        vcf="results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz",
        idx="results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz.tbi",
    log:
        "logs/preprocess/extract_pair_data.{species}.{pair}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools view {input.vcf} -S <(awk '{{print $1}}' {input.samples}) --force-samples | \
        bcftools view -i "AC>0 && AC<AN" | \
        bcftools annotate -x ^FORMAT/GT | \
        bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule extract_1pop_exonic_data:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        anno=rules.annotate_biallelic_snps.output.txt,
    output:
        vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    params:
        condition=lambda wildcards: (
            "$9~/^synonymous/"
            if wildcards.mut_type == "syn"
            else "$9~/^nonsynonymous/"
        ),
    resources:
        mem_gb=32,
    log:
        "logs/preprocess/extract_1pop_exonic_data.{species}.{ppl}.chr{i}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools view {input.vcf} -R <(awk '{params.condition}{{print $1"\\t"$2}}' {input.anno}) |\
            bgzip -c > {output.vcf} ) 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_1pop_exonic_data:
    input:
        vcfs=expand(
            "results/processed_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    log:
        "logs/preprocess/concat_1pop_exonic_data.{species}.{ppl}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule extract_2pop_exonic_data:
    input:
        vcf=rules.extract_pair_data.output.vcf,
        anno=rules.annotate_biallelic_snps.output.txt,
    output:
        vcf="results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    params:
        condition=lambda wildcards: (
            "$9~/^synonymous/"
            if wildcards.mut_type == "syn"
            else "$9~/^nonsynonymous/"
        ),
    resources:
        mem_gb=32,
    log:
        "logs/preprocess/extract_2pop_exonic_data.{species}.{pair}.chr{i}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools view {input.vcf} -R <(awk '{params.condition}{{print $1"\\t"$2}}' {input.anno}) |\
            bgzip -c > {output.vcf} ) 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_2pop_exonic_data:
    input:
        vcfs=expand(
            "results/processed_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        vcf="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    log:
        "logs/preprocess/concat_2pop_exonic_data.{species}.{pair}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule convert_ncbi_gtf:
    input:
        gtf=main_config["genome_annotation"],
    output:
        gtf="results/annotated_data/{species}.gowinda.gtf",
    log:
        "logs/gene_enrichment/convert_ncbi_gtf.{species}.log",
    conda:
        "../envs/selscape-env.yaml",
    script:
        "../scripts/ncbi_gtf2gowinda.py"


rule convert_ncbi_go:
    input:
        gtf=main_config["genome_annotation"],
        gene2go=main_config["gene2go"]
    output:
        go2gene="results/annotated_data/{species}.gowinda.go2gene",
    params:
        tax_id=main_config["tax_id"],
    log:
        "logs/gene_enrichment/convert_ncbi_go.{species}.log",
    conda:
        "../envs/selscape-env.yaml",
    script:
        "../scripts/ncbi_go2gowinda.py"
