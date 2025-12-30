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


rule polarize_1pop:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        anc_alleles=get_anc_allele_bed,
    output:
        vcf=temp("results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz"),
        idx=temp("results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.vcf.gz.tbi"),
    resources:
        mem_gb=32,
    log:
        "logs/polarization/polarize_1pop.{species}.{ppl}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/polarize_snps.py"


rule polarize_2pop:
    input:
        vcf=rules.extract_pair_data.output.vcf,
        anc_alleles=get_anc_allele_bed,
    output:
        vcf=temp("results/polarized_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz"),
        idx=temp("results/polarized_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.snps.vcf.gz.tbi"),
    resources:
        mem_gb=32,
    log:
        "logs/polarization/polarize_2pop.{species}.{pair}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/polarize_snps.py"


rule polarize_1pop_exonic_data:
    input:
        vcf=rules.extract_1pop_exonic_data.output.vcf,
        anc_alleles=get_anc_allele_bed,
    output:
        vcf=temp("results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz"),
        idx=temp("results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi"),
    resources:
        mem_gb=32,
    log:
        "logs/polarization/polarize_1pop_exonic_data.{species}.{ppl}.chr{i}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/polarize_snps.py"


rule polarize_2pop_exonic_data:
    input:
        vcf=rules.extract_2pop_exonic_data.output.vcf,
        anc_alleles=get_anc_allele_bed,
    output:
        vcf=temp("results/polarized_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz"),
        idx=temp("results/polarized_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi"),
    resources:
        mem_gb=32,
    log:
        "logs/polarization/polarize_2pop_exonic_data.{species}.{pair}.chr{i}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/polarize_snps.py"


rule concat_polarized_1pop_exonic_data:
    input:
        vcfs=expand(
            "results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        vcf="results/polarized_data/{species}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/polarized_data/{species}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    log:
        "logs/polarization/concat_polarized_1pop_exonic_data.{species}.{ppl}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_polarized_2pop_exonic_data:
    input:
        vcfs=expand(
            "results/polarized_data/{species}/2pop/{pair}/{pair}.chr{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
            i=main_config["chromosomes"],
            allow_missing=True,
        ),
    output:
        vcf="results/polarized_data/{species}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/polarized_data/{species}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    log:
        "logs/polarization/concat_polarized_2pop_exonic_data.{species}.{pair}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """

rule test_hwe_polarized:
    input:
        vcf=rules.polarize_1pop.output.vcf,
    output:
        hwe_outliers=temp(
            "results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.hwe.outliers"
        ),
    params:
        output_prefix="results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps",
        hwe_threshold=main_config["hwe_pvalue"],
    log:
        "logs/polarization/test_hwe_polarized.{species}.{ppl}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        plink --vcf {input.vcf} --hardy --out {params.output_prefix} --set-missing-var-ids @:# 2> {log}
        ( awk '$7>$8' {params.output_prefix}.hwe | \
        sed '1d' | \
        awk -v threshold={params.hwe_threshold} '$9<threshold{{print $2}}' | \
        awk -F ":" '{{print $1"\\t"$2}}' > {output.hwe_outliers} ) 2>> {log}
        """


rule remove_repeats_polarized:
    input:
        vcf=rules.polarize_1pop.output.vcf,
        hwe_outliers=rules.test_hwe_polarized.output.hwe_outliers,
    output:
        vcf=temp("results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.repeats.removed.vcf.gz"),
        idx=temp("results/polarized_data/{species}/1pop/{ppl}/{ppl}.chr{i}.biallelic.snps.repeats.removed.vcf.gz.tbi"),
    params:
        rmsk=main_config["rmsk"] if "rmsk" in main_config else "",
        seg_dup=main_config["seg_dup"] if "seg_dup" in main_config else "",
        sim_rep=main_config["sim_rep"] if "sim_rep" in main_config else "",
    log:
        "logs/polarization/remove_repeats_polarized.{species}.{ppl}.chr{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        set -e
        exec 2> {log}

        cmd="bcftools view {input.vcf}"
        if [ -n "{params.rmsk}" ] && [ -s "{params.rmsk}" ]; then
            cmd="$cmd | bcftools view -T ^{params.rmsk}"
        fi
        if [ -n "{params.sim_rep}" ] && [ -s "{params.sim_rep}" ]; then
            cmd="$cmd | bcftools view -T ^{params.sim_rep}"
        fi
        if [ -n "{params.seg_dup}" ] && [ -s "{params.seg_dup}" ]; then
            cmd="$cmd | bcftools view -T ^{params.seg_dup}"
        fi
        if [ -s "{input.hwe_outliers}" ]; then
            cmd="$cmd | bcftools view -T ^{input.hwe_outliers}"
        fi
        cmd="$cmd | bgzip -c > {output.vcf}"
        eval "$cmd"

        tabix -p vcf {output.vcf}
        """
