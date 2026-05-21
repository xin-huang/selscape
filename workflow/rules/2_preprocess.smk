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
        vcf=temp("results/processed_data/{species}/{dataset}/all/{i}.biallelic.snps.vcf.gz"),
        idx=temp("results/processed_data/{species}/{dataset}/all/{i}.biallelic.snps.vcf.gz.tbi"),
    log:
        "logs/preprocess/extract_biallelic_snps.{species}.{dataset}.{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule create_pair_info:
    input:
        metadata=get_metadata,
    output:
        pair_info="results/samples/{species}/{dataset}/{pair}/{pair}.list",
    log:
        "logs/preprocess/create_pair_info.{species}.{dataset}.{pair}.log",
    conda:
        "../envs/selscape-env.yaml"
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
        idx=rules.extract_biallelic_snps.output.idx,
        ref_gene=rules.download_annovar_db.output.ref_gene,
    output:
        avinput=temp("results/annotated_data/{species}/{dataset}/all/{i}.biallelic.snps.{ref_genome}.avinput"),
        txt="results/annotated_data/{species}/{dataset}/all/{i}.biallelic.snps.{ref_genome}_multianno.txt",
    resources:
        cpus=8,
        mem_mb=32000,
    params:
        output_prefix="results/annotated_data/{species}/{dataset}/all/{i}.biallelic.snps",
        db_dir="resources/tools/annovar/{ref_genome}_db",
    log:
        "logs/preprocess/annotate_biallelic_snps.{species}.{dataset}.{i}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
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
        idx=rules.extract_biallelic_snps.output.idx,
        metadata=get_metadata,
    output:
        vcf=temp("results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.snps.vcf.gz"),
        idx=temp("results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.snps.vcf.gz.tbi"),
    log:
        "logs/preprocess/extract_pop_data.{species}.{dataset}.{ppl}.{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools view {input.vcf} -S <(sed '1d' {input.metadata} | grep -w {wildcards.ppl} | awk '{{print $1}}') --force-samples | \
        bcftools view -i "AC>0 && AC<AN" | \
        bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule extract_pair_data:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        idx=rules.extract_biallelic_snps.output.idx,
        samples=rules.create_pair_info.output.pair_info,
    output:
        vcf=temp("results/processed_data/{species}/{dataset}/2pop/{pair}/{pair}.{i}.biallelic.snps.vcf.gz"),
        idx=temp("results/processed_data/{species}/{dataset}/2pop/{pair}/{pair}.{i}.biallelic.snps.vcf.gz.tbi"),
    log:
        "logs/preprocess/extract_pair_data.{species}.{dataset}.{pair}.{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools view {input.vcf} -S <(awk '{{print $1}}' {input.samples}) --force-samples | \
        bcftools view -i "AC>0 && AC<AN" | \
        bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule extract_1pop_exonic_data:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        idx=rules.extract_pop_data.output.idx,
        anno=rules.annotate_biallelic_snps.output.txt,
    output:
        vcf=temp("results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz"),
        idx=temp("results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi"),
    params:
        condition=lambda wildcards: (
            "$9~/^synonymous/"
            if wildcards.mut_type == "syn"
            else "$9~/^nonsynonymous/"
        ),
    resources:
        mem_mb=32000,
    log:
        "logs/preprocess/extract_1pop_exonic_data.{species}.{dataset}.{ppl}.{i}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools view {input.vcf} -R <(awk '{params.condition}{{print $1"\\t"$2}}' {input.anno}) |\
            bgzip -c > {output.vcf} ) 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_1pop_exonic_data:
    input:
        vcfs=lambda wc: expand(
            rules.extract_1pop_exonic_data.output.vcf,
            i=get_chromosomes(wc),
            allow_missing=True,
        ),
        idxs=lambda wc: expand(
            rules.extract_1pop_exonic_data.output.idx,
            i=get_chromosomes(wc),
            allow_missing=True,
        ),
    output:
        vcf="results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    log:
        "logs/preprocess/concat_1pop_exonic_data.{species}.{dataset}.{ppl}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule extract_2pop_exonic_data:
    input:
        vcf=rules.extract_pair_data.output.vcf,
        idx=rules.extract_pair_data.output.idx,
        anno=rules.annotate_biallelic_snps.output.txt,
    output:
        vcf=temp("results/processed_data/{species}/{dataset}/2pop/{pair}/{pair}.{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz"),
        idx=temp("results/processed_data/{species}/{dataset}/2pop/{pair}/{pair}.{i}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi"),
    params:
        condition=lambda wildcards: (
            "$9~/^synonymous/"
            if wildcards.mut_type == "syn"
            else "$9~/^nonsynonymous/"
        ),
    resources:
        mem_mb=32000,
    log:
        "logs/preprocess/extract_2pop_exonic_data.{species}.{dataset}.{pair}.{i}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools view {input.vcf} -R <(awk '{params.condition}{{print $1"\\t"$2}}' {input.anno}) |\
            bgzip -c > {output.vcf} ) 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_2pop_exonic_data:
    input:
        vcfs=lambda wc: expand(
            rules.extract_2pop_exonic_data.output.vcf,
            i=get_chromosomes(wc),
            allow_missing=True,
        ),
        idxs=lambda wc: expand(
            rules.extract_2pop_exonic_data.output.idx,
            i=get_chromosomes(wc),
            allow_missing=True,
        ),
    output:
        vcf="results/processed_data/{species}/{dataset}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz",
        idx="results/processed_data/{species}/{dataset}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.{ref_genome}.vcf.gz.tbi",
    log:
        "logs/preprocess/concat_2pop_exonic_data.{species}.{dataset}.{pair}.{mut_type}.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule test_hwe:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        idx=rules.extract_pop_data.output.idx,
    output:
        hwe_outliers=temp(
            "results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.snps.hwe.outliers"
        ),
    params:
        output_prefix="results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.snps",
        hwe_threshold=get_hwe_pvalue,
    log:
        "logs/preprocess/test_hwe.{species}.{dataset}.{ppl}.{i}.log",
    conda:
        "../envs/selscape-env.yaml"
    shell:
        """
        ( bcftools annotate --set-id '%CHROM:%POS' {input.vcf} | \
        plink --vcf /dev/stdin --hardy --out {params.output_prefix} ) 2> {log}
        ( awk '$7>$8' {params.output_prefix}.hwe | \
        sed '1d' | \
        awk -v threshold={params.hwe_threshold} '$9<threshold{{print $2}}' | \
        awk -F ":" '{{print $1"\\t"$2}}' > {output.hwe_outliers} ) 2>> {log}
        """


rule remove_repeats:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        idx=rules.extract_pop_data.output.idx,
        hwe_outliers=rules.test_hwe.output.hwe_outliers,
    output:
        vcf=temp("results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.snps.repeats.removed.vcf.gz"),
        idx=temp("results/processed_data/{species}/{dataset}/1pop/{ppl}/{ppl}.{i}.biallelic.snps.repeats.removed.vcf.gz.tbi"),
    params:
        rmsk=get_rmsk,
        seg_dup=get_seg_dup,
        sim_rep=get_sim_rep,
    log:
        "logs/preprocess/remove_repeats.{species}.{dataset}.{ppl}.{i}.log",
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


rule convert_ncbi_gtf:
    input:
        gtf=get_genome_annotation,
    output:
        gtf="results/annotated_data/{species}/{dataset}.gowinda.gtf",
    log:
        "logs/gene_enrichment/convert_ncbi_gtf.{species}.{dataset}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/preprocessing/ncbi_gtf2gowinda.py"


rule convert_ncbi_go:
    input:
        gtf=get_genome_annotation,
        gene2go=get_gene2go,
    output:
        go2gene="results/annotated_data/{species}/{dataset}.gowinda.go2gene",
    params:
        tax_id=get_tax_id,
    log:
        "logs/gene_enrichment/convert_ncbi_go.{species}.{dataset}.log",
    conda:
        "../envs/selscape-env.yaml"
    script:
        "../scripts/preprocessing/ncbi_go2gowinda.py"
