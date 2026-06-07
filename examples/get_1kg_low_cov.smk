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



include: "get_shared_data.smk"

example_pops_low = ["YRI", "CHS"]
example_chrs_low = ["21"]


rule create_examples_low:
    input:
        expand("examples/data/Human/1kg_low_cov/{i}.vcf.gz", i=example_chrs_low),
        expand(
            "examples/data/Human/1kg_low_cov/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor.{i}.bed.gz",
            i=example_chrs_low,
        ),
        expand(
            "examples/data/Human/1kg_low_cov/repeats/hg19.{type}.autosomes.bed",
            type=["rmsk", "seg.dups", "simple.repeats"],
        ),
        "examples/data/Human/1kg_low_cov/metadata.txt",
        "examples/data/Human/1kg_low_cov/annotation/Human.hg19.gtf.gz",
        rules.download_gene2go.output.gene2go,
        "examples/data/Human/genome/hg19.chrom.sizes.bed",
        "examples/data/Human/genome/hg19.cytoBand.txt.gz",
        "examples/data/Human/genome/hg19.chrom.sizes.bed",
        "examples/data/Human/genome/hg19.cytoBand.txt.gz",

rule download_1kg_low_vcf:
    output:
        vcf="examples/data/Human/1kg_low_cov/{i}.vcf.gz",
        idx="examples/data/Human/1kg_low_cov/{i}.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
            -O {output.vcf}
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi \
            -O {output.idx}
        """


rule create_1kg_low_metadata:
    input:
        samples=rules.download_1KG_info.output.samples,
    output:
        metadata="examples/data/Human/1kg_low_cov/metadata.txt",
    params:
        pop1=example_pops_low[0],
        pop2=example_pops_low[1],
    shell:
        """
        grep -w {params.pop1} {input.samples} | awk -v pop={params.pop1} '{{print $1"\\t"pop}}' > {output.metadata}
        grep -w {params.pop2} {input.samples} | awk -v pop={params.pop2} '{{print $1"\\t"pop}}' >> {output.metadata}
        sed -i '1iSample\\tPopulation' {output.metadata}
        """


rule download_ncbi_annotation_hg19:
    output:
        gtf="examples/data/Human/1kg_low_cov/annotation/Human.hg19.gtf.gz",
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz \
            -O {output.gtf}
        """


rule download_ensembl_ancestral_alleles_hg19:
    output:
        anc_alleles=temp("examples/data/Human/1kg_low_cov/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2"),
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2 \
            -O {output.anc_alleles}
        tar -xjf {output.anc_alleles} -C examples/data/Human/1kg_low_cov/ancestral_alleles
        """


rule extract_anc_info_hg19:
    input:
        anc_alleles=rules.download_ensembl_ancestral_alleles_hg19.output.anc_alleles,
    output:
        bed=temp("examples/data/Human/1kg_low_cov/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor.{i}.bed"),
    params:
        fasta="examples/data/Human/1kg_low_cov/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{i}.fa",
    run:
        import pysam
        import re

        fasta = pysam.FastaFile(params.fasta)
        with open(output.bed, "wt") as out:
            for raw_chrom in fasta.references:
                match = re.search(r"GRCh\d+:(\d+|X|Y)", raw_chrom)
                if not match:
                    print(f"Skipping unrecognized chromosome name: {raw_chrom}")
                    continue
                chrom = match.group(1)
                seq = fasta.fetch(raw_chrom).upper()
                for pos, base in enumerate(seq):
                    if base in "ACGT":
                        out.write(f"{chrom}\t{pos}\t{pos+1}\t{base}\n")
        fasta.close()


rule compress_anc_info_hg19:
    input:
        bed=rules.extract_anc_info_hg19.output.bed,
    output:
        bed="examples/data/Human/1kg_low_cov/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor.{i}.bed.gz",
    shell:
        """
        bgzip -c {input.bed} > {output.bed}
        tabix -p bed {output.bed}
        """


rule download_hg19_repeat_files:
    output:
        rmsk="examples/data/Human/1kg_low_cov/repeats/hg19.rmsk.txt.gz",
        segdup="examples/data/Human/1kg_low_cov/repeats/hg19.genomicSuperDups.txt.gz",
        simrep="examples/data/Human/1kg_low_cov/repeats/hg19.simpleRepeat.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz -O {output.rmsk}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz -O {output.segdup}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz -O {output.simrep}
        """


rule convert_hg19_repeat_files:
    input:
        rmsk=rules.download_hg19_repeat_files.output.rmsk,
        segdup=rules.download_hg19_repeat_files.output.segdup,
        simrep=rules.download_hg19_repeat_files.output.simrep,
    output:
        rmsk="examples/data/Human/1kg_low_cov/repeats/hg19.rmsk.autosomes.bed",
        segdup="examples/data/Human/1kg_low_cov/repeats/hg19.seg.dups.autosomes.bed",
        simrep="examples/data/Human/1kg_low_cov/repeats/hg19.simple.repeats.autosomes.bed",
    shell:
        """
        zcat {input.rmsk}   | awk 'BEGIN{{OFS="\\t"}}$6!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $6,$7,$8,$11,$2,$10}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.rmsk}
        zcat {input.segdup} | awk 'BEGIN{{OFS="\\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$6,$7}}'   | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.segdup}
        zcat {input.simrep} | awk 'BEGIN{{OFS="\\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$11}}'     | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.simrep}
        """


rule download_hg19_genome_files:
    output:
        chrom_sizes=temp("examples/data/Human/genome/hg19.chrom.sizes"),
        cytoband="examples/data/Human/genome/hg19.cytoBand.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O {output.chrom_sizes}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz -O {output.cytoband}
        """

rule convert_hg19_chrom_sizes_to_bed:
    input:
        chrom_sizes="examples/data/Human/genome/hg19.chrom.sizes",
    output:
        bed="examples/data/Human/genome/hg19.chrom.sizes.bed",
    shell:
        """
        awk 'BEGIN{{OFS="\\t"}}{{print $1, 0, $2}}' {input.chrom_sizes} > {output.bed}
        """
