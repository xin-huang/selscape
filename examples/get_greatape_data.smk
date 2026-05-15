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



example_chr_ape = ["21"]


rule create_examples_greatape:
    input:
        expand("examples/data/greatape/Pan/{i}.filteranno.vcf.gz", i=example_chr_ape),
        expand("examples/data/greatape/Pan/{i}.filteranno.vcf.gz.tbi", i=example_chr_ape),
        "examples/data/greatape/metadata.txt",


rule download_greatape_vcf:
    output:
        vcf="examples/data/greatape/Pan/{i}.filteranno.vcf.gz",
        idx="examples/data/greatape/Pan/{i}.filteranno.vcf.gz.tbi",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/merged_segregating/Pan/Pan_all_filtered/chr{wildcards.i}.filteranno.vcf.gz \
            -O {output.vcf}
        bcftools index -t {output.vcf}
        """


rule download_greatape_metadata:
    output:
        raw="examples/data/greatape/metadata_full.txt",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/metadata_full.txt \
            -O {output.raw}
        """


rule create_greatape_metadata:
    input:
        metadata=rules.download_greatape_metadata.output.raw,
    output:
        metadata="examples/data/greatape/metadata.txt",
    shell:
        """
        grep -v captive {input.metadata} | awk 'NR>1 {{print $4"\\t"$2}}' > {output.metadata}
        sed -i '1iSample\\tPopulation' {output.metadata}
        """
