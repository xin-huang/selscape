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


rule download_ref_genome:
    output:
        ref_genome="resources/data/{ref_genome}/{ref_genome}.fa",
    params:
        dir="resources/data/ref_genome",
    log:
        "logs/download/download_ref_genome.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenpath/{wildcards.ref_genome}/bigZips/{wildcards.ref_genome}.fa.gz -O {params.dir}/{wildcards.ref_genome}.fa.gz > {log} 2>&1
        gzip -d {params.dir}/{wildcards.ref_genome}.fa.gz > {log} 2>&1
        """


rule download_repeats:
    output:
        rmsk="resources/data/{ref_genome}/{ref_genome}.rmsk.txt.gz",
        segdup="resources/data/{ref_genome}/{ref_genome}.genomicSuperDups.txt.gz",
        simrep="resources/data/{ref_genome}/{ref_genome}.simpleRepeat.txt.gz",
    log:
        "logs/download/download_repeats.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        for repeat_type in rmsk genomicSuperDups simpleRepeat; do
            url="https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.ref_genome}/database/${{repeat_type}}.txt.gz"
            output_file="resources/data/{wildcards.ref_genome}/{wildcards.ref_genome}.${{repeat_type}}.txt.gz"
            
            if wget -q --spider "$url" 2>>{log}; then
                wget -c "$url" -O "$output_file" >>{log} 2>&1
            else
                touch "$output_file"
            fi
        done
        """


rule download_selscan:
    output:
        selscan="resources/tools/selscan/selscan-2.0.3",
        norm="resources/tools/selscan/norm",
    log:
        "logs/download/download_selscan.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        rm -rf selscan-2.0.3* >> {log} 2>&1
        wget -c https://github.com/szpiech/selscan/archive/refs/tags/v2.0.3.tar.gz >> {log} 2>&1
        tar -xvzf v2.0.3.tar.gz >> {log} 2>&1
        mv selscan-2.0.3/bin/linux/selscan-2.0.3 {output.selscan} >> {log} 2>&1
        mv selscan-2.0.3/bin/linux/norm {output.norm} >> {log} 2>&1
        rm -rf v2.0.3.tar.gz selscan-2.0.3/ >> {log} 2>&1
        """


rule download_betascan:
    output:
        betascan="resources/tools/betascan/BetaScan.py",
    log:
        "logs/download/download_betascan.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        git clone https://github.com/ksiewert/BetaScan > {log} 2>&1
        mv BetaScan/BetaScan.py {output.betascan} > {log} 2>&1
        rm -rf BetaScan > {log} 2>&1
        """


rule download_annovar_db:
    output:
        ref_gene="resources/tools/annovar/{ref_genome}_db/{ref_genome}_refGene.txt",
        ref_seq="resources/tools/annovar/{ref_genome}_db/{ref_genome}_refGeneMrna.fa",
    params:
        db_dir="resources/tools/annovar/{ref_genome}_db",
        seq_dir="resources/tools/annovar/{ref_genome}_db/seq",
    log:
        "logs/download/download_annovar_db.{ref_genome}.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        resources/tools/annovar/annotate_variation.pl -buildver {wildcards.ref_genome} -downdb refGene {params.db_dir} > {log} 2>&1
        resources/tools/annovar/annotate_variation.pl -buildver {wildcards.ref_genome} -downdb seq {params.seq_dir} > {log} 2>&1
        resources/tools/annovar/retrieve_seq_from_fasta.pl {output.ref_gene} -seqfile {params.seq_dir}/{wildcards.ref_genome}.fa -format refGene -outfile {output.ref_seq} > {log} 2>&1
        """


rule download_gowinda:
    output:
        gowinda="resources/tools/gowinda/Gowinda-1.12.jar",
    log:
        "logs/download/download_gowinda.log",
    conda:
        "../envs/selscape-env.yaml",
    shell:
        """
        wget -c https://sourceforge.net/projects/gowinda/files/Gowinda-1.12.jar/download -O {output.gowinda} > {log} 2>&1
        """
