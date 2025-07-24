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


import cyvcf2
import gzip
import subprocess
import sys
from cyvcf2 import VCF, Writer


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

input_vcf = snakemake.input.vcf
input_anc_alleles = snakemake.input.anc_alleles
output_vcf = snakemake.output.vcf

# Load ancestral alleles
anc_alleles = dict()
with gzip.open(input_anc_alleles, "rt") as f:
    for l in f:
        e = l.rstrip().split("\t")
        anc_alleles[f"{e[0]}:{e[2]}"] = e[3]

# Process VCF
vcf = VCF(input_vcf)

# Add AA field to header
vcf.add_info_to_header({'ID': 'AA', 'Description': 'Ancestral Allele', 'Type': 'String', 'Number': '1'})

w = Writer(output_vcf, vcf, mode="wz")
for v in vcf:
    key = f"{v.CHROM}:{v.POS}"
    if key in anc_alleles.keys():
        ancestral_base = anc_alleles[key]
        
        # Skip if neither REF nor ALT matches ancestral
        if (v.REF != ancestral_base) and (v.ALT[0] != ancestral_base):
            continue
            
        # Add ancestral allele annotation
        v.INFO["AA"] = ancestral_base
        
        # Flip if REF doesn't match ancestral
        if v.REF != ancestral_base:
            tmp = v.ALT[0]
            v.ALT = [v.REF]
            v.REF = tmp
            v.genotypes = [[int(not g[0]), int(not g[1]), True] for g in v.genotypes]
        
        w.write_record(v)

w.close()
vcf.close()
subprocess.run(["tabix", "-p", "vcf", output_vcf], check=True)
