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


import gzip
import re
from collections import defaultdict


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

gtf_file = snakemake.input.gtf
gene2go_file = snakemake.input.gene2go
tax_id_target = str(snakemake.params.tax_id)

geneid2go = defaultdict(set)
go_term_description = {}

with gzip.open(gene2go_file, "rt") as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        tax_id, geneid, go_id, go_term = fields[0], fields[1], fields[2], " ".join(fields[4:-2])
        if tax_id == tax_id_target:
            geneid2go[geneid].add(go_id)
            go_term_description[go_id] = go_term

go2genes = defaultdict(set)

with gzip.open(gtf_file, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "gene":
            continue

        gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
        gene_id = gene_id_match.group(1) if gene_id_match else "NA"

        geneid_match = re.search(r'db_xref "GeneID:(\d+)"', fields[8])
        geneid = geneid_match.group(1) if geneid_match else "NA"

        if geneid in geneid2go:
            for go_id in geneid2go[geneid]:
                go2genes[go_id].add(gene_id)

with open(snakemake.output.go2gene, "w") as out:
    for go_id in sorted(go2genes.keys()):
        description = go_term_description.get(go_id, "NA")
        genes_str = " ".join(sorted(go2genes[go_id]))
        out.write(f"{go_id}\t{description}\t{genes_str}\n")
