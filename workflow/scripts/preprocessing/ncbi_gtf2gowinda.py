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


import re
import gzip


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh


def convert_chromosome_name(name):
    """
    Flexible chromosome name conversion for any species.
    Handles common NCBI RefSeq name formats.
    """
    base_name = name.split('.')[0]
    
    if base_name.startswith('NC_'):
        try:
            chr_num = base_name.split('_')[1].lstrip('0')
            
            if chr_num == '':  
                chr_num = '0'
            
            return f'chr{chr_num}'
            
        except (IndexError, ValueError):
            return name
    
    elif name.startswith('chr'):
        return name  
    
    elif re.match(r'^[0-9]+$|^[IVX]+$', name):
        return f'chr{name}'
    
    else:
        return name


def keep_chromosome(chr_name):
    """
    Flexible filtering for chromosomes to include.
    Excludes common unwanted sequences but keeps species-specific ones.
    """
    chr_lower = chr_name.lower()
    
    exclude_patterns = [
        'un',
        'random',
        'alt',
        'fix',
        'scaffold',
        'contig',
    ]
    
    for pattern in exclude_patterns:
        if pattern in chr_lower:
            return False
    
    return True


input_file = snakemake.input.gtf
output_file = snakemake.output.gtf

with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        chrom_name = fields[0]
        feature_type = fields[2]
        attr_field = fields[8]

        chr_name = convert_chromosome_name(chrom_name)
        
        if not keep_chromosome(chr_name):
            continue
            
        if feature_type != "gene":
            continue

        gene_id_match = re.search(r'gene_id "[^"]+";', attr_field)
        if not gene_id_match:
            continue
            
        gene_id_str = gene_id_match.group(0)
        
        fields[0] = chr_name
        fields[2] = 'exon'  
        out_line = '\t'.join(fields[:8]) + '\t' + gene_id_str + '\n'
        outfile.write(out_line)
