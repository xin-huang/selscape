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


import allel
import numpy as np
import pandas as pd

# Snakemake inputs
vcf_file = snakemake.input.vcf
pair_info_file = snakemake.input.pair_info
output_scores = snakemake.output.scores
window_size = int(snakemake.params.window_size)
step_size_ratio = float(snakemake.params.step_size_ratio)
step_size = int(step_size_ratio * window_size)

# Read pair info to get population assignments
pair_info = pd.read_csv(pair_info_file, sep='\t')
pop1_name = pair_info.iloc[0, 1]
pop2_name = pair_info.iloc[-1, 1]

pop1_samples = pair_info[pair_info.iloc[:, 1] == pop1_name].iloc[:, 0].tolist()
pop2_samples = pair_info[pair_info.iloc[:, 1] == pop2_name].iloc[:, 0].tolist()

# Read VCF
callset = allel.read_vcf(vcf_file, fields=['variants/POS', 'calldata/GT', 'samples'])
gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']
samples = callset['samples']

# Get sample indices
pop1_indices = [i for i, s in enumerate(samples) if s in pop1_samples]
pop2_indices = [i for i, s in enumerate(samples) if s in pop2_samples]

# Calculate allele counts for each population
ac1 = gt.count_alleles(subpop=pop1_indices)
ac2 = gt.count_alleles(subpop=pop2_indices)

# Calculate Delta Tajima's D
delta = allel.moving_delta_tajima_d(ac1, ac2, size=window_size, step=step_size)

# Window positions
window_start = allel.moving_statistic(pos, np.min, size=window_size, step=step_size)
window_end = allel.moving_statistic(pos, np.max, size=window_size, step=step_size)

# Remove NaN values
valid_mask = ~np.isnan(delta)
window_start = window_start[valid_mask]
window_end = window_end[valid_mask]
delta = delta[valid_mask]

# Save results
results_df = pd.DataFrame({
    'window_start': window_start.astype(int),
    'window_end': window_end.astype(int),
    'n_snps': window_size,
    'delta_tajima_d': delta,
})
results_df.to_csv(output_scores, sep='\t', index=False)
