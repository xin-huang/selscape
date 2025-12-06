import allel
import numpy as np
import pandas as pd


# Snakemake inputs
vcf_file = snakemake.input.vcf
output_scores = snakemake.output.scores
window_size = snakemake.params.window_size
step_size = snakemake.params.step_size

# Read VCF
callset = allel.read_vcf(vcf_file, fields=['variants/POS', 'calldata/GT'])
gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']

# Split populations
n_samples = gt.shape[1]
pop1_indices = list(range(n_samples // 2))
pop2_indices = list(range(n_samples // 2, n_samples))

# Calculate allele counts
ac1 = gt.count_alleles(subpop=pop1_indices)
ac2 = gt.count_alleles(subpop=pop2_indices)

# Calculate Delta Tajima's D
d1 = allel.moving_tajima_d(ac1, size=window_size, step=step_size)
d2 = allel.moving_tajima_d(ac2, size=window_size, step=step_size)
delta = d1 - d2

# Window positions
window_pos = allel.moving_statistic(pos, np.mean, size=window_size, step=step_size)

# Remove NaN values
valid_mask = ~np.isnan(delta)
window_pos = window_pos[valid_mask]
delta = delta[valid_mask]

delta_z = (delta - np.mean(delta)) / np.std(delta)

# Save results
results_df = pd.DataFrame({
    'position': window_pos.astype(int),
    'delta_tajima_d': delta,
    'delta_tajima_d_standardized': delta_z
})
results_df.to_csv(output_scores, sep='\t', index=False)
