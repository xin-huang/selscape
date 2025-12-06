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

# Calculate allele counts
ac = gt.count_alleles()

# Calculate Tajima's D
d = allel.moving_tajima_d(ac, size=window_size, step=step_size)

# Window positions
window_pos = allel.moving_statistic(pos, np.mean, size=window_size, step=step_size)

# Remove NaN values
valid_mask = ~np.isnan(d)
window_pos = window_pos[valid_mask]
d = d[valid_mask]

# Save results
results_df = pd.DataFrame({
    'position': window_pos.astype(int),
    'tajima_d': d
})
results_df.to_csv(output_scores, sep='\t', index=False)
