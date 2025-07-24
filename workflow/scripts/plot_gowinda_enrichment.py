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


import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

# Read Gowinda enrichment results
try:
    enrichment_data = pd.read_csv(
        snakemake.input['enrichment'],
        sep='\t',
        header=0,
    )
except pd.errors.EmptyDataError:
    enrichment_data = pd.DataFrame()

if enrichment_data.empty:
    plt.figure()
    plt.title("No data available")
    plt.savefig(snakemake.output['count_plot'], dpi=300, bbox_inches='tight')
    plt.savefig(snakemake.output['qscore_plot'], dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)


# Take top 20 terms by p_adjusted 
plot_data = enrichment_data.sort_values('p_adjusted').head(20).copy()

# Calculate qscore
plot_data['qscore'] = -np.log10(plot_data['p_adjusted'])

# Shorten descriptions
plot_data['short_description'] = plot_data['description'].apply(
    lambda x: x[:47] + "..." if len(x) > 50 else x
)

# Create individual colors for each bar based on p_adjusted
def get_color(p_val):
    if p_val < 0.001:
        return '#0000FF'  # Blue
    elif p_val < 0.01:
        return '#800080'  # Purple
    elif p_val < 0.05:
        return '#FF0000'  # Red
    else:
        return '#808080'  # Gray

plot_data['color'] = plot_data['p_adjusted'].apply(get_color)

# Create custom legend
legend_elements = [
    Patch(facecolor='#0000FF', label='< 0.001'),
    Patch(facecolor='#800080', label='< 0.01'),
    Patch(facecolor='#FF0000', label='< 0.05'),
    Patch(facecolor='#808080', label='>= 0.05')
]

# ==== PLOT 1: COUNT PLOT ====
# Sort by genes_found for plotting 
count_data = plot_data.sort_values('genes_found', ascending=True).copy()

# Create the plot
fig, ax = plt.subplots(figsize=(12, max(8, len(count_data) * 0.4)))

# Create horizontal bar plot using genes_found (count)
bars = ax.barh(range(len(count_data)), count_data['genes_found'], color=count_data['color'])

# Customize the plot
ax.set_yticks(range(len(count_data)))
ax.set_yticklabels(count_data['short_description'], fontsize=9)
ax.set_xlabel('Count', fontsize=11)

# Center the title better
plt.suptitle(f'GO Enrichment Analysis (Gowinda) - Top {len(count_data)} terms',
            fontsize=14, fontweight='bold', x=0.45)

# Remove spines and grid
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)

ax.legend(handles=legend_elements, title='p.adjust', loc='lower right')

plt.tight_layout()
plt.savefig(snakemake.output['count_plot'], dpi=300, bbox_inches='tight')
plt.close()

# ==== PLOT 2: QSCORE PLOT ====
# Sort by qscore for plotting
qscore_data = plot_data.sort_values('qscore', ascending=True).copy()

# Create the plot
fig, ax = plt.subplots(figsize=(12, max(8, len(qscore_data) * 0.4)))

# Create horizontal bar plot using qscore
bars = ax.barh(range(len(qscore_data)), qscore_data['qscore'], color=qscore_data['color'])

# Customize the plot
ax.set_yticks(range(len(qscore_data)))
ax.set_yticklabels(qscore_data['short_description'], fontsize=9)
ax.set_xlabel('qscore (-log10(p.adjusted))', fontsize=11)

# Center the title
plt.suptitle(f'GO Enrichment Analysis (Gowinda) - Top {len(qscore_data)} terms',
            fontsize=14, fontweight='bold', x=0.45)

# Remove spines and grid
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)

ax.legend(handles=legend_elements, title='p.adjust', loc='lower right')

plt.tight_layout()
plt.savefig(snakemake.output['qscore_plot'], dpi=300, bbox_inches='tight')
plt.close()
