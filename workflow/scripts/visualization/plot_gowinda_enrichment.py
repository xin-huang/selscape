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
from matplotlib.ticker import MaxNLocator

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

plot_title = snakemake.params['title']

# Read Gowinda enrichment results
try:
    data = pd.read_csv(
        snakemake.input['enrichment'],
        sep='\t',
        header=0,
    )
except pd.errors.EmptyDataError:
    data = pd.DataFrame()

if data.empty:
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, 'No results', ha='center', va='center', fontsize=16, color='#666')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    plt.savefig(snakemake.output['plot'], dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Sort by p_adjusted (ascending), then by genes_found (descending)
plot_data = data.sort_values(['p_adjusted', 'genes_found'], ascending=[True, False]).head(20).copy()

plot_data['short_description'] = plot_data['description'].apply(
    lambda x: x[:47] + "..." if len(x) > 50 else x
)

# Color mapping
def get_color(p_val):
    if p_val < 0.001:
        return '#FF0000'  # Red
    elif p_val < 0.01:
        return '#800080'  # Purple
    elif p_val < 0.05:
        return '#0000FF'  # Blue
    else:
        return '#808080'  # Gray

plot_data['color'] = plot_data['p_adjusted'].apply(get_color)

# custom legend
legend_elements = [
    Patch(facecolor='#FF0000', label='< 0.001'),
    Patch(facecolor='#800080', label='< 0.01'),
    Patch(facecolor='#0000FF', label='< 0.05'),
    Patch(facecolor='#808080', label='≥ 0.05')
]

# Create plot
fig, ax = plt.subplots(figsize=(12, max(8, len(plot_data) * 0.4)))

plot_data_display = plot_data[::-1]

ax.barh(range(len(plot_data_display)), plot_data_display['genes_found'], color=plot_data_display['color'])
ax.set_yticks(range(len(plot_data_display)))
ax.set_yticklabels(plot_data_display['short_description'], fontsize=9)
ax.set_xlabel('Count', fontsize=11)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_title(f'{plot_title} (Top {len(plot_data)} terms)', fontsize=14, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)
ax.legend(handles=legend_elements, title='FDR adjusted p-value', loc='lower right')


plt.tight_layout()
plt.savefig(snakemake.output['plot'], dpi=300, bbox_inches='tight')
plt.close()
