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
    plt.savefig(snakemake.output['count_plot'], dpi=300, bbox_inches='tight')
    plt.savefig(snakemake.output['qscore_plot'], dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Prepare plot data
plot_data = data.sort_values('p_adjusted').head(20).copy()
plot_data['qscore'] = -np.log10(plot_data['p_adjusted'])
plot_data['short_description'] = plot_data['description'].apply(
    lambda x: x[:47] + "..." if len(x) > 50 else x
)

#  Color mapping based on p-value
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

# Function to create enrichment plot
def create_plot(data_sorted, value_col, xlabel, output_file):
    fig, ax = plt.subplots(figsize=(12, max(8, len(data_sorted) * 0.4)))
    ax.barh(range(len(data_sorted)), data_sorted[value_col], color=data_sorted['color'])
    ax.set_yticks(range(len(data_sorted)))
    ax.set_yticklabels(data_sorted['short_description'], fontsize=9)
    ax.set_xlabel(xlabel, fontsize=11)
    plt.suptitle(f'{plot_title} (Top {len(data_sorted)} terms)',
                fontsize=14, fontweight='bold', x=0.45)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.3)
    ax.legend(handles=legend_elements, title='p.adjust', loc='lower right')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# Create count plot
create_plot(
    plot_data.sort_values('genes_found', ascending=True),
    'genes_found',
    'Count',
    snakemake.output['count_plot']
)

# Create qscore plot
create_plot(
    plot_data.sort_values('qscore', ascending=True),
    'qscore',
    'qscore (-log10(p.adjusted))',
    snakemake.output['qscore_plot']
)
