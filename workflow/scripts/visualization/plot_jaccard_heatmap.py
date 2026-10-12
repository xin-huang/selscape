# Copyright 2026 Xin Huang and Simon Chen
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
from itertools import combinations
import matplotlib
 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
 
log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh
 
gene_files = snakemake.input.genes
labels = snakemake.params.labels
population_groups = snakemake.params.population_groups or {}
plot_title = snakemake.params.title
 
 
def jaccard(a, b):
    union = a | b
    if len(union) == 0:
        return float("nan")
    return len(a & b) / len(union)
 
 
def no_results(plot_path, table_path, message):
    with open(table_path, "w") as out:
        out.write("")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.55, "No results", ha="center", va="center", fontsize=16, color="#666")
    ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=9, color="#999")
    ax.axis("off")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()
    print(f"no results ({plot_path}): {message}")
 
 
def plot_heatmap(matrix, title, output_path):
    n_rows = len(matrix.index)
    n_cols = len(matrix.columns)
    fig, ax = plt.subplots(figsize=(0.6 * n_cols + 3, 0.35 * n_rows + 2), dpi=300)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e8e8e8")
    masked = np.ma.masked_invalid(matrix.values.astype(float))
    im = ax.imshow(masked, vmin=0, vmax=1, cmap=cmap, aspect="auto")
 
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(matrix.index)
 
    for i in range(n_rows):
        for j in range(n_cols):
            val = matrix.values[i, j]
            if np.isnan(val):
                continue
            color = "white" if val < 0.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=8)
 
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard similarity")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
 
 
# (dataset, population) -> set(genes); labels are "species:dataset:population", same order as input
gene_sets = {}
datasets = []
species_of = {}
for path, label in zip(gene_files, labels):
    species, dataset, population = label.split(":", 2)
    species_of[dataset] = species
    if dataset not in datasets:
        datasets.append(dataset)
    with open(path) as handle:
        gene_sets[(dataset, population)] = {g for g in (line.strip() for line in handle) if g and g != "Gene"}
 
# dataset level: genes pooled (union) over all populations of a dataset
pooled = {dataset: set() for dataset in datasets}
for (dataset, population), genes in gene_sets.items():
    pooled[dataset].update(genes)
 
if len(datasets) < 2:
    no_results(snakemake.output.dataset_plot, snakemake.output.dataset_table, "fewer than two datasets")
else:
    matrix = pd.DataFrame(float("nan"), index=datasets, columns=datasets)
    for dataset in datasets:
        matrix.loc[dataset, dataset] = jaccard(pooled[dataset], pooled[dataset])
    for dataset_a, dataset_b in combinations(datasets, 2):
        value = jaccard(pooled[dataset_a], pooled[dataset_b])
        matrix.loc[dataset_a, dataset_b] = value
        matrix.loc[dataset_b, dataset_a] = value
    matrix.to_csv(snakemake.output.dataset_table, sep="\t")
    plot_heatmap(matrix, f"{plot_title}: datasets", snakemake.output.dataset_plot)
    print(f"dataset level: {len(datasets)} datasets")
 
# population level: same population compared across datasets, no pooling;
# only populations present in at least two datasets
present = {}
for dataset, population in gene_sets:
    present.setdefault(population, []).append(dataset)
shared = {population for population, found in present.items() if len(found) >= 2}
 
populations = [
    population
    for group in population_groups
    for population in (population_groups[group] or {}).get("populations", [])
    if population in shared
]
populations += sorted(shared - set(populations))
 
pairs = [
    (dataset_a, dataset_b)
    for dataset_a, dataset_b in combinations(datasets, 2)
    if species_of[dataset_a] == species_of[dataset_b]
    and any((dataset_a, p) in gene_sets and (dataset_b, p) in gene_sets for p in populations)
]
 
if not pairs:
    no_results(snakemake.output.population_plot, snakemake.output.population_table,
               "no population shared between datasets")
else:
    columns = [f"{dataset_a} vs {dataset_b}" for dataset_a, dataset_b in pairs]
    table = pd.DataFrame(float("nan"), index=populations, columns=columns)
    for population in populations:
        for (dataset_a, dataset_b), column in zip(pairs, columns):
            if (dataset_a, population) in gene_sets and (dataset_b, population) in gene_sets:
                table.loc[population, column] = jaccard(
                    gene_sets[(dataset_a, population)], gene_sets[(dataset_b, population)]
                )
    table.to_csv(snakemake.output.population_table, sep="\t")
    plot_heatmap(table.T, f"{plot_title}: populations", snakemake.output.population_plot)
    print(f"population level: {len(populations)} populations x {len(pairs)} dataset pairs")
 
log_fh.close()
