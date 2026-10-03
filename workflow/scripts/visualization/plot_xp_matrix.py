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


import os
import sys
import matplotlib
 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle
 
log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh
 
gene_files = snakemake.input.genes
population_groups = snakemake.params.population_groups or {}
plot_title = snakemake.params.title
 
FALLBACK_COLORS = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#D55E00"]
 
 
def no_results(message):
    with open(snakemake.output.table, "w") as out:
        out.write("focal\treference\tn_genes\n")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.55, "No results", ha="center", va="center", fontsize=16, color="#666")
    ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=9, color="#999")
    ax.axis("off")
    plt.savefig(snakemake.output.plot, bbox_inches="tight")
    plt.close()
    print(f"no results: {message}")
    log_fh.close()
    sys.exit(0)
 
 
counts = {}
for path in gene_files:
    name = os.path.basename(path)
    pop1, pop2 = name.split(".", 1)[0].split("_", 1)
    focal = name.rsplit(".focal_", 1)[1].split(".")[0]
    reference = pop2 if focal == pop1 else pop1
    with open(path) as handle:
        genes = {g for g in (line.strip() for line in handle) if g and g != "Gene"}
    counts[(focal, reference)] = len(genes)
 
present = {pop for pair in counts for pop in pair}
group_of = {
    population: group
    for group, cfg in population_groups.items()
    for population in (cfg or {}).get("populations", [])
}
 
populations = [
    population
    for group in population_groups
    for population in (population_groups[group] or {}).get("populations", [])
    if population in present
]
populations += sorted(present - set(populations))
 
if len(populations) < 2:
    no_results("fewer than two populations")
 
index = {population: i for i, population in enumerate(populations)}
n = len(populations)
matrix = np.full((n, n), np.nan)
for (focal, reference), n_genes in counts.items():
    matrix[index[focal], index[reference]] = n_genes
 
with open(snakemake.output.table, "w") as out:
    out.write("focal\treference\tn_genes\n")
    for focal in populations:
        for reference in populations:
            if focal != reference:
                out.write(f"{focal}\t{reference}\t{counts.get((focal, reference), 0)}\n")
 
groups = [group for group in population_groups if group in set(group_of.get(p) for p in populations)]
group_color = {
    group: (population_groups[group] or {}).get("color") or FALLBACK_COLORS[i % len(FALLBACK_COLORS)]
    for i, group in enumerate(groups)
}
 
vmax = np.nanmax(matrix) if np.isfinite(matrix).any() else 1
cmap = matplotlib.colormaps["Blues"].copy()
cmap.set_bad("#F2F4F6")
 
side = max(5.0, 0.34 * n + 2.6)
fig, ax = plt.subplots(figsize=(side + 1.2, side), dpi=200)
image = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=max(vmax, 1), aspect="equal")
 
# white gaps between group blocks
for i in range(1, n):
    if group_of.get(populations[i]) != group_of.get(populations[i - 1]):
        ax.axhline(i - 0.5, color="#FFFFFF", linewidth=2.5)
        ax.axvline(i - 0.5, color="#FFFFFF", linewidth=2.5)
 
if n <= 12:
    for row in range(n):
        for column in range(n):
            if np.isfinite(matrix[row, column]):
                value = int(matrix[row, column])
                text_color = "#FFFFFF" if value > max(vmax, 1) * 0.6 else "#333333"
                ax.text(column, row, str(value), ha="center", va="center", fontsize=7.5, color=text_color)
 
label_size = 8 if n <= 20 else 7
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(populations, rotation=0, fontsize=label_size, color="#333333")
ax.set_yticklabels(populations, fontsize=label_size, color="#333333")
 
spans = []
for i, population in enumerate(populations):
    group = group_of.get(population)
    if spans and spans[-1][0] == group:
        spans[-1][2] = i
    else:
        spans.append([group, i, i])
 
BAND, GAP = 0.16, 0.14
edge = 0.5 + GAP + BAND
for group, start, end in spans:
    if group is None:
        continue
    color = group_color[group]
    ax.add_patch(Rectangle((start - 0.5, -edge), end - start + 1, BAND, color=color, clip_on=False))
    ax.add_patch(Rectangle((-edge, start - 0.5), BAND, end - start + 1, color=color, clip_on=False))
    ax.text((start + end) / 2, -edge - 0.12, group, ha="center", va="bottom",
            fontsize=label_size + 1, color="#333333", fontweight="bold")
 
ax.set_xlim(-edge - 0.05, n - 0.5)
ax.set_ylim(n - 0.5, -edge - 0.95)
ax.set_xlabel("Reference population", fontsize=9)
ax.set_ylabel("Focal population", fontsize=9)
ax.legend(
    handles=[Patch(facecolor=group_color[g], edgecolor="none", label=g) for g in groups],
    loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=min(len(groups), 6),
    frameon=False, fontsize=label_size, handlelength=1.1, handleheight=0.9, columnspacing=1.4,
)
ax.tick_params(length=0)
for spine in ax.spines.values():
    spine.set_visible(False)
 
bar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
bar.set_label("Outlier genes", fontsize=9)
bar.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
bar.ax.tick_params(labelsize=8)
bar.outline.set_visible(False)
 
ax.set_title(plot_title, fontsize=11, fontweight="bold", pad=10)
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
 
print(f"{n} populations, {len(counts)} filled cells, max {int(max(vmax, 0))} genes")
log_fh.close()
