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
from upsetplot import UpSet, from_contents
 
log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh
 
gene_files = snakemake.input.genes
population_groups = snakemake.params.population_groups or {}
max_intersections = int(snakemake.params.max_intersections)
plot_title = snakemake.params.title
 
BAR = "#0072B2"
FALLBACK_COLORS = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#D55E00"]
 
 
def no_results(message):
    with open(snakemake.output.table, "w") as out:
        out.write("groups\tdegree\tn_genes\n")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.55, "No results", ha="center", va="center", fontsize=16, color="#666")
    ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=9, color="#999")
    ax.axis("off")
    plt.savefig(snakemake.output.plot, bbox_inches="tight")
    plt.close()
    print(f"no results: {message}")
    log_fh.close()
    sys.exit(0)
 
 
group_of = {
    population: group
    for group, cfg in population_groups.items()
    for population in (cfg or {}).get("populations", [])
}
 
group_genes = {}
ungrouped = set()
for path in gene_files:
    population = os.path.basename(path).rsplit(".focal_", 1)[1].split(".")[0]
    group = group_of.get(population)
    if group is None:
        ungrouped.add(population)
        continue
    with open(path) as handle:
        genes = {g for g in (line.strip() for line in handle) if g and g != "Gene"}
    group_genes.setdefault(group, set()).update(genes)
 
if ungrouped:
    print("ignored, not in population_groups:", ", ".join(sorted(ungrouped)))
 
groups = [group for group in population_groups if group_genes.get(group)]
if len(groups) < 2:
    no_results("fewer than two population groups have candidates")
 
memberships = from_contents({group: group_genes[group] for group in groups})
intersections = memberships.index.value_counts()
 
with open(snakemake.output.table, "w") as out:
    out.write("groups\tdegree\tn_genes\n")
    for pattern, count in intersections.items():
        members = [group for group, present in zip(groups, pattern) if present]
        out.write(f"{'&'.join(members)}\t{len(members)}\t{count}\n")
 
upset = UpSet(
    memberships,
    subset_size="count",
    sort_by="cardinality",
    sort_categories_by=None,
    max_subset_rank=max_intersections,
    show_counts=True,
    facecolor=BAR,
    element_size=None,
)
for i, group in enumerate(groups):
    color = (population_groups[group] or {}).get("color") or FALLBACK_COLORS[i % len(FALLBACK_COLORS)]
    upset.style_categories(group, bar_facecolor=color)
 
fig = plt.figure(figsize=(max(6.0, 0.55 * min(len(intersections), max_intersections) + 3.0), 4.5), dpi=200)
upset.plot(fig=fig)
fig.suptitle(plot_title, fontsize=11, fontweight="bold")
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
 
print(f"{len(intersections)} intersections, showing at most {max_intersections}")
log_fh.close()
