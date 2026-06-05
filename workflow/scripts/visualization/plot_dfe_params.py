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
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import pandas as pd

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh

df = pd.read_csv(snakemake.input.data, sep="\t")

populations = snakemake.params.populations
mu_ylim = snakemake.params.mu_ylim
sigma_ylim = snakemake.params.sigma_ylim
population_groups = snakemake.params.population_groups

pop_colors = []
pop_labels = []
for pop in populations:
    label = next(
        (lbl for lbl, grp in population_groups.items() if pop in grp["populations"]),
        pop,
    )
    color = next(
        (
            grp["color"]
            for grp in population_groups.values()
            if pop in grp["populations"]
        ),
        None,
    )
    pop_labels.append(label)
    pop_colors.append(color)

# Fallback: auto colors wenn None
n = len(populations)
assert n > 0, "No populations provided; cannot generate DFE parameter plot."
auto_colors = [cm.tab20(i / n) for i in range(n)]
colors = [c if c is not None else auto_colors[i] for i, c in enumerate(pop_colors)]


assert len(df) == len(
    populations
), f"Row count mismatch: TSV has {len(df)} rows, expected {len(populations)}"
assert (
    list(df["Pop"]) == populations
), f"Population order mismatch: TSV has {list(df['Pop'])}, expected {populations}"

x = list(range(len(populations)))
x_tick_labels = [f"{dataset}_{pop}" for pop, dataset in zip(df["Pop"], df["Dataset"])]

fig, axs = plt.subplots(
    nrows=2, ncols=2, constrained_layout=True, figsize=(10, 4), dpi=350
)
gridspec = axs[0, 0].get_subplotspec().get_gridspec()
for a in axs[:, 1]:
    a.remove()


def plot_param(ax, param, ylim):
    lb, ub = f"{param}_lb", f"{param}_ub"
    ax.scatter(x, df[param].values, color=colors, zorder=2)
    ax.scatter(x, df[lb].values, marker="_", color="grey")
    ax.scatter(x, df[ub].values, marker="_", color="grey")
    for i in range(len(populations)):
        ax.plot(
            [i, i],
            [df[param].values[i], df[ub].values[i]],
            linestyle="dashed",
            color="grey",
            zorder=1,
        )
        ax.plot(
            [i, i],
            [df[lb].values[i], df[param].values[i]],
            linestyle="dashed",
            color="grey",
            zorder=1,
        )
    if param == "sigma":
        current = ax.get_ylim()
        ax.set_ylim(bottom=0, top=current[1])
    if ylim is not None:
        if isinstance(ylim, (list, tuple)) and len(ylim) == 2:
            ax.set_ylim(ylim[0], ylim[1])
        else:
            raise ValueError(
                f"Invalid ylim for {param}: {ylim}. Expected [min, max] or null."
            )
    ax.set_xticks(x, x_tick_labels, rotation=90)


plot_param(axs[0, 0], "mu", mu_ylim)
plot_param(axs[1, 0], "sigma", sigma_ylim)
axs[0, 0].set_ylabel(r"$\mu$")
axs[1, 0].set_ylabel(r"$\sigma$")

# Legend — unique groups from pop_labels
subfig = fig.add_subfigure(gridspec[:, 1])
seen = {}
for label, color in zip(pop_labels, colors):
    if label not in seen:
        seen[label] = color
handles = [mpatches.Patch(color=color, label=label) for label, color in seen.items()]
subfig.legend(handles=handles, fontsize=8, handlelength=1.5, loc="upper left")

fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0, wspace=0.1)
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()

log_fh.close()
