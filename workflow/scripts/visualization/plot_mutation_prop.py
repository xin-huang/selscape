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


import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

input_file = snakemake.input.dfe_popt
plot_title = snakemake.params.title

with open(input_file, "r") as f:
    lines = f.readlines()

start_idx = None
for i, line in enumerate(lines):
    if line.strip().startswith("# Converged results"):
        start_idx = i + 1
        break

if start_idx is None:
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.text(0.5, 0.5, 'No results', ha='center', va='center', fontsize=16, color='#666')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)

for line in lines[start_idx:]:
    if line.strip() == "" or line.strip().startswith("#"):
        continue
    cols = line.strip().split()
    log_mu = float(cols[1])
    log_sigma = float(cols[2])
    break

ps = stats.lognorm.cdf([1, 10, 100], log_sigma, scale=np.exp(log_mu))
props = [ps[0], ps[1] - ps[0], ps[2] - ps[1], 1 - ps[2]]

fig = plt.figure(figsize=(5, 4))
plt.bar([0, 1, 2, 3], props, alpha=0.7, color="steelblue")
plt.ylabel("Proportion")
plt.xlabel("2N|s|")
plt.xticks(
    [0, 1, 2, 3],
    ["<1", "1–10", "10–100", ">100"],
)
plt.grid(alpha=0.3)
plt.title(plot_title, fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
