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


import os
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset

matplotlib.use("Agg")

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh

population = snakemake.params.population
ref_genome = snakemake.params.ref_genome
tracks_cfg = snakemake.params.tracks  # list of dicts from rule params

STATS = {
    "iHS":  "iHS",
    "nSL":  "nSL",
    "B1":   r"$\beta^{(1)}$",
    "mtjd": r"$D_{Tm}$",
    "wtjd": r"$D_{Tw}$",
}


def load_score_file(score_file, score_column):
    df = pd.read_csv(score_file, sep="\t")
    df["chr"] = "chr" + df["CHR"].astype(str)
    df["pos"] = df["BP"].astype(int) if "BP" in df.columns else df["window_start"].astype(int)
    df["score"] = df[score_column].astype(float)
    df = df[~df["chr"].isin(["chrX", "chrY"])]
    df = df.dropna(subset=["score"])
    return df[["chr", "pos", "score"]]


# Load all track data from params
track_data = []
for t in tracks_cfg:
    score_file = getattr(snakemake.input, t["file"])
    df = load_score_file(score_file, t["score_col"])
    vmin, vmax = df["score"].min(), df["score"].max()
    track_data.append((t["name"], df, t["r_range"], t["color"], vmin, vmax))

# Load chromosome bed
chr_bed_file = snakemake.input.chr_bed
cytoband_file = snakemake.input.cytoband
# chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset(ref_genome)
chr_bed_df = pd.read_csv(chr_bed_file, sep="\t", header=None)
autosomes = [f"chr{i}" for i in range(1, 23)]
chr_bed_df = chr_bed_df[chr_bed_df[0].isin(autosomes)]

temp_dir = os.path.dirname(snakemake.output.plot)
os.makedirs(temp_dir, exist_ok=True)
filtered_chr_bed = os.path.join(temp_dir, f"{population}_chr.bed")
chr_bed_df.to_csv(filtered_chr_bed, sep="\t", header=False, index=False)

# Setup circos
circos = Circos.initialize_from_bed(filtered_chr_bed, start=-80, end=260, space=3, endspace=False)
circos.add_cytoband_tracks((81, 85), cytoband_file)
first_sector_name = circos.sectors[0].name

for sector in circos.sectors:
    axis_track = sector.add_track((81, 85), r_pad_ratio=0.0)
    axis_track.axis(lw=0.6)
    sector.text(sector.name, r=88, size=8)

    for stat_name, df, (r_inner, r_outer), color, vmin, vmax in track_data:
        track = sector.add_track((r_inner, r_outer), r_pad_ratio=0.0)
        track.axis(lw=0.3)
        chr_data = df[df["chr"] == sector.name].sort_values("pos")
        if chr_data.empty:
            continue
        track.line(
            chr_data["pos"].values,
            chr_data["score"].values,
            color=color,
            lw=0.3,
            vmin=vmin - 0.5,
            vmax=vmax + 0.5,
        )
        if sector.name == first_sector_name:
            track.yticks(
                [vmin, vmax],
                vmin=vmin - 0.5,
                vmax=vmax + 0.5,
                labels=[f"{round(vmin, 1):.1f}", f"{round(vmax, 1):.1f}"],
                side="left",
                tick_length=1.2,
                label_size=6,
                label_margin=0.5,
            )

    if sector.name == first_sector_name:
        circos.text("Chr", r=axis_track.r_center, deg=-90, size=8)
        for stat_name, df, (r_inner, r_outer), color, vmin, vmax in track_data:
            circos.text(
                STATS.get(stat_name, stat_name),
                r=(r_inner + r_outer) / 2,
                deg=-90,
                color=color,
                size=6,
                va="center",
            )

circos.text(population, size=12, r=0, weight="bold")
fig = circos.plotfig()
fig.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")
plt.close()

log_fh.close()
