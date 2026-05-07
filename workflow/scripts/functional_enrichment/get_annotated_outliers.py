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


import pandas as pd
import os
import sys

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

outlier_file = snakemake.input.outliers
annotation_files = snakemake.input.annotation
output_file = snakemake.output.annotated_outliers

with open(output_file, "w"):
    pass

outlier_df = pd.read_csv(outlier_file, sep="\t")

if outlier_df.empty:
    log_fh.close()
    sys.exit(0)

multianno_dict = {
    os.path.basename(f).split(".")[0].removeprefix("chr"): f
    for f in annotation_files
}

outlier_df["CHR"] = outlier_df["CHR"].astype(str).str.removeprefix("chr")
outlier_df["BP"] = outlier_df["BP"].astype(int)

filtered_variants = []
for chrom, group in outlier_df.groupby("CHR"):
    print(f"Processing chromosome {chrom}...")
    if chrom not in multianno_dict:
        print(f"No annotation file for chromosome: {chrom}")
        continue

    print(f"Reading annotation file: {multianno_dict[chrom]}")
    multianno_df = pd.read_csv(multianno_dict[chrom], sep="\t")
    multianno_df["Chr"] = multianno_df["Chr"].astype(str).str.removeprefix("chr")
    multianno_df["Start"] = multianno_df["Start"].astype(int)
    multianno_df["End"] = multianno_df["End"].astype(int)

    variants = multianno_df[multianno_df["Start"].isin(set(group["BP"]))]
    print(f"Filtered {len(variants)} variants for chromosome {chrom}.")
    if not variants.empty:
        filtered_variants.append(variants.iloc[:, :11].copy())

if filtered_variants:
    result_df = pd.concat(filtered_variants, ignore_index=True)
    result_df["Chr"] = result_df["Chr"].astype(int)
    result_df = result_df.sort_values(by=["Chr", "Start", "End"])
    result_df.to_csv(output_file, sep="\t", index=False)

log_fh.close()
