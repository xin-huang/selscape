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


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

outlier_file = snakemake.input.outliers
annotation_files = snakemake.input.annotation
output_file = snakemake.output.annotated_candidates

with open(output_file, "w"):
    pass

outlier_df = pd.read_csv(outlier_file, sep="\t")

#value_col_name = outlier_df.columns[3]

multianno_dict = {
    os.path.basename(f).split("chr")[1].split(".")[0]: f for f in annotation_files
}

candidate_positions = set()
#position_to_value = {}

for _, row in outlier_df.iterrows():
    chrom = str(row["CHR"]).removeprefix("chr")
    pos = int(row["BP"])
    #value = row[value_col_name]
    candidate_positions.add((chrom, pos))
    #position_to_value[(chrom, pos)] = value

filtered_variants = []

for chrom in set(chrom for chrom, _ in candidate_positions):
    print(f"Processing chromosome {chrom}...")

    if chrom not in multianno_dict:
        print(f"No annotation file for chromosome: {chrom}")
        continue

    multianno_file = multianno_dict[chrom]
    print(f"Reading annotation file: {multianno_file}")

    multianno_df = pd.read_csv(multianno_file, sep="\t")

    multianno_df["Chr"] = multianno_df["Chr"].astype(str).apply(lambda x: x.removeprefix("chr"))
    multianno_df["Start"] = multianno_df["Start"].astype(int)
    multianno_df["End"] = multianno_df["End"].astype(int)

    positions_for_chrom = {pos for c, pos in candidate_positions if c == chrom}

    variants_in_region = multianno_df[multianno_df["Start"].isin(positions_for_chrom)]

    print(f"Filtered {len(variants_in_region)} variants for chromosome {chrom}.")

    if not variants_in_region.empty:
        variants_in_region = variants_in_region.iloc[:, :11].copy()
        #variants_in_region[value_col_name] = variants_in_region["Start"].apply(
        #    lambda pos: position_to_value.get((chrom, pos), None)
        #)
        filtered_variants.append(variants_in_region)

if filtered_variants:
    result_df = pd.concat(filtered_variants, ignore_index=True)
    result_df["Chr"] = result_df["Chr"].astype(int)
    result_df = result_df.sort_values(by=["Chr", "Start", "End"], ascending=[True, True, True])
    result_df.to_csv(output_file, sep="\t", index=False)
