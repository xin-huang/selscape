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

import pandas as pd

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

input_file = snakemake.input.outliers
output_file = snakemake.output.scores

score_column = snakemake.params.score_column
pop1_sign = str(snakemake.params.pop1_sign).strip().lower()

pair = snakemake.wildcards.pair
focal = snakemake.wildcards.focal
pop1_name, pop2_name = pair.split("_", 1)

assert pop1_sign in ("positive", "negative"), f"pop1_sign must be positive or negative, got {pop1_sign!r}"
assert focal in (pop1_name, pop2_name), f"focal {focal!r} is not part of pair {pair!r}"


def write_empty(path):
    """manhattan.R writes a zero-byte file when it has nothing to report; emit a
    header so downstream pd.read_csv sees an empty table instead of an error."""
    with open(path, "w") as out:
        out.write(f"SNP\tCHR\tBP\t{score_column}\n")

try:
    outliers = pd.read_csv(input_file, sep="\t")
except pd.errors.EmptyDataError:
    print(f"{input_file} is empty -- writing a header-only output", file=log_fh)
    write_empty(output_file)
    log_fh.close()
    sys.exit(0)

if score_column not in outliers.columns:
    raise KeyError(
        f"score column {score_column!r} not found in {input_file}. "
        f"Columns present: {list(outliers.columns)}. "
        f"manhattan.R keeps the signed score alongside the abs_ column."
    )

scores = pd.to_numeric(outliers[score_column], errors="coerce")
keep_positive = (pop1_sign == "positive") == (focal == pop1_name)
mask = scores > 0 if keep_positive else scores < 0

outliers[mask].to_csv(output_file, sep="\t", index=False)

print(
    f"pair={pair} focal={focal} pop1_sign={pop1_sign} score_column={score_column} "
    f"keeping {'positive' if keep_positive else 'negative'} scores",
    file=log_fh,
)
print(
    f"{int(mask.sum())} of {len(outliers)} outliers assigned to {focal}, "
    f"{int((scores == 0).sum())} dropped (score exactly 0), "
    f"{int(scores.isna().sum())} dropped (not numeric)",
    file=log_fh,
)

log_fh.close()
