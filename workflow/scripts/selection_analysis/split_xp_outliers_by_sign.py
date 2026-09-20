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
output_pop1 = snakemake.output.pop1
output_pop2 = snakemake.output.pop2

score_column = snakemake.params.score_column
pop1_sign = str(snakemake.params.pop1_sign).strip().lower()

pair = snakemake.wildcards.pair
pop1_name, pop2_name = pair.split("_", 1)

assert pop1_sign in ("positive", "negative"), f"pop1_sign must be positive or negative, got {pop1_sign!r}"


def write_empty(path):
    with open(path, "w"):
        pass


try:
    outliers = pd.read_csv(input_file, sep="\t")
except pd.errors.EmptyDataError:
    print(f"{input_file} is empty -- writing empty outputs", file=log_fh)
    write_empty(output_pop1)
    write_empty(output_pop2)
    log_fh.close()
    sys.exit(0)

if score_column not in outliers.columns:
    raise KeyError(
        f"score column {score_column!r} not found in {input_file}. "
        f"Columns present: {list(outliers.columns)}. "
        f"manhattan.R keeps the signed score alongside the abs_ column."
    )

scores = pd.to_numeric(outliers[score_column], errors="coerce")

if pop1_sign == "positive":
    mask_pop1 = scores > 0
    mask_pop2 = scores < 0
else:
    mask_pop1 = scores < 0
    mask_pop2 = scores > 0

n_zero = int((scores == 0).sum())
n_nan = int(scores.isna().sum())

outliers[mask_pop1].to_csv(output_pop1, sep="\t", index=False)
outliers[mask_pop2].to_csv(output_pop2, sep="\t", index=False)

print(
    f"pair={pair} pop1_sign={pop1_sign} score_column={score_column}",
    file=log_fh,
)
print(
    f"{int(mask_pop1.sum())} outliers assigned to pop1={pop1_name}, "
    f"{int(mask_pop2.sum())} to pop2={pop2_name}, "
    f"{n_zero} dropped (score exactly 0), {n_nan} dropped (not numeric)",
    file=log_fh,
)

log_fh.close()
