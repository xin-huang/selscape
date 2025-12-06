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
import numpy as np


log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

file_list = snakemake.input.scores
usecols = snakemake.params.usecols
is_abs = snakemake.params.is_abs
col_names = snakemake.params.col_names

# Read all files and concatenate into a single DataFrame
df = pd.concat(
    [pd.read_csv(file, delim_whitespace=True, header=None, usecols=usecols) for file in file_list],
    ignore_index=True
)

if is_abs: df.iloc[:, 2] = df.iloc[:, 2].abs()
df_sorted = df.sort_values(by=df.columns[2], ascending=False)
df_sorted.columns = col_names
df_sorted.to_csv(snakemake.output.scores, sep="\t", index=False, header=True)
