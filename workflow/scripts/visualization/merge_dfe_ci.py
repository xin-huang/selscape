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
import pandas as pd
import re

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh
sys.stdout = log_fh

# get input files from Snakemake
populations   = snakemake.params.populations
bestfit_files = snakemake.input.bestfit_files
ci_files      = snakemake.input.ci_files
output_file   = snakemake.output.merged
datasets = snakemake.params.datasets

def parse_ci_file(ci_file):
    with open(ci_file, 'r') as f:
        lines = f.readlines()

    lower_bounds = None
    upper_bounds = None

    for line in reversed(lines):
        if line.startswith("Lower bounds"):
            match = re.search(r'\[(.*?)\]', line)
            if match:
                lower_bounds = [float(x) for x in match.group(1).split()]
        elif line.startswith("Upper bounds"):
            match = re.search(r'\[(.*?)\]', line)
            if match:
                upper_bounds = [float(x) for x in match.group(1).split()]
        if lower_bounds and upper_bounds:
            break

    return lower_bounds, upper_bounds


def parse_bestfit_file(bestfit_file):
    with open(bestfit_file, 'r') as f:
        for line in f:
            if line.startswith('# Converged results'):
                next(f)
                values = next(f).strip().split()
                log_mu    = float(values[1])
                log_sigma = float(values[2])
                misid     = float(values[3]) if len(values) > 3 else None
                return log_mu, log_sigma, misid
    return None, None, None

data = []
for pop, dataset, bestfit_file, ci_file in zip(populations, datasets, bestfit_files, ci_files):
    log_mu, log_sigma, misid = parse_bestfit_file(bestfit_file)
    lower_bounds, upper_bounds = parse_ci_file(ci_file)
    data.append({
        'Pop':      pop,
        'Dataset': dataset, 
        'mu':       log_mu,
        'mu_lb':    lower_bounds[0],
        'mu_ub':    upper_bounds[0],
        'sigma':    log_sigma,
        'sigma_lb': lower_bounds[1],
        'sigma_ub': upper_bounds[1],
        'misid':    misid,
        'misid_lb': lower_bounds[2] if len(lower_bounds) > 2 else None,
        'misid_ub': upper_bounds[2] if len(upper_bounds) > 2 else None,
    })

df = pd.DataFrame(data)
df.to_csv(output_file, sep='\t', index=False)
