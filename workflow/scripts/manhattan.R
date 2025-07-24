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



suppressMessages(library(qqman))

logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile, type = "message")

input_file <- snakemake@input[["scores"]]
score_column <- snakemake@params[["score_column"]]
cutoff <- as.numeric(snakemake@params[["cutoff"]])
use_absolute <- tolower(snakemake@params[["use_absolute"]]) %in% c("true", "1", "yes")
plot_width <- as.numeric(snakemake@params[["width"]])
plot_height <- as.numeric(snakemake@params[["height"]])
output_table <- snakemake@output[["candidates"]]
color1 <- snakemake@params[["color1"]]
color2 <- snakemake@params[["color2"]]
output_plot <- snakemake@output[["plot"]]

data <- read.table(input_file, header = TRUE)

score_used_col <- if (use_absolute) {
  abs_col <- paste0("abs_", score_column)
  data[[abs_col]] <- abs(data[[score_column]])
  abs_col
} else {
  score_column
}

data_sorted <- data[order(-data[[score_used_col]]), ]

top_n <- as.integer(nrow(data) * cutoff)
if (top_n < 1) stop("cutoff too small; no candidates selected.")
top_candidates <- data_sorted[1:top_n, ]
write.table(top_candidates, file = output_table, quote = FALSE, sep = "\t", row.names = FALSE)

threshold <- min(top_candidates[[score_used_col]])

png(output_plot, width = plot_width, height = plot_height, units = "px")
manhattan(data, p = score_used_col, logp = FALSE, genomewideline = threshold,
          suggestiveline = FALSE, col = c(color1, color2),
          ylab = if (use_absolute) paste0("|", score_column, "|") else score_column)
dev.off()
