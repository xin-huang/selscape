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


library(DT)
library(data.table)
library(htmlwidgets)

dt <- fread(snakemake@input[[1]], header = TRUE)
tbl <- DT::datatable(
  dt,
  filter = "top", 
  extensions = "Buttons",
  options = list(
    pageLength = 50, scrollX = TRUE,
    dom = "Bfrtip",
    buttons = c("copy", "csv", "excel")
  ),
  rownames = FALSE
)

saveWidget(tbl, snakemake@output[[1]], selfcontained = TRUE)
