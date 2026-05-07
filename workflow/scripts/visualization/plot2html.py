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

import base64

plot_file = snakemake.input.plot
output_html = snakemake.output.html
title = snakemake.params.title

# Read image and convert to base64
with open(plot_file, 'rb') as f:
    img_data = base64.b64encode(f.read()).decode('utf-8')

# Create HTML
html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            text-align: center;
        }}
        h1 {{
            font-size: 18px;
            font-weight: bold;
            margin-bottom: 20px;
            color: #333;
        }}
        img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            padding: 10px;
        }}
    </style>
</head>
<body>
    <h1>{title}</h1>
    <img src="data:image/png;base64,{img_data}" alt="{title}">
</body>
</html>
"""

with open(output_html, 'w') as f:
    f.write(html)

