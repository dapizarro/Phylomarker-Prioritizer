"""Informe HTML."""
from __future__ import annotations

import html

import pandas as pd
import yaml

from .layout import OutputLayout


def create_html_report(
    layout: OutputLayout,
    runs: pd.DataFrame,
    warnings: pd.DataFrame,
    metrics: pd.DataFrame,
    scored: pd.DataFrame,
    config: dict,
) -> None:
    layout.report_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    eligible_count = int(
        scored["eligible"].sum()
    )

    excluded_count = int(
        (~scored["eligible"]).sum()
    )

    top_complete = scored.sort_values(
        "cell_occupancy",
        ascending=False,
    ).head(15)

    table_rows = "\n".join(
        (
            "<tr>"
            f"<td>{html.escape(str(row.gene_id))}</td>"
            f"<td>{row.cell_occupancy:.3f}</td>"
            f"<td>{row.parsimony_informative_sites}</td>"
            f"<td>{row.gap_fraction:.3f}</td>"
            f"<td>{row.composition_variability:.4f}</td>"
            "</tr>"
        )
        for row in top_complete.itertuples(
            index=False
        )
    )

    warning_count = len(warnings)

    report = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PhyloMarker Select report</title>
<style>
body {{
    max-width: 1100px;
    margin: 40px auto;
    padding: 0 24px;
    font-family: system-ui, sans-serif;
    color: #17202a;
    line-height: 1.5;
}}
h1, h2 {{
    color: #17324d;
}}
.cards {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 16px;
}}
.card {{
    border: 1px solid #d9e2ec;
    border-radius: 10px;
    padding: 16px;
    background: #f8fafc;
}}
.value {{
    font-size: 2rem;
    font-weight: 700;
}}
table {{
    border-collapse: collapse;
    width: 100%;
}}
th, td {{
    border-bottom: 1px solid #d9e2ec;
    padding: 8px;
    text-align: left;
}}
.notice {{
    border-left: 4px solid #d97706;
    padding: 12px;
    background: #fff7ed;
}}
pre {{
    white-space: pre-wrap;
    background: #f1f5f9;
    padding: 16px;
}}
</style>
</head>
<body>
<h1>PhyloMarker Select</h1>

<p>
Evolution-aware marker characterization and panel optimization.
PCA is exploratory and is not used as a biological quality ranking.
</p>

<div class="cards">
<div class="card">
<div class="value">{len(runs)}</div>
Validated BUSCO runs
</div>

<div class="card">
<div class="value">{len(metrics)}</div>
Aligned markers
</div>

<div class="card">
<div class="value">{eligible_count}</div>
Eligible markers
</div>

<div class="card">
<div class="value">{excluded_count}</div>
Excluded markers
</div>

<div class="card">
<div class="value">{warning_count}</div>
Validation warnings
</div>
</div>

<h2>Scientific warning</h2>

<div class="notice">
Marker suitability depends on taxonomic sampling and evolutionary objective.
High occupancy does not guarantee phylogenetic signal, and many informative
sites do not guarantee an unbiased or correct gene tree. Final panel quality
must be evaluated using independent phylogenetic analyses.
</div>

<h2>Most complete markers</h2>

<table>
<thead>
<tr>
<th>Gene</th>
<th>Cell occupancy</th>
<th>PIS</th>
<th>Gap fraction</th>
<th>Composition variability</th>
</tr>
</thead>
<tbody>
{table_rows}
</tbody>
</table>

<h2>Interpretation</h2>

<ul>
<li>PCA describes multivariate structure; it does not rank biological quality.</li>
<li>Singleton taxonomic groups are reported separately from replicated groups.</li>
<li>Eligibility thresholds are hard constraints and cannot be compensated by a high score.</li>
<li>Panel optimization penalizes redundant genes.</li>
<li>Random and single-criterion panels are generated as controls.</li>
</ul>

<h2>Resolved configuration</h2>

<pre>{html.escape(yaml.safe_dump(config, sort_keys=False))}</pre>
</body>
</html>
"""

    layout.report_index.write_text(
        report,
        encoding="utf-8",
    )
