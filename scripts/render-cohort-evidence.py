#!/usr/bin/env python3
"""Render the static, offline-compatible evidence tables from aggregate JSON."""
import argparse
import json
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NAMES = {"age": "age", "sex": "sex", "bmi": "BMI", "race": "ethnicity/race",
         "education": "education", "smoking": "smoking", "alcohol": "alcohol use",
         "physical_activity": "physical activity"}


def render(data):
    cohorts = data["cohorts"]
    headers = ''.join(f'<th scope="col">{escape(c["label"])}<span class="row-subtitle">'
                      f'{c["marker_count"]} / 10 markers</span></th>' for c in cohorts)
    heat_rows = []
    for marker in data["markers"]:
        cells = []
        for c in cohorts:
            present = marker in c["markers"]
            status = "Available" if present else "Unavailable in analyzed wave"
            css = "heat-cell" + ("" if present else " unavailable")
            title = escape(f'{c["label"]} · {marker}: {status}')
            cells.append(f'<td><span class="{css}" aria-label="{status}" title="{title}">'
                         f'{"✓" if present else "—"}</span></td>')
        heat_rows.append(f'<tr><th scope="row">{escape(marker)}</th>{"".join(cells)}</tr>')

    rows = []
    for c in cohorts:
        cells = []
        for metric in ["DISCO", "DM"]:
            e = c["estimates"][metric]
            for value, lo, hi, places in [("c_index", "c_low", "c_high", 3), ("hr", "hr_low", "hr_high", 2)]:
                cells.append(f'<td data-metric="{metric}" data-estimate="{value}">{e[value]:.{places}f}'
                             f'<span class="ci">({e[lo]:.{places}f}–{e[hi]:.{places}f})</span></td>')
        rows.append(f'<tr data-cohort="{c["cohort"]}"><th scope="row">{escape(c["label"])}'
                    f'<span class="row-subtitle">{c["marker_count"]} markers</span></th>'
                    f'<td>{c["n"]:,}<span class="ci">{c["events"]:,} deaths</span></td>{"".join(cells)}</tr>')
    covariates = ''.join(f'<li><strong>{escape(c["label"])}:</strong> '
                         f'{", ".join(NAMES[x] for x in c["covariates"])}.</li>' for c in cohorts)
    waves = '; '.join(f'{c["label"]}: {c["wave"]}' for c in cohorts)
    ph = []
    for c in cohorts:
        for metric, e in c["estimates"].items():
            if e["ph_test_p"] < .05:
                ph.append(f'{c["label"]} {metric} (p = {e["ph_test_p"]:.3g})')
    return f'''  <section class="cohort-evidence" aria-labelledby="evidenceHeading">
    <div class="evidence-header"><div><div class="eyebrow">Cohort evidence</div><h2 id="evidenceHeading">Biomarker coverage &amp; mortality prediction</h2><p>Availability of the calculator's 10 biomarkers in the analyzed cohort waves. Performance uses each cohort's available panel and the fixed UK Biobank reference matrix.</p></div></div>
    <div class="evidence-legend"><span class="legend-item"><span class="legend-square" aria-hidden="true"></span>Available</span><span class="legend-item"><span class="legend-square unavailable" aria-hidden="true"></span>Unavailable in analyzed wave</span></div>
    <div class="table-scroll" role="region" aria-label="Biomarker availability by cohort" tabindex="0"><table class="evidence-table availability-table"><thead><tr><th scope="col">Biomarker</th>{headers}</tr></thead><tbody>
      {chr(10).join(heat_rows)}
    </tbody></table></div>
    <p class="evidence-footnote">CHARLS lacks RBC and albumin in the analyzed panel; RuLAS wave 2 lacks urea. Their results use 8- and 9-marker subsets, respectively. “Available” describes cohort-level collection, not complete observations for every participant. CLHLS uses hs-CRP. Tables scroll horizontally on small screens.</p>
    <h3 class="performance-heading">All-cause mortality · Fully adjusted models</h3>
    <div class="table-scroll" role="region" aria-label="Mortality prediction performance" tabindex="0"><table class="evidence-table performance-table"><thead><tr><th rowspan="2" scope="col">Cohort / panel</th><th rowspan="2" scope="col">Participants / deaths</th><th colspan="2" class="metric-group" scope="colgroup">DISCO</th><th colspan="2" class="metric-group" scope="colgroup">DM</th></tr><tr><th scope="col">C-index (95% CI)</th><th scope="col">HR per SD (95% CI)</th><th scope="col">C-index (95% CI)</th><th scope="col">HR per SD (95% CI)</th></tr></thead><tbody id="performanceRows">
      {chr(10).join(rows)}
    </tbody></table></div>
    <p class="evidence-footnote"><strong>C-index is for the score + covariates model</strong>, fitted and evaluated in the same cohort (apparent performance). HR is per 1 SD within that cohort's analysis sample. DISCO uses its natural-log scale; DM is untransformed. Covariates and panel size vary by cohort, so these estimates are not a direct ranking of cohorts or panels.</p>
    <p class="evidence-footnote">All models adjust for age, sex, BMI, education, smoking and alcohol use; UKB and NHANES also include ethnicity/race, and UKB includes physical activity. Some score associations vary with follow-up time; their single HR should not be interpreted as constant over time.</p>
    <details class="evidence-methods"><summary>Analysis methods &amp; cohort-specific covariates</summary><div id="evidenceMethods">
      <p><strong>Model and samples.</strong> Separate Cox models with Efron ties for DISCO and DM, using identical complete-case participants within each cohort. Only binary mortality status (0/1) and positive follow-up are accepted; lost-to-follow-up codes are excluded. No additional imputation is performed on the supplied cleaned biomarker inputs. Age and BMI are continuous; other covariates are categorical. “Fully adjusted” refers to the available adjustment set listed below, not every potential confounder.</p>
      <ul>{covariates}</ul>
      <p><strong>Analyzed waves.</strong> {escape(waves)}.</p>
      <p><strong>Scoring.</strong> Biomarkers are harmonized to the calculator's units, CRP is transformed with ln(CRP + 1), and the fixed UKB reference is applied. CHARLS and RuLAS use matching submatrices with the original weights, without renormalization. Their subset scores are distinct from full 10-marker scores. The calculator's relative z-scores and percentiles also use the matching cohort panel.</p>
      <p><strong>Estimation.</strong> Harrell's C is calculated from each full model's linear predictor. HR and C-index intervals are 95% Wald intervals; C-index bounds are restricted to [0, 1]. These models are unweighted, including NHANES; its complex survey design is not incorporated. Estimates have not been corrected for optimism or validated as a fixed mortality prediction model in held-out data. The calculator does not produce an absolute mortality risk.</p>
      <p><strong>Proportional-hazards check.</strong> Score-specific Schoenfeld-residual tests give p &lt; 0.05 for {escape('; '.join(ph))}. These results indicate departures from proportional hazards; the displayed HRs summarize the fitted associations and warrant time-varying-effect sensitivity analyses. Other test p-values are in the downloadable results.</p>
      <p><strong>Provenance.</strong> Recomputed {escape(data['generated_utc'][:10])} using the existing cohort biomarker and survival files and available covariates. <a href="data/cohort-evidence.json" download>Download aggregate results (JSON)</a> · <a href="scripts/rebuild-cohort-evidence.R">Reproduction script</a>. Only aggregate results are published; participant-level inputs remain outside this repository.</p>
    </div></details>
  </section>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Fail if the embedded tables are stale')
    args = parser.parse_args()
    page = ROOT / 'index.html'
    html = page.read_text()
    start = html.index('  <section class="cohort-evidence"')
    end = html.index('  </section>', start) + len('  </section>')
    section = render(json.loads((ROOT / 'data/cohort-evidence.json').read_text()))
    updated = html[:start] + section + html[end:]
    if args.check:
        if updated != html:
            raise SystemExit('Cohort tables are stale; run scripts/render-cohort-evidence.py')
        print('Cohort tables match aggregate results.')
    else:
        page.write_text(updated)
        print('Rendered cohort evidence.')


if __name__ == '__main__':
    main()
