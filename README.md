# DISCO Calculator

Web calculator for **DISCO (Distance of Covariance)** and **DM (Mahalanobis
distance)**, measures of homeostatic dysregulation that reveal organ-system
interconnections underlying mortality and disease risk.

Live: https://jerryhaom.github.io/DISCO-Calculator/

## What it computes

DISCO quantifies how much one person's biomarker profile perturbs the
correlation structure of a young, healthy reference population:

```
DISCO_i = log( n_ref² · Σ_j,k  w_jk · (PCC_ref_jk − PCC_ref+i_jk)² )
```

where `PCC_ref` is the reference correlation matrix, `PCC_ref+i` is the
correlation matrix after adding individual *i*, `w` is an age-based weight
matrix, and `n_ref` is the reference sample size.

DM (Mahalanobis distance) quantifies the standardized deviation of the
biomarker vector from the reference mean:

```
DM_i = log( √( (x_i − μ)ᵀ S⁻¹ (x_i − μ) ) )
```

where `μ` and `S` are the mean vector and covariance matrix of the reference
population.

## How values are estimated — no raw data required

Both metrics are computed entirely from a **standard reference matrix** of
aggregate statistics:

- `μ` — reference mean vector (10 values)
- `ss` — cross-product matrix, from which the covariance `S = ss/(n_ref−1)`
  and correlation matrix `PCC` are derived (10×10)
- `w` — age-based weight matrix (10×10)
- `n_ref` — reference sample size

No individual-level reference data is stored or required. Given a person's own
10 biomarker values, DISCO is obtained by the rank-one update
`ss_new = ss + (n_ref/(n_ref+1))·δδᵀ` (where `δ = x − μ`), and DM by the usual
quadratic form `(x−μ)ᵀS⁻¹(x−μ)`. The reference matrix is derived from UK
Biobank participants aged ≤40 (n = 5,282) and is freely redistributable.

## Inputs (10 biomarkers)

CRP, glucose, total cholesterol, HDL cholesterol, creatinine, uric acid, white
blood cell count, urea, red blood cell count, albumin.

CRP is entered as raw mg/L and log-transformed (`log(CRP+1)`) automatically.

## Interpretation

- Higher DISCO / DM = greater homeostatic dysregulation.
- Because absolute values are cohort-specific (population structure + assay
  calibration), the calculator also reports a **cohort-adjusted z-score** and
  percentile relative to a selected reference cohort (UK Biobank, NHANES,
  CHARLS, CLHLS, or RuLAS).

In validation cohorts, this standard reference matrix preserves 96–97% of
predictive power across cohorts.

## Reference

Hao M, et al. *Distance of Covariance (DISCO): a novel measure of network
homeostatic dysregulation revealing organ system interconnections underlying
mortality and disease risk.* (manuscript)

## Contact

haombio@gmail.com
