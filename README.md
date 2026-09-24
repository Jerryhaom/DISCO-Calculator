# DISCO Calculator

Web calculator for **DISCO (Distance of Covariance)**, a measure of network
homeostatic dysregulation that reveals organ-system interconnections underlying
mortality and disease risk.

Live: https://jerryhaom.github.io/DISCO-Calculator/

## What it computes

DISCO quantifies how much one person's biomarker profile perturbs the correlation
structure of a young, healthy reference population:

```
DISCO_i = log( n_ref² · Σ_jk w_jk (PCC_ref_jk − PCC^{ref+i}_jk)² )
```

The calculator uses a **standard reference matrix** (mean vector, cross-product
matrix, and age weights) derived from UK Biobank participants aged ≤40
(n = 5,282). The matrix contains only aggregate statistics — no individual-level
data — and is freely redistributable.

## Inputs (10 biomarkers)

CRP, glucose, total cholesterol, HDL cholesterol, creatinine, uric acid, white
blood cell count, urea, red blood cell count, albumin.

CRP is entered as raw mg/L and log-transformed (`log(CRP+1)`) automatically.

## Interpretation

- Higher DISCO = greater homeostatic dysregulation.
- Because absolute DISCO is cohort-specific (population structure + assay
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
