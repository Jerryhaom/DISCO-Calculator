# DISCO Calculator

Web calculator for **DISCO (Distance of Covariance)**, a measure of entropy in
large ensembles of biological information. DISCO demonstrates that organs and
systems exhibit interconnected increased entropy with age, and it predicts
mortality, frailty, and age-related chronic disease across multiple data types
(clinical biomarkers, proteomics, metabolomics, microbiomes) in five cohorts
(UK Biobank, NHANES, and three Chinese cohorts of older adults). The companion
**DM (Mahalanobis distance)** is also provided.

Live: https://jerryhaom.github.io/DISCO-Calculator/

## What it computes

DISCO quantifies how much one person's biomarker profile perturbs the
correlation structure of a young, healthy reference population:

$$
\mathrm{DISCO}_i = \log\Bigg[ n_{\mathrm{ref}}^{2} \sum_{j \neq k} w_{jk} \Big( \rho_{jk}^{\mathrm{ref}} - \rho_{jk}^{\mathrm{ref}+i} \Big)^{2} \Bigg]
$$

where $\rho^{\mathrm{ref}}$ is the reference correlation matrix,
$\rho^{\mathrm{ref}+i}$ is the correlation matrix after adding individual
$i$, $w$ is an age-based weight matrix, and $n_{\mathrm{ref}}$ is the reference
sample size.

DM (Mahalanobis distance) quantifies the standardized deviation of the
biomarker vector from the reference mean:

$$
\mathrm{DM}_i = \log\Bigg[ \sqrt{ (x_i - \mu)^{\top} \Sigma^{-1} (x_i - \mu) } \Bigg]
$$

Here $x_i$ is the vector after CRP transformation and reference standardization,
and $\mu$ is the mean vector on that standardized scale. The implemented DM
matrix is $\Sigma = ss/(n_{\text{ref}}-1) + 10^{-8}I$.
All logarithms are natural logarithms. The main DM result is **log DM**;
the untransformed distance is also displayed as **raw DM**.

## How values are estimated

Both metrics are computed entirely from a **standard reference matrix** of
aggregate statistics:

- `ref_mean` and `ref_sd`: reference means and standard deviations used to
  standardize the 10 markers after transforming CRP
- $\mu$: reference mean vector on the standardized scale (10 values)
- $ss$: centered cross-product matrix on the standardized scale, from which
  the covariance and correlation matrices are derived ($10\times10$)
- $w$: age-based weight matrix ($10\times10$)
- $n_{\text{ref}}$: reference sample size

No individual-level reference data is stored or required. Given a person's own
10 biomarker values, DISCO is obtained by the rank-one update
$ss_{\text{new}} = ss + \frac{n_{\text{ref}}}{n_{\text{ref}}+1}\,\delta\delta^{\top}$
(where $\delta = x - \mu$), and DM by the usual quadratic form
$(x-\mu)^{\top} \Sigma^{-1} (x-\mu)$. The reference matrix is derived from UK
Biobank participants aged $\le 40$ (n = 5,282) and is freely redistributable.

Given these aggregate statistics, the rank-one update is algebraically exact
up to parameter rounding and floating-point precision; individual reference
records are not required. The full off-diagonal DISCO sum uses symmetric
weights normalized to sum to one over the full matrix.

## Inputs (10 biomarkers)

CRP, glucose, total cholesterol, HDL cholesterol, creatinine, uric acid, white
blood cell count, urea, red blood cell count, albumin.

CRP is entered as raw mg/L and log-transformed ($\log(\text{CRP}+1)$)
automatically.

All 10 values must be finite and non-negative, in the displayed units. Missing
values are not imputed. The displayed normal ranges are informational and do
not enter the score. Invalid inputs and non-finite calculation results produce
an error instead of a result or interpretation. Editing an input or changing
the cohort hides the previous results until the calculator is run again.

## Interpretation

- Higher DISCO / DM = greater homeostatic dysregulation.
- DISCO and DM always use the same UK Biobank reference matrix. Selecting
  UK Biobank, NHANES, CHARLS, CLHLS, or RuLAS changes only the DISCO z-score
  and estimated percentile.
- The DISCO z-score is `(DISCO - cohort mean) / cohort standard deviation`.
  It describes standing relative to the selected cohort, without adjusting
  for the individual's age or sex.
- Percentiles are linearly interpolated between stored z-score thresholds at
  the 5th, 10th, 25th, 50th, 75th, 90th, and 95th percentiles. Within that
  range, rounded estimates are marked with `≈`. Below or above the available
  thresholds, the result is `<5th` or `>95th`; an exact tail rank is unknown.
  Values at the endpoints correspond to the 5th and 95th percentiles. A
  tolerance of `1e-12` in z-score units preserves knots against floating-point
  error when reconstructing z-scores from scores.
- The result label compares DISCO with the cohort mean. It does not assign
  clinical risk categories or estimate disease probability.

## Parameter provenance and validation

`STANDARD_REF` and `COHORT_STATS` are bundled in `index.html`. This repository
does not currently include the underlying data or the parameter export scripts.
The cohort statistics must correspond to exactly the same 10 markers,
transformations, reference matrix, weight normalization, and log-score scale
used here. Their provenance and external calibration cannot be established
from this repository alone. A parameter export script and versioned validation
results are needed to substantiate predictive-performance claims or clinical
cutoffs; no such cutoffs are used by the interface.

## Checks

Run the dependency-free regression checks with Node.js 18 or newer:

```sh
node --test tests/calculator.test.cjs
```

The checks cover fixed score examples, percentile knots and tails for all five
cohorts, interpolation and ordinal formatting, invalid input and overflow
handling, and result invalidation after edits. They verify implementation
behavior, not clinical validity or the source of the bundled parameters.

## References

1. Hao M, et al. *Human aging reflects increases in entropy across organ
   networks.* (manuscript).

2. Cohen AA, et al. *A novel statistical approach shows evidence for
   multi-system physiological dysregulation during aging.* Mech Ageing Dev.
   2013;134(3-4):110–7.

## Contact

haombio@gmail.com
