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
correlation structure of a young reference population:

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
\mathrm{DM}_i = \sqrt{ (x_i - \mu)^{\top} \Sigma^{-1} (x_i - \mu) }
$$

Here $x_i$ is the vector after CRP transformation and reference standardization,
and $\mu$ is the mean vector on that standardized scale. The implemented DM
matrix is $\Sigma = ss/(n_{\text{ref}}-1) + 10^{-8}I$.
DISCO and the CRP preprocessing step use natural logarithms. **DM is reported
as the untransformed Mahalanobis distance**, without a logarithm.

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
(where $\delta = x - \mu$), and DM by the square root of the quadratic form
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
- The main DISCO and DM results always use all 10 markers and the same UK
  Biobank reference matrix. Selecting a cohort changes the relative comparison.
  CHARLS's stored distribution uses 8 markers (excluding RBC/ALB); RuLAS's
  uses 9 (excluding UREA). Their comparison DISCO is therefore recalculated
  from the corresponding subset and labeled with its panel size. The other
  cohorts use the full 10-marker score. Subset matrices retain the original
  weights without renormalization, matching the reference-distribution source.
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

`STANDARD_REF` and `COHORT_STATS` are bundled in `index.html`. The local source
analyses are `19_rebuild_logcrp_plus1.R` (reference) and
`21_disco_zscore_percentiles.R` (cohort distributions). The latter's exported
means, SDs and quantile knots match the bundled values, and confirm the 10/8/9
panel sizes. The distributions use all participants with a finite score,
whereas the mortality analyses below additionally require complete covariates
and valid survival outcomes. These are different analysis samples.

The original parameter-export scripts and participant data are not included
in this repository. The new aggregate evidence and its reconstruction script
make the mortality estimates auditable, but do not establish clinical cutoffs,
absolute risk calibration, or held-out predictive performance.

## Biomarker coverage and mortality evidence

The page includes an accessible 10-by-5 availability heatmap and separate
DISCO/DM mortality results. UKB, NHANES and CLHLS use 10 markers; CHARLS uses
8; RuLAS wave 2 uses 9. Coverage is specific to the analyzed waves, not every
wave of each cohort. The accompanying aggregate results are in
[`data/cohort-evidence.json`](data/cohort-evidence.json).

- Separate Cox models use the same complete-case participants for DISCO and
  untransformed DM within each cohort. HRs are per 1 SD in that analysis sample.
- All models include age, sex, BMI, education, smoking and alcohol use; UKB and
  NHANES additionally include ethnicity/race, and UKB includes physical activity.
  Age and BMI are continuous; other covariates are categorical.
- C-index means Harrell's C for the **full score + covariates model**, fitted
  and evaluated in the same cohort. It is apparent performance, without
  optimism correction. These are not stand-alone biomarker C-indices or
  held-out validation of one fixed mortality model.
- Models use Efron ties and 95% Wald confidence intervals. NHANES is unweighted;
  survey-design inference is not implemented. No additional biomarker
  imputation is performed on the supplied cleaned inputs.
- Score-specific proportional-hazards tests detect departures for UKB DISCO
  and DM, NHANES DISCO and CHARLS DISCO. Their fitted HRs should not be read as
  time-constant effects; time-varying sensitivity analyses remain necessary.
  Test p-values are included in the JSON.

To regenerate with authorized access to the original controlled data, install
R packages `survival`, `jsonlite`, and `openxlsx`, then run:

```sh
export DISCO_DATA_ROOT=/path/to/controlled/DISCO
export DISCO_RULAS_QUESTIONNAIRE=/path/to/controlled/wave2.xlsx
Rscript scripts/rebuild-cohort-evidence.R
python3 scripts/render-cohort-evidence.py
```

The R script checks join-key uniqueness, reconstructs unit conversions,
excludes invalid outcomes (including CLHLS lost-to-follow-up codes), and checks
vectorized DISCO against literal matrix updates. Detailed variable mappings
and covariate recodes are in that script. Only aggregate JSON is written;
the controlled source files must remain outside this repository. The Python
renderer embeds the tables in HTML so they also work when opened as a local file.

## Checks

Run the dependency-free regression checks with Node.js 18 or newer:

```sh
node --test tests/*.test.cjs
python3 scripts/render-cohort-evidence.py --check
```

The checks cover fixed score examples, percentile knots and tails for all five
cohorts, interpolation and ordinal formatting, invalid input and overflow
handling, result invalidation after edits, panel-matched cohort comparisons,
and aggregate evidence integrity. They verify implementation behavior, not
clinical validity.

## References

1. Hao M, et al. *Human aging reflects increases in entropy across organ
   networks.* (manuscript).

2. Cohen AA, et al. *A novel statistical approach shows evidence for
   multi-system physiological dysregulation during aging.* Mech Ageing Dev.
   2013;134(3-4):110–7.

## Contact

haombio@gmail.com
