# DISCO Calculator

Web calculator for **DISCO (Distance of Covariance)** and **DM (Mahalanobis
distance)**, measures of homeostatic dysregulation that reveal organ-system
interconnections underlying mortality and disease risk.

Live: https://jerryhaom.github.io/DISCO-Calculator/

## What it computes

DISCO quantifies how much one person's biomarker profile perturbs the
correlation structure of a young, healthy reference population:

$$
\text{DISCO}_i = \log\!\Big( n_{\text{ref}}^2 \sum_{j,k} w_{jk}\, \big( \text{PCC}_{jk}^{\text{ref}} - \text{PCC}_{jk}^{\text{ref}+i} \big)^2 \Big)
$$

where $\text{PCC}^{\text{ref}}$ is the reference correlation matrix,
$\text{PCC}^{\text{ref}+i}$ is the correlation matrix after adding individual
$i$, $w$ is an age-based weight matrix, and $n_{\text{ref}}$ is the reference
sample size.

DM (Mahalanobis distance) quantifies the standardized deviation of the
biomarker vector from the reference mean:

$$
\text{DM}_i = \log\!\Big( \sqrt{ (x_i - \mu)^\top \, S^{-1} \, (x_i - \mu) } \Big)
$$

where $\mu$ and $S$ are the mean vector and covariance matrix of the reference
population.

## How values are estimated — no raw data required

Both metrics are computed entirely from a **standard reference matrix** of
aggregate statistics:

- $\mu$ — reference mean vector (10 values)
- $ss$ — cross-product matrix, from which the covariance $S = ss/(n_{\text{ref}}-1)$
  and correlation matrix $\text{PCC}$ are derived ($10\times10$)
- $w$ — age-based weight matrix ($10\times10$)
- $n_{\text{ref}}$ — reference sample size

No individual-level reference data is stored or required. Given a person's own
10 biomarker values, DISCO is obtained by the rank-one update
$ss_{\text{new}} = ss + \frac{n_{\text{ref}}}{n_{\text{ref}}+1}\,\delta\delta^{\top}$
(where $\delta = x - \mu$), and DM by the usual quadratic form
$(x-\mu)^{\top} S^{-1} (x-\mu)$. The reference matrix is derived from UK
Biobank participants aged $\le 40$ (n = 5,282) and is freely redistributable.

## Inputs (10 biomarkers)

CRP, glucose, total cholesterol, HDL cholesterol, creatinine, uric acid, white
blood cell count, urea, red blood cell count, albumin.

CRP is entered as raw mg/L and log-transformed ($\log(\text{CRP}+1)$)
automatically.

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
