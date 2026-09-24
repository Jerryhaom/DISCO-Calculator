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

where $\mu$ and $\Sigma$ are the mean vector and covariance matrix of the
reference population.

## How values are estimated

Both metrics are computed entirely from a **standard reference matrix** of
aggregate statistics:

- $\mu$: reference mean vector (10 values)
- $ss$: cross-product matrix, from which the covariance $\Sigma = ss/(n_{\text{ref}}-1)$
  and correlation matrix $\rho$ are derived ($10\times10$)
- $w$: age-based weight matrix ($10\times10$)
- $n_{\text{ref}}$: reference sample size

No individual-level reference data is stored or required. Given a person's own
10 biomarker values, DISCO is obtained by the rank-one update
$ss_{\text{new}} = ss + \frac{n_{\text{ref}}}{n_{\text{ref}}+1}\,\delta\delta^{\top}$
(where $\delta = x - \mu$), and DM by the usual quadratic form
$(x-\mu)^{\top} \Sigma^{-1} (x-\mu)$. The reference matrix is derived from UK
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

## References

1. Hao M, et al. *Human aging reflects increases in entropy across organ
   networks.* (manuscript).

2. Cohen AA, et al. *A novel statistical approach shows evidence for
   multi-system physiological dysregulation during aging.* Mech Ageing Dev.
   2013;134(3-4):110–7.

## Contact

haombio@gmail.com
