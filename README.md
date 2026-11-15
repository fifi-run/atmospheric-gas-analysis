# Atmospheric Gas Data Analysis: Missing Data, PCA and Clustering

## Overview
This coursework project explores a monthly atmospheric-gas dataset with Date and five gas measurements. The submitted report describes 183 observations spanning 2000–2019. The work focuses on data quality, missing-value handling, PCA and unsupervised clustering.

This is exploratory multivariate analysis. It does **not** implement time-series forecasting or explicitly model temporal autocorrelation.

## Workflow
1. Explore distributions and pairwise correlations.
2. Flag potential outliers using Tukey's 1.5-IQR rule. The coursework workflow replaces flagged values with missing values; this is a sensitivity choice requiring scientific justification.
3. Inspect missingness and run an MCAR diagnostic. The report records p = 0.856; this means the test fails to reject MCAR, not that MCAR has been proven.
4. Apply MICE with predictive mean matching (PMM), m = 10.
5. Run PCA and select the number of components reaching at least 85% cumulative variance.
6. Compare K-means and Ward.D2 hierarchical clustering, with k = 3 retained from the coursework analysis.

## Findings reported in the coursework
- 183 observations and six columns (Date plus five gas variables).
- CO and CO2 are identified as containing missing observations.
- MCAR-test p-value reported as 0.856.
- MICE configuration: m = 10, PMM.
- Three clusters selected; reported K-means between/total sum-of-squares ratio is approximately 75.77%.

**PCA discrepancy:** the report's methods section says two components reach at least 85% cumulative variance, while its discussion states 97%. Reconcile these figures against the original PCA output before citing either as definitive.

## Code and data
- `R/atmospheric_gas_analysis.R`: cleaned, project-relative workflow. It writes diagnostic tables, PCA eigenvalues, cluster assignments and plots to `outputs/`.
- `data/README.md`: input provenance and redistribution note.

The raw CSV is not included because the supplied archive did not establish its exact upstream source or redistribution licence. Add an authorised copy at `data/MaunaLoa_miss.csv` only after confirming permission.

## Run
Install R and the packages `readr`, `dplyr`, `naniar`, `mice`, `FactoMineR`, `factoextra` and `cluster`. From the repository root run:

```r
source("R/atmospheric_gas_analysis.R")
```

The workflow has not been verified end-to-end in this update because the source CSV and an R package lockfile were not available in the repository.

## Limitations and next steps
- PCA and clustering use the first completed imputed dataset, not pooled results across imputations.
- MCAR testing does not prove a missingness mechanism.
- Extreme atmospheric observations may be scientifically meaningful; automatic outlier replacement can remove signal.
- Temporal structure is not modelled.
- Future work should compare imputation strategies, assess cluster stability across imputations and model temporal dependence explicitly.
