# Atmospheric Gas Data Analysis: Missing Data, PCA and Clustering

An exploratory multivariate analysis of monthly atmospheric gas measurements from January 2000 to March 2015? No: the supplied source data spans January 2000 to March 2015? 

## Project overview

This project was completed as a data visualisation and analysis coursework assignment. It examines five measured variables (CO, CO2, Methane, NitrousOx and CFC11) across 183 monthly observations, using missing-data diagnostics, multiple imputation, principal component analysis (PCA), K-means and hierarchical clustering.

The work is **exploratory**. It does not estimate causal effects, identify verified pollution sources, or model temporal autocorrelation / forecast future concentrations.

## Research questions

1. What are the distributions, pairwise relationships and missingness patterns in the five measurements?
2. How do imputation and the treatment of extreme observations affect the multivariate representation?
3. How many principal components are needed to retain at least 85% of variance?
4. Do K-means and Ward.D2 hierarchical clustering produce comparable groupings in PCA space?

## Dataset

The source file supplied with the coursework is named `MaunaLoa_miss.csv`. It contains 183 monthly records and six columns: Date plus five numeric gas measurements. The date range is January 2000 to March 2015? 

Original missingness in the supplied CSV:
| Variable | Missing observations |
|---|---:|
| CO | 21 |
| CO2 | 14 |
| Methane | 0 |
| NitrousOx | 0 |
| CFC11 | 0 |
| **Total** | **35** |

The source archive did not contain a reliable upstream URL, version or redistribution licence. The raw CSV is therefore not committed. See [data/README.md](data/README.md) before obtaining or redistributing the data.

## Analysis workflow

The revised R workflow deliberately separates two pipelines:

- **Primary pipeline:** preserves all observed values and imputes only the 35 originally missing cells using MICE with predictive mean matching (PMM), 10 imputations.
- **IQR sensitivity pipeline:** additionally flags observations outside Tukey's 1.5-IQR bounds and treats them as missing before MICE. This reproduces the original coursework's outlier-handling decision, but does not assume that extreme atmospheric observations are erroneous.

Both pipelines produce PCA, clustering diagnostics, K-means assignments, Ward.D2 hierarchical assignments and cluster-stability comparisons across imputations. The primary workflow also runs the pooled regression specified in the coursework (`CO ~ CO2 + Methane + NitrousOx + CFC11`).

## Key findings in the original coursework report

- The report records an MCAR-test p-value of 0.856. This is a failure to reject the MCAR null hypothesis, not proof that missingness is completely at random.
- The coursework applies IQR outlier replacement before imputation, resulting in 55 missing cells (35 originally missing plus 20 IQR-flagged observations).
- The report selects three clusters and reports a K-means between/total sum-of-squares ratio of approximately 75.77%.
- The report's PCA variance statements are inconsistent: one section reports PC1 at 77% and PC2 at 20.4% (about 97.4% combined), while another says the two retained components explain 85%. The revised script recalculates and exports the eigenvalues rather than choosing one conflicting figure.

These are results reported in the original coursework, not outputs from a verified execution of the revised script.

## Repository structure

```text
.
├── R/
│   └── atmospheric_gas_analysis.R
├── data/
│   └── README.md
├── docs/
│   ├── ANALYSIS_DECISIONS.md
│   └── DATA_DICTIONARY.md
├── outputs/
│   └── .gitkeep
├── .gitignore
└── README.md
```

## Reproduce the analysis

1. Install R (4.2 or later recommended).
2. Obtain the dataset from its original authorised source and save it as `data/MaunaLoa_miss.csv`.
3. Install the required packages:

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "ggplot2", "naniar",
  "mice", "FactoMineR", "factoextra", "cluster"
))
```

4. From the repository root, run:

```r
source("R/atmospheric_gas_analysis.R")
```

The script writes tables, plots, cluster assignments and `sessionInfo()` to `outputs/`.

**Execution status:** the revised script has been committed but has not yet been run in an R environment. Do not cite newly generated metrics until the script has been executed with the authorised input file and its outputs have been inspected.

## Methodological limitations

- MCAR testing does not prove the missingness mechanism.
- IQR is a univariate screening rule; it does not establish that an observation is an error.
- PCA and clustering are sensitive to scaling and imputation choices.
- Cluster labels are descriptive groupings, not validated atmospheric mechanisms or pollution-source categories.
- The analysis does not account for serial dependence, seasonality or external drivers.
- Cluster stability across imputations is assessed with adjusted Rand index, but this does not replace external validation or domain interpretation.

See [Analysis decisions](docs/ANALYSIS_DECISIONS.md) and the [data dictionary](docs/DATA_DICTIONARY.md) for more detail.
