# Analysis decisions and interpretation

## 1. Original missing values

The supplied CSV has 35 missing cells: 21 in CO and 14 in CO2. The original coursework runs an MCAR test and reports p = 0.856.

A p-value above 0.05 means the test does not provide sufficient evidence to reject the MCAR null hypothesis. It does not prove MCAR, and it does not establish that missingness is unrelated to unobserved values or the data-collection process. The missingness pattern and subject-matter context should also be considered.

## 2. Outlier treatment

The coursework applies Tukey's 1.5-IQR rule to each gas variable and replaces flagged observations with NA before MICE. The supplied data yield 20 additional IQR flags across the five variables. Combined with the 35 original missing cells, this gives 55 cells to impute in the coursework's cleaned dataset.

The IQR rule is a screening heuristic, not an error detector. In environmental measurements, extreme values may be meaningful. Therefore the revised script runs:
- a primary analysis that retains observed values and imputes only original missing cells;
- an IQR sensitivity analysis that additionally treats flagged observations as missing.

The two analyses should be compared rather than silently choosing the result that appears more convenient.

## 3. Multiple imputation

The original coursework uses MICE with predictive mean matching (PMM), m = 10, and takes the first completed dataset for PCA and clustering. The revised workflow retains m = 10 and PMM but repeats PCA and clustering across all completed datasets to inspect the sensitivity of cluster assignments.

The pooled linear model `CO ~ CO2 + Methane + NitrousOx + CFC11` is fitted using the primary imputation object and pooled with Rubin's rules. Coefficients are conditional associations in this dataset; they should not be described as causal effects.

## 4. PCA

PCA is run on the imputed gas variables, excluding Date, with unit-variance scaling. The number of retained components is selected as the first component count reaching at least 85% cumulative explained variance.

The coursework report contains an internal inconsistency:
- it reports PC1 = 77% and PC2 = 20.4%, which sum to approximately 97.4%;
- elsewhere it says the two retained components explain 85%.

The revised script exports the eigenvalue table and scree plot for both pipelines. The exact percentage should be taken from those generated outputs after execution, not copied from the inconsistent report text.

## 5. Clustering

The coursework selects k = 3 from an elbow-plot interpretation and reports a K-means between/total sum-of-squares ratio of 75.77403%. The revised script retains k = 3 for comparability but exports WSS and mean silhouette diagnostics for k = 2 to 8.

Ward.D2 hierarchical clustering is also cut at k = 3. Adjusted Rand index (ARI) is used to compare K-means and hierarchical assignments, and to compare K-means assignments across the ten imputed datasets. ARI measures agreement between partitions; it does not prove that clusters correspond to real atmospheric processes.

## 6. Temporal interpretation

Although the data are monthly, the coursework does not model autocorrelation, seasonality or temporal trends as part of the PCA/clustering workflow. The revised script visualises measurements over time but does not fit a time-series model. Any interpretation of clusters as pollution sources, atmospheric mechanisms or extreme events would require external validation and additional covariates.
