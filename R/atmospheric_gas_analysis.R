# Atmospheric Gas Data Analysis
# Reproducible coursework workflow with an explicit outlier sensitivity analysis.
#
# Input: data/MaunaLoa_miss.csv
# Run from the repository root:
#   source("R/atmospheric_gas_analysis.R")
#
# Primary analysis: retain observed values and impute only originally missing cells.
# Sensitivity analysis: additionally flag 1.5-IQR observations as missing, matching
# the original coursework preprocessing decision. This is not automatically
# scientifically preferable: extreme atmospheric values may be real observations.

required <- c("readr", "dplyr", "tidyr", "ggplot2", "naniar",
              "mice", "FactoMineR", "factoextra", "cluster")
missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install required packages: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(naniar)
  library(mice)
  library(FactoMineR)
  library(factoextra)
  library(cluster)
})

set.seed(123)
input_file <- file.path("data", "MaunaLoa_miss.csv")
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       ". See data/README.md for data provenance requirements.")
}

raw <- readr::read_csv(input_file, show_col_types = FALSE)
if (!"Date" %in% names(raw)) stop("The input must contain a Date column.")
raw <- raw |>
  mutate(Date = as.Date(Date)) |>
  arrange(Date)

gas_names <- setdiff(names(raw), "Date")
if (length(gas_names) < 2 ||
    !all(vapply(raw[gas_names], is.numeric, logical(1)))) {
  stop("All columns other than Date must be numeric gas measurements.")
}
if (anyNA(raw$Date)) stop("Date contains missing or unparseable values.")
if (anyDuplicated(raw$Date)) stop("Date contains duplicate observations.")

gas_raw <- raw |> select(all_of(gas_names))
n_obs <- nrow(gas_raw)

# -------------------------------------------------------------------------
# 1. Data profile and exploratory analysis
# -------------------------------------------------------------------------
data_profile <- data.frame(
  variable = names(gas_raw),
  observations = n_obs,
  missing_n = colSums(is.na(gas_raw)),
  missing_pct = round(colMeans(is.na(gas_raw)) * 100, 2),
  minimum = vapply(gas_raw, min, numeric(1), na.rm = TRUE),
  median = vapply(gas_raw, median, numeric(1), na.rm = TRUE),
  maximum = vapply(gas_raw, max, numeric(1), na.rm = TRUE)
)
write.csv(data_profile, file.path(output_dir, "data_profile.csv"),
          row.names = FALSE)

long_raw <- raw |>
  pivot_longer(cols = all_of(gas_names),
               names_to = "variable", values_to = "value")

p_distribution <- ggplot(long_raw, aes(x = value)) +
  geom_histogram(bins = 30, na.rm = TRUE, fill = "steelblue") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of atmospheric gas measurements",
       x = "Observed value", y = "Count")
ggsave(file.path(output_dir, "distributions.png"), p_distribution,
       width = 10, height = 6, dpi = 300)

p_timeseries <- ggplot(long_raw, aes(x = Date, y = value)) +
  geom_line(na.rm = TRUE, linewidth = 0.4) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Monthly observations over time",
       x = "Date", y = "Observed value")
ggsave(file.path(output_dir, "monthly_series.png"), p_timeseries,
       width = 10, height = 9, dpi = 300)

correlation <- cor(gas_raw, use = "pairwise.complete.obs")
write.csv(correlation, file.path(output_dir, "pairwise_correlations.csv"))

# -------------------------------------------------------------------------
# 2. Missingness diagnostics on the original data (before outlier treatment)
# -------------------------------------------------------------------------
missing_pattern <- mice::md.pattern(gas_raw, plot = FALSE)
write.csv(missing_pattern, file.path(output_dir, "missingness_pattern.csv"))

mcar_result <- tryCatch(
  naniar::mcar_test(gas_raw),
  error = function(e) data.frame(error = conditionMessage(e))
)
write.csv(as.data.frame(mcar_result),
          file.path(output_dir, "mcar_test_original_missingness.csv"),
          row.names = FALSE)

# A non-significant MCAR test is a failure to reject MCAR, not proof that
# the data are MCAR. Interpret alongside the missingness pattern and context.

# -------------------------------------------------------------------------
# 3. Outlier diagnostics and sensitivity dataset
# -------------------------------------------------------------------------
iqr_bounds <- lapply(gas_raw, function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  spread <- IQR(x, na.rm = TRUE)
  c(lower = q[1] - 1.5 * spread,
    upper = q[2] + 1.5 * spread)
})

outlier_flags <- lapply(names(gas_raw), function(nm) {
  bounds <- iqr_bounds[[nm]]
  x <- gas_raw[[nm]]
  !is.na(x) & (x < bounds["lower"] | x > bounds["upper"])
})
names(outlier_flags) <- names(gas_raw)

outlier_summary <- data.frame(
  variable = names(outlier_flags),
  iqr_flagged_n = vapply(outlier_flags, sum, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(outlier_summary, file.path(output_dir, "iqr_outlier_summary.csv"),
          row.names = FALSE)

outlier_records <- bind_rows(lapply(names(outlier_flags), function(nm) {
  idx <- which(outlier_flags[[nm]])
  if (length(idx) == 0) return(NULL)
  data.frame(Date = raw$Date[idx], variable = nm,
             value = gas_raw[[nm]][idx])
}))
write.csv(outlier_records, file.path(output_dir, "iqr_outlier_records.csv"),
          row.names = FALSE)

gas_sensitivity <- gas_raw
for (nm in names(gas_sensitivity)) {
  gas_sensitivity[[nm]][outlier_flags[[nm]]] <- NA_real_
}

missingness_comparison <- data.frame(
  variable = names(gas_raw),
  original_missing_n = colSums(is.na(gas_raw)),
  iqr_flagged_n = vapply(outlier_flags, sum, integer(1)),
  missing_after_iqr_n = colSums(is.na(gas_sensitivity))
)
write.csv(missingness_comparison,
          file.path(output_dir, "missingness_before_after_iqr.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 4. Multiple imputation (primary and sensitivity pipelines)
# -------------------------------------------------------------------------
run_mice <- function(data, seed) {
  mice::mice(data, m = 10, method = "pmm", seed = seed,
             printFlag = FALSE, maxit = 10)
}

imp_primary <- run_mice(gas_raw, 123)
imp_sensitivity <- run_mice(gas_sensitivity, 456)

# Save imputation diagnostics and density plots.
png(file.path(output_dir, "imputation_density_primary.png"),
    width = 1400, height = 900, res = 150)
print(plot(imp_primary))
dev.off()

png(file.path(output_dir, "imputation_density_iqr_sensitivity.png"),
    width = 1400, height = 900, res = 150)
print(plot(imp_sensitivity))
dev.off()

# Pooled regression follows the coursework model specification.
pooled_model <- mice::pool(
  with(imp_primary, lm(CO ~ CO2 + Methane + NitrousOx + CFC11))
)
pooled_summary <- summary(pooled_model, conf.int = TRUE)
write.csv(pooled_summary,
          file.path(output_dir, "pooled_regression_primary.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 5. PCA and clustering helper functions
# -------------------------------------------------------------------------
select_pca <- function(completed_data, variance_threshold = 85) {
  pca <- FactoMineR::PCA(completed_data, scale.unit = TRUE, graph = FALSE)
  cumulative <- cumsum(pca$eig[, 2])
  n_components <- which(cumulative >= variance_threshold)[1]
  if (is.na(n_components)) n_components <- ncol(pca$ind$coord)
  list(model = pca, n_components = n_components,
       scores = as.data.frame(
         pca$ind$coord[, seq_len(n_components), drop = FALSE]
       ))
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  index <- sum(choose2(tab))
  row_sum <- sum(choose2(rowSums(tab)))
  col_sum <- sum(choose2(colSums(tab)))
  expected <- row_sum * col_sum / choose2(n)
  maximum <- (row_sum + col_sum) / 2
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

analyse_imputations <- function(imputed, prefix) {
  completed_sets <- lapply(seq_len(imputed$m), function(i) {
    mice::complete(imputed, action = i)
  })

  pca_results <- lapply(completed_sets, select_pca)
  first <- pca_results[[1]]
  pca <- first$model

  eigen_table <- as.data.frame(pca$eig)
  eigen_table$component <- seq_len(nrow(eigen_table))
  write.csv(eigen_table,
            file.path(output_dir, paste0(prefix, "_pca_eigenvalues.csv")),
            row.names = FALSE)

  loading_table <- as.data.frame(pca$var$coord)
  loading_table$variable <- rownames(loading_table)
  write.csv(loading_table,
            file.path(output_dir, paste0(prefix, "_pca_loadings.csv")),
            row.names = FALSE)

  p_scree <- factoextra::fviz_eig(pca, addlabels = TRUE) +
    ggtitle(paste(prefix, "PCA scree plot"))
  ggsave(file.path(output_dir, paste0(prefix, "_pca_scree.png")),
         p_scree, width = 8, height = 5, dpi = 300)

  p_biplot <- factoextra::fviz_pca_biplot(
    pca, repel = TRUE, col.var = "black", col.ind = "steelblue"
  ) + ggtitle(paste(prefix, "PCA biplot"))
  ggsave(file.path(output_dir, paste0(prefix, "_pca_biplot.png")),
         p_biplot, width = 8, height = 6, dpi = 300)

  scores <- first$scores
  k_values <- 2:min(8, nrow(scores) - 1)
  wss <- vapply(k_values, function(k) {
    set.seed(123)
    kmeans(scores, centers = k, nstart = 25)$tot.withinss
  }, numeric(1))
  silhouette <- vapply(k_values, function(k) {
    set.seed(123)
    fit <- kmeans(scores, centers = k, nstart = 25)
    mean(cluster::silhouette(fit$cluster, dist(scores))[, "sil_width"])
  }, numeric(1))
  write.csv(data.frame(k = k_values, total_within_ss = wss,
                       mean_silhouette = silhouette),
            file.path(output_dir, paste0(prefix, "_cluster_diagnostics.csv")),
            row.names = FALSE)

  # k = 3 is retained as the original coursework choice; diagnostics are
  # provided so readers can inspect rather than assume it is optimal.
  k <- 3L
  set.seed(123)
  km <- kmeans(scores, centers = k, nstart = 25)
  km_sil <- cluster::silhouette(km$cluster, dist(scores))
  km_summary <- data.frame(
    k = k,
    between_total_ss_pct = 100 * km$betweenss / km$totss,
    mean_silhouette = mean(km_sil[, "sil_width"]),
    n_components = first$n_components,
    cumulative_variance_pct = sum(pca$eig[seq_len(first$n_components), 2])
  )
  write.csv(km_summary,
            file.path(output_dir, paste0(prefix, "_kmeans_summary.csv")),
            row.names = FALSE)
  write.csv(data.frame(observation = seq_len(nrow(scores)),
                       cluster = km$cluster),
            file.path(output_dir, paste0(prefix, "_kmeans_assignments.csv")),
            row.names = FALSE)

  hc <- hclust(dist(scores), method = "ward.D2")
  hc_cluster <- cutree(hc, k = k)
  write.csv(data.frame(observation = seq_len(nrow(scores)),
                       cluster = hc_cluster),
            file.path(output_dir, paste0(prefix, "_hierarchical_assignments.csv")),
            row.names = FALSE)
  write.csv(data.frame(adjusted_rand_index =
                         adjusted_rand_index(km$cluster, hc_cluster)),
            file.path(output_dir, paste0(prefix, "_kmeans_vs_hierarchical_ari.csv")),
            row.names = FALSE)

  p_km <- factoextra::fviz_cluster(
    list(data = scores, cluster = km$cluster)
  ) + ggtitle(paste(prefix, "K-means clustering (k = 3)"))
  ggsave(file.path(output_dir, paste0(prefix, "_kmeans_clusters.png")),
         p_km, width = 8, height = 6, dpi = 300)

  png(file.path(output_dir, paste0(prefix, "_hierarchical_dendrogram.png")),
      width = 1400, height = 900, res = 150)
  plot(hc, labels = FALSE, main = paste(prefix, "Ward.D2 dendrogram"))
  rect.hclust(hc, k = k, border = 2:4)
  dev.off()

  # Compare k-means assignments across all ten imputed datasets. ARI is
  # invariant to arbitrary cluster label numbering.
  all_km <- lapply(pca_results, function(res) {
    set.seed(123)
    kmeans(res$scores, centers = k, nstart = 25)$cluster
  })
  stability <- data.frame(
    imputation = seq_along(all_km),
    ari_vs_imputation_1 = vapply(all_km, function(z) {
      adjusted_rand_index(all_km[[1]], z)
    }, numeric(1))
  )
  write.csv(stability,
            file.path(output_dir, paste0(prefix, "_imputation_cluster_stability.csv")),
            row.names = FALSE)

  list(pca = pca, kmeans = km, hierarchical = hc,
       n_components = first$n_components)
}

primary_results <- analyse_imputations(imp_primary, "primary")
sensitivity_results <- analyse_imputations(imp_sensitivity, "iqr_sensitivity")

# Compare the primary and IQR-sensitivity clustering solutions.
primary_assignments <- read.csv(
  file.path(output_dir, "primary_kmeans_assignments.csv")
)$cluster
sensitivity_assignments <- read.csv(
  file.path(output_dir, "iqr_sensitivity_kmeans_assignments.csv")
)$cluster
write.csv(
  data.frame(adjusted_rand_index_primary_vs_iqr =
               adjusted_rand_index(primary_assignments, sensitivity_assignments)),
  file.path(output_dir, "primary_vs_iqr_sensitivity_ari.csv"),
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))

message("Analysis completed. Review outputs/ and compare primary results with the IQR sensitivity analysis.")
