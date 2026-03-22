#' Identify outliers using Isolation Forest
#'
#' Fits an Isolation Forest model through `isotree::isolation.forest()` while
#' treating `frequency` as observation multiplicity through
#' `sample_weights`. With `weightsAsSampleProb = FALSE`, higher frequencies are
#' interpreted as duplicated observations rather than higher sub-sampling
#' probability. Outliers are flagged from the anomaly score returned by
#' `predict(..., type = "score")`.
#'
#' Assumptions: feature columns are numeric, frequencies are non-negative
#' integer counts interpreted as case multiplicities, and the anomaly structure
#' is well captured by recursively isolating sparse observations. Unlike LOF,
#' this method does not rely on pairwise distances and is usually less
#' sensitive to feature scaling.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param featureColumns Character vector of feature column names. If `NULL`,
#'   all columns except `valueColumn` and `frequencyColumn` are used.
#' @param ntrees A single positive integer giving the number of trees.
#' @param sampleSize Optional positive integer sub-sample size used for each
#'   tree. If `NULL`, uses the number of positive-frequency rows.
#' @param ndim A single positive integer passed to `isotree` to control the
#'   split type. `1` corresponds to standard axis-aligned Isolation Forest.
#' @param scoreCutoff Optional positive numeric cutoff for the anomaly score.
#'   If `NULL`, it is derived from `tailProb` using the frequency-weighted
#'   empirical score distribution.
#' @param tailProb Upper-tail probability in `(0, 1)` used to derive
#'   `scoreCutoff` when `scoreCutoff` is `NULL`.
#' @param weightsAsSampleProb Logical; passed to `isotree` to determine whether
#'   `frequency` acts as row-sampling probability (`TRUE`) or duplicated-row
#'   multiplicity (`FALSE`).
#' @param standardizeData Logical; passed to `isotree::isolation.forest()`.
#' @param seed Integer random seed passed to `isotree`.
#' @param nthreads Positive integer number of threads passed to `isotree`.
#'
#' @return A data.frame with `isolationForestScore`, `scoreCutoff`,
#'   `outlierProportion`, and `isOutlier` columns appended.
#' @export
isolationForestOutliers <- function(df,
                                    valueColumn = "value",
                                    frequencyColumn = "frequency",
                                    featureColumns = NULL,
                                    ntrees = 500,
                                    sampleSize = NULL,
                                    ndim = 1,
                                    scoreCutoff = NULL,
                                    tailProb = 0.99,
                                    weightsAsSampleProb = FALSE,
                                    standardizeData = FALSE,
                                    seed = 1,
                                    nthreads = 1) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!requireNamespace("isotree", quietly = TRUE)) {
    stop("Package `isotree` is required for `isolationForestOutliers()`. Install it first.")
  }
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  if (!is.numeric(ntrees) || length(ntrees) != 1 || ntrees < 1 || ntrees %% 1 != 0) {
    stop("`ntrees` must be a single positive integer.")
  }
  if (!is.null(sampleSize) && (!is.numeric(sampleSize) || length(sampleSize) != 1 || sampleSize < 1 || sampleSize %% 1 != 0)) {
    stop("`sampleSize` must be NULL or a single positive integer.")
  }
  if (!is.numeric(ndim) || length(ndim) != 1 || ndim < 1 || ndim %% 1 != 0) {
    stop("`ndim` must be a single positive integer.")
  }
  if (!is.null(scoreCutoff) && (!is.numeric(scoreCutoff) || length(scoreCutoff) != 1 || scoreCutoff <= 0)) {
    stop("`scoreCutoff` must be NULL or a single positive number.")
  }
  if (!is.numeric(tailProb) || length(tailProb) != 1 || tailProb <= 0 || tailProb >= 1) {
    stop("`tailProb` must be a single number in (0, 1).")
  }
  if (!is.logical(weightsAsSampleProb) || length(weightsAsSampleProb) != 1) {
    stop("`weightsAsSampleProb` must be TRUE or FALSE.")
  }
  if (!is.logical(standardizeData) || length(standardizeData) != 1) {
    stop("`standardizeData` must be TRUE or FALSE.")
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed %% 1 != 0) {
    stop("`seed` must be a single integer.")
  }
  if (!is.numeric(nthreads) || length(nthreads) != 1 || nthreads < 1 || nthreads %% 1 != 0) {
    stop("`nthreads` must be a single positive integer.")
  }

  if (is.null(featureColumns)) {
    featureColumns <- setdiff(names(df), c(valueColumn, frequencyColumn))
  }
  if (length(featureColumns) < 1) {
    stop("No feature columns available. Provide `featureColumns` explicitly.")
  }
  if (!all(featureColumns %in% names(df))) {
    missing <- setdiff(featureColumns, names(df))
    stop(paste0("Missing feature columns: ", paste(missing, collapse = ", "), "."))
  }

  featuresDf <- df[, featureColumns, drop = FALSE]
  if (!all(vapply(featuresDf, is.numeric, logical(1)))) {
    stop("All `featureColumns` must be numeric.")
  }

  frequencies <- as.numeric(df[[frequencyColumn]])
  positiveIdx <- which(frequencies > 0)
  if (length(positiveIdx) < 2) {
    stop("Need at least 2 rows with positive frequency for Isolation Forest.")
  }
  if (anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  x <- as.data.frame(featuresDf[positiveIdx, , drop = FALSE])
  w <- frequencies[positiveIdx]
  if (is.null(sampleSize)) {
    sampleSize <- nrow(x)
  }
  if (sampleSize > nrow(x)) {
    stop("`sampleSize` must be <= the number of rows with positive frequency.")
  }

  fit <- isotree::isolation.forest(
    data = x,
    sample_size = sampleSize,
    ntrees = as.integer(ntrees),
    ndim = as.integer(ndim),
    sample_weights = w,
    weights_as_sample_prob = weightsAsSampleProb,
    standardize_data = standardizeData,
    seed = as.integer(seed),
    nthreads = as.integer(nthreads)
  )

  scorePositive <- as.numeric(
    stats::predict(
      fit,
      newdata = x,
      type = "score",
      nthreads = as.integer(nthreads)
    )
  )

  cutoff <- if (is.null(scoreCutoff)) {
    .quantileFromFreq(scorePositive, w, tailProb, interpolate = FALSE, preSorted = FALSE)
  } else {
    scoreCutoff
  }

  isolationForestScore <- rep(NA_real_, nrow(df))
  outlierProportion <- rep(0, nrow(df))
  isolationForestScore[positiveIdx] <- scorePositive
  outlierProportion[positiveIdx] <- as.numeric(scorePositive > cutoff)

  out <- df
  out$isolationForestScore <- isolationForestScore
  out$scoreCutoff <- cutoff
  out$outlierProportion <- outlierProportion
  out$isOutlier <- outlierProportion > 0

  out
}
