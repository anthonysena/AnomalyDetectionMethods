#' Summarize a value-frequency distribution
#'
#' Computes descriptive statistics for a distribution represented by values
#' with integer frequencies. The method sorts by value, uses frequency-weighted
#' calculations for the mean and median, and computes non-expanding quartiles
#' using inverse-ECDF count ranks. The mode is reported as a comma-separated
#' list of all values tied for the highest frequency.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#'
#' @return A data.frame with summary statistics: `n`, `mean`, `median`, `mode`,
#'   `q1`, `q3`, and `iqr`.
#' @export
summarizeValueFrequency <- function(df,
                                   valueColumn = "value",
                                   frequencyColumn = "frequency") {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  totalN <- sum(frequencies)

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]

  meanValue <- sum(values * frequencies) / totalN

  cumulativeFrequencies <- cumsum(frequencies)
  medianValue <- values[which(cumulativeFrequencies >= totalN / 2)[1]]

  maxFrequency <- max(frequencies)
  modes <- values[frequencies == maxFrequency]
  modeValue <- paste(modes, collapse = ", ")

  q1 <- values[which(cumulativeFrequencies >= 0.25 * totalN)[1]]
  q3 <- values[which(cumulativeFrequencies >= 0.75 * totalN)[1]]
  iqrValue <- q3 - q1

  data.frame(
    n = totalN,
    mean = meanValue,
    median = medianValue,
    mode = modeValue,
    q1 = q1,
    q3 = q3,
    iqr = iqrValue,
    stringsAsFactors = FALSE
  )
}

#' Identify outliers using Tukey fences
#'
#' Applies Tukey's method to a value-frequency distribution. Quartiles and IQR
#' are computed from the frequency-weighted distribution, then lower and upper
#' fences are calculated as `Q1 - k * IQR` and `Q3 + k * IQR`. Values outside
#' these fences are flagged as outliers. Assumptions: the method is
#' distribution-free (no normality assumption) and requires numeric values with
#' non-negative integer frequencies.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param fenceMultiplier A single positive number `k` used to scale the IQR.
#'
#' @return A data.frame with `lowerFence`, `upperFence`, and `isOutlier`
#'   columns appended.
#' @export
tukeyFences <- function(df,
                        valueColumn = "value",
                        frequencyColumn = "frequency",
                        fenceMultiplier = 1.5) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.numeric(fenceMultiplier) || length(fenceMultiplier) != 1 || fenceMultiplier <= 0) {
    stop("`fenceMultiplier` must be a single positive number.")
  }

  stats <- summarizeValueFrequency(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )

  lowerFence <- stats$q1 - fenceMultiplier * stats$iqr
  upperFence <- stats$q3 + fenceMultiplier * stats$iqr

  out <- df
  out$lowerFence <- lowerFence
  out$upperFence <- upperFence
  out$isOutlier <- (out[[valueColumn]] < lowerFence) | (out[[valueColumn]] > upperFence)

  out
}

#' Compute quantile-based outlier thresholds
#'
#' Computes lower and upper thresholds from a value-frequency distribution
#' using cumulative frequencies. Quantiles are located by the first value whose
#' cumulative frequency meets or exceeds the target rank; when `interpolate`
#' is TRUE, thresholds are linearly interpolated within the frequency bin.
#' Assumptions: the method is distribution-free (no normality assumption) and
#' requires numeric values with non-negative integer frequencies.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param lowerProb Lower quantile probability in `[0, 1)`.
#' @param upperProb Upper quantile probability in `(0, 1]`, greater than
#'   `lowerProb`.
#' @param interpolate Logical; whether to linearly interpolate within the
#'   selected frequency bin.
#'
#' @return A data.frame with `lowerThreshold`, `upperThreshold`, and
#'   `isOutlier` columns appended.
#' @export
quantileThresholds <- function(df,
                               valueColumn = "value",
                               frequencyColumn = "frequency",
                               lowerProb = 0.01,
                               upperProb = 0.99,
                               interpolate = FALSE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.numeric(lowerProb) || length(lowerProb) != 1 || lowerProb < 0 || lowerProb > 1) {
    stop("`lowerProb` must be a single number in [0, 1].")
  }
  if (!is.numeric(upperProb) || length(upperProb) != 1 || upperProb < 0 || upperProb > 1) {
    stop("`upperProb` must be a single number in [0, 1].")
  }
  if (lowerProb >= upperProb) {
    stop("`lowerProb` must be < `upperProb`.")
  }
  if (!is.logical(interpolate) || length(interpolate) != 1) {
    stop("`interpolate` must be TRUE or FALSE.")
  }

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]

  lowerThreshold <- .quantileFromFreq(values, frequencies, lowerProb, interpolate, preSorted = TRUE)
  upperThreshold <- .quantileFromFreq(values, frequencies, upperProb, interpolate, preSorted = TRUE)

  out <- df
  out$lowerThreshold <- lowerThreshold
  out$upperThreshold <- upperThreshold
  out$isOutlier <- (out[[valueColumn]] < lowerThreshold) | (out[[valueColumn]] > upperThreshold)

  out
}

#' Identify outliers using z-scores
#'
#' Computes z-scores for a value-frequency distribution using the
#' frequency-weighted mean and standard deviation. Points with absolute
#' z-scores exceeding the cutoff are flagged as outliers. Assumptions: this
#' method is best suited to roughly symmetric, bell-shaped distributions and
#' requires numeric values with non-negative integer frequencies.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param zCutoff A single positive number indicating the absolute z-score
#'   threshold used to flag outliers.
#'
#' @return A data.frame with `zScore` and `isOutlier` columns appended.
#' @export
zScoreOutliers <- function(df,
                           valueColumn = "value",
                           frequencyColumn = "frequency",
                           zCutoff = 3) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.numeric(zCutoff) || length(zCutoff) != 1 || zCutoff <= 0) {
    stop("`zCutoff` must be a single positive number.")
  }

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]

  stats <- .weightedMeanSd(values, frequencies)
  if (stats$sd <= 0) {
    stop("Standard deviation is 0; z-scores are undefined.")
  }

  zScore <- (values - stats$mean) / stats$sd

  out <- df
  out$zScore <- zScore
  out$isOutlier <- abs(zScore) > zCutoff

  out
}

#' Identify outliers using modified z-scores
#'
#' Computes modified z-scores for a value-frequency distribution using the
#' median and MAD (median absolute deviation). The constant 0.6745 rescales
#' MAD so that modified z-scores match standard z-scores under normality.
#' Points with absolute modified z-scores exceeding the cutoff are flagged as
#' outliers. Assumptions: distribution-free and robust to skew/outliers, but
#' requires numeric values with non-negative integer frequencies and a non-zero
#' MAD.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param zCutoff A single positive number indicating the absolute modified
#'   z-score threshold used to flag outliers.
#' @param interpolate Logical; whether to linearly interpolate quantiles when
#'   locating the median and MAD from cumulative frequencies.
#'
#' @return A data.frame with `modifiedZScore` and `isOutlier` columns appended.
#' @export
modifiedZScoreOutliers <- function(df,
                                   valueColumn = "value",
                                   frequencyColumn = "frequency",
                                   zCutoff = 3.5,
                                   interpolate = FALSE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.numeric(zCutoff) || length(zCutoff) != 1 || zCutoff <= 0) {
    stop("`zCutoff` must be a single positive number.")
  }
  if (!is.logical(interpolate) || length(interpolate) != 1) {
    stop("`interpolate` must be TRUE or FALSE.")
  }

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]

  medianValue <- .quantileFromFreq(values, frequencies, 0.5, interpolate, preSorted = TRUE)
  absDeviations <- abs(values - medianValue)

  orderingMad <- order(absDeviations)
  absDeviations <- absDeviations[orderingMad]
  madFrequencies <- frequencies[orderingMad]

  madValue <- .quantileFromFreq(
    absDeviations,
    madFrequencies,
    0.5,
    interpolate,
    preSorted = TRUE
  )
  if (madValue == 0) {
    stop("MAD is 0; modified z-scores are undefined.")
  }

  modifiedZScore <- 0.6745 * (values - medianValue) / madValue
  modifiedZScore <- modifiedZScore[order(ordering)]

  out <- df
  out$modifiedZScore <- modifiedZScore
  out$isOutlier <- abs(modifiedZScore) > zCutoff

  out
}

#' Identify outliers using robust Mahalanobis distance (MCD) via rrcovHD
#'
#' Fits a robust multivariate outlier detector with `rrcovHD::OutlierMahdist`
#' using MCD (`control = "mcd"` by default). Because the package methods in this
#' package accept value-frequency input, this function expands rows by
#' `frequency` to fit and score observations, then aggregates scores back to the
#' original rows.
#'
#' Assumptions: feature columns are numeric and represent a multivariate
#' observation profile for each row; frequencies are non-negative integer counts.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param featureColumns Character vector of column names to use as multivariate
#'   features. If `NULL`, all columns except `valueColumn` and `frequencyColumn`
#'   are used.
#' @param control Estimator choice/control for `rrcovHD::OutlierMahdist`.
#'   Defaults to `"mcd"`.
#' @param trace Logical; passed to `rrcovHD::OutlierMahdist`.
#' @param ... Additional arguments passed to `rrcovHD::OutlierMahdist`.
#'
#' @return A data.frame with `robustDistance`, `distanceCutoff`,
#'   `outlierProportion`, and `isOutlier` columns appended.
#' @export
mcdOutliersRrcovHD <- function(df,
                               valueColumn = "value",
                               frequencyColumn = "frequency",
                               featureColumns = NULL,
                               control = "mcd",
                               trace = FALSE,
                               ...) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!requireNamespace("rrcovHD", quietly = TRUE)) {
    stop("Package `rrcovHD` is required. Please install it.")
  }
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  if (!is.logical(trace) || length(trace) != 1) {
    stop("`trace` must be TRUE or FALSE.")
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

  frequencies <- df[[frequencyColumn]]
  positiveIdx <- which(frequencies > 0)
  if (length(positiveIdx) < 2) {
    stop("Need at least 2 rows with positive frequency for multivariate fitting.")
  }

  if (anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  expandedRowIndex <- rep(seq_len(nrow(df)), frequencies)
  expandedFeatures <- as.matrix(featuresDf[expandedRowIndex, , drop = FALSE])

  fit <- rrcovHD::OutlierMahdist(
    x = expandedFeatures,
    control = control,
    trace = trace,
    ...
  )

  distanceExpanded <- as.numeric(rrcovHD::getDistance(fit))
  cutoffExpanded <- as.numeric(rrcovHD::getCutoff(fit))
  if (length(cutoffExpanded) == 1) {
    cutoffExpanded <- rep(cutoffExpanded, length(distanceExpanded))
  }
  flagExpanded <- as.numeric(distanceExpanded > cutoffExpanded)

  robustDistance <- rep(NA_real_, nrow(df))
  distanceCutoff <- rep(NA_real_, nrow(df))
  outlierProportion <- rep(0, nrow(df))

  robustDistanceByRow <- tapply(distanceExpanded, expandedRowIndex, mean)
  distanceCutoffByRow <- tapply(cutoffExpanded, expandedRowIndex, mean)
  outlierPropByRow <- tapply(flagExpanded, expandedRowIndex, mean)

  idx <- as.integer(names(robustDistanceByRow))
  robustDistance[idx] <- as.numeric(robustDistanceByRow)
  distanceCutoff[idx] <- as.numeric(distanceCutoffByRow)
  outlierProportion[idx] <- as.numeric(outlierPropByRow)

  out <- df
  out$robustDistance <- robustDistance
  out$distanceCutoff <- distanceCutoff
  out$outlierProportion <- outlierProportion
  out$isOutlier <- outlierProportion > 0

  out
}

#' Identify outliers using studentized residuals from a linear model
#'
#' Fits a weighted linear regression model using `frequency` as case weights and
#' computes studentized residuals. Rows are flagged as outliers when the
#' absolute studentized residual exceeds `residualCutoff`.
#'
#' Assumptions: response and feature columns are numeric, frequencies are
#' non-negative integer counts used as case weights, and the linear model is
#' appropriate for the response-feature relationship.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param responseColumn Name of the response column for the regression model.
#'   Defaults to `valueColumn`.
#' @param featureColumns Character vector of predictor column names. If `NULL`,
#'   all columns except `valueColumn`, `frequencyColumn`, and `responseColumn`
#'   are used.
#' @param residualCutoff A single positive number used to flag outliers by
#'   `abs(studentizedResidual) > residualCutoff`.
#' @param deleted Logical; if TRUE uses externally studentized residuals
#'   (`stats::rstudent`), otherwise internally studentized residuals
#'   (`stats::rstandard`).
#'
#' @return A data.frame with `studentizedResidual`, `residualCutoff`,
#'   `outlierProportion`, and `isOutlier` columns appended.
#' @export
studentizedResidualOutliers <- function(df,
                                        valueColumn = "value",
                                        frequencyColumn = "frequency",
                                        responseColumn = valueColumn,
                                        featureColumns = NULL,
                                        residualCutoff = 3,
                                        deleted = TRUE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!responseColumn %in% names(df)) {
    stop(paste0("`responseColumn` not found in `df`: ", responseColumn, "."))
  }
  if (!is.numeric(df[[responseColumn]])) {
    stop("`responseColumn` must be numeric.")
  }
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  if (!is.numeric(residualCutoff) || length(residualCutoff) != 1 || residualCutoff <= 0) {
    stop("`residualCutoff` must be a single positive number.")
  }
  if (!is.logical(deleted) || length(deleted) != 1) {
    stop("`deleted` must be TRUE or FALSE.")
  }

  if (is.null(featureColumns)) {
    featureColumns <- setdiff(names(df), c(valueColumn, frequencyColumn, responseColumn))
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

  frequencies <- df[[frequencyColumn]]
  positiveIdx <- which(frequencies > 0)
  if (length(positiveIdx) < 3) {
    stop("Need at least 3 rows with positive frequency for linear-model fitting.")
  }

  if (anyNA(df[[responseColumn]][positiveIdx]) || anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("`responseColumn` and `featureColumns` cannot contain NA for rows with positive frequency.")
  }

  formula <- stats::as.formula(paste("response ~", paste(featureColumns, collapse = " + ")))
  modelDf <- df[positiveIdx, c(responseColumn, featureColumns), drop = FALSE]
  names(modelDf)[1] <- "response"
  modelWeights <- frequencies[positiveIdx]
  fit <- stats::lm(formula = formula, data = modelDf, weights = modelWeights)

  studentizedPositive <- if (deleted) {
    stats::rstudent(fit)
  } else {
    stats::rstandard(fit)
  }
  if (anyNA(studentizedPositive)) {
    stop("Studentized residuals contain NA; model may be singular or under-identified.")
  }

  studentizedResidual <- rep(NA_real_, nrow(df))
  studentizedResidual[positiveIdx] <- as.numeric(studentizedPositive)
  outlierProportion <- rep(0, nrow(df))
  outlierProportion[positiveIdx] <- as.numeric(abs(studentizedPositive) > residualCutoff)

  out <- df
  out$studentizedResidual <- studentizedResidual
  out$residualCutoff <- residualCutoff
  out$outlierProportion <- outlierProportion
  out$isOutlier <- outlierProportion > 0

  out
}

#' Identify outliers using generalized ESD
#'
#' Applies the generalized ESD (Extreme Studentized Deviate) test to a
#' value-frequency distribution, iteratively removing the most extreme point
#' and comparing the test statistic to a critical value derived from the
#' t-distribution. The number of outliers is the largest iteration where the
#' statistic exceeds the critical value. Assumptions: approximately normal
#' data, numeric values with non-negative integer frequencies.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param maxOutliers Maximum number of outliers to test for.
#' @param alpha Significance level for the ESD test.
#'
#' @return A data.frame with `outlierCount` and `isOutlier` columns appended.
#' @export
generalizedESDOutliers <- function(df,
                                   valueColumn = "value",
                                   frequencyColumn = "frequency",
                                   maxOutliers = 10,
                                   alpha = 0.05) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.numeric(maxOutliers) || length(maxOutliers) != 1 || maxOutliers < 1 || maxOutliers %% 1 != 0) {
    stop("`maxOutliers` must be a single positive integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  totalN <- sum(frequencies)
  if (totalN < 3) {
    stop("`df` must have at least 3 total observations for ESD.")
  }

  remainingFrequencies <- frequencies
  outlierCount <- rep(0, length(values))
  R <- numeric(maxOutliers)
  lambda <- numeric(maxOutliers)

  for (i in seq_len(maxOutliers)) {
    n <- sum(remainingFrequencies)
    if (n < 3) {
      R <- R[seq_len(i - 1)]
      lambda <- lambda[seq_len(i - 1)]
      break
    }

    stats <- .weightedMeanSd(values, remainingFrequencies, sample = TRUE)
    if (!is.finite(stats$sd) || stats$sd == 0) {
      R <- R[seq_len(i - 1)]
      lambda <- lambda[seq_len(i - 1)]
      break
    }

    deviations <- abs(values - stats$mean)
    idx <- which.max(deviations)
    R[i] <- deviations[idx] / stats$sd

    p <- 1 - alpha / (2 * n)
    tCrit <- stats::qt(p, df = n - 2)
    lambda[i] <- ((n - 1) * tCrit) / sqrt((n - 2 + tCrit^2) * n)

    outlierCount[idx] <- outlierCount[idx] + 1
    remainingFrequencies[idx] <- remainingFrequencies[idx] - 1
    if (remainingFrequencies[idx] < 0) {
      remainingFrequencies[idx] <- 0
    }
  }

  k <- if (length(R) == 0) 0 else max(which(R > lambda), 0)
  if (k < length(outlierCount)) {
    outlierCount <- if (k == 0) rep(0, length(outlierCount)) else pmin(outlierCount, 1)
  }

  out <- df
  out$outlierCount <- outlierCount
  out$isOutlier <- outlierCount > 0

  out
}
