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
