#' Validate a value-frequency data frame
#'
#' Performs shared validation for data frames that contain a value column and a
#' frequency column used by the cross-sectional methods. Checks that the input
#' is a data.frame, required columns exist, values are numeric without missing
#' entries, frequencies are non-negative integers, and the total frequency is
#' non-zero.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#'
#' @keywords internal
.validateValueFrequencyDf <- function(df,
  valueColumn = "value",
  frequencyColumn = "frequency") {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }

  requiredColumns <- c(valueColumn, frequencyColumn)
  if (!all(requiredColumns %in% names(df))) {
    stop(
      paste0(
        "`df` must contain columns: ",
        paste(requiredColumns, collapse = ", "),
        "."
      )
    )
  }
  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]

  if (length(values) == 0) {
    stop("`df` must have at least one row.")
  }
  if (anyNA(values) || anyNA(frequencies)) {
    stop("`value` and `frequency` columns cannot contain NA.")
  }
  if (!is.numeric(values)) {
    stop("`value` column must be numeric.")
  }
  if (!is.numeric(frequencies)) {
    stop("`frequency` column must be numeric.")
  }
  if (any(frequencies < 0)) {
    stop("`frequency` values must be >= 0.")
  }
  if (any(frequencies %% 1 != 0)) {
    stop("`frequency` values must be integers (counts).")
  }

  totalN <- sum(frequencies)
  if (totalN == 0) {
    stop("Total frequency is 0; nothing to summarize.")
  }
}
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

  # Non-expanding quartiles (type = "inverse-ECDF" with count ranks):
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

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]

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
  totalN <- sum(frequencies)

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]

  cumulativeFrequencies <- cumsum(frequencies)

  quantileFromFreq <- function(p) {
    target <- p * totalN
    if (target <= 0) return(values[1])
    if (target >= totalN) return(values[length(values)])
    idx <- which(cumulativeFrequencies >= target)[1]
    if (!interpolate) {
      return(values[idx])
    }
    prevCount <- if (idx > 1) cumulativeFrequencies[idx - 1] else 0
    if (frequencies[idx] == 0 || target == prevCount) return(values[idx])
    prop <- (target - prevCount) / frequencies[idx]
    prevValue <- if (idx > 1) values[idx - 1] else values[1]
    currValue <- values[idx]
    prevValue + prop * (currValue - prevValue)
  }

  lowerThreshold <- quantileFromFreq(lowerProb)
  upperThreshold <- quantileFromFreq(upperProb)

  out <- df
  out$lowerThreshold <- lowerThreshold
  out$upperThreshold <- upperThreshold
  out$isOutlier <- (out[[valueColumn]] < lowerThreshold) | (out[[valueColumn]] > upperThreshold)

  out
}
