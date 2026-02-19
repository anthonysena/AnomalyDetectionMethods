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

#' Compute a quantile from value-frequency data
#'
#' Computes the quantile for a distribution represented by values and
#' frequencies. Values are optionally sorted by value, then the quantile is
#' located by the first value whose cumulative frequency meets or exceeds the
#' target rank. When `interpolate` is TRUE, linearly interpolates within the
#' frequency bin. Assumptions: numeric values with non-negative integer
#' frequencies and a positive total count.
#'
#' @param values Numeric vector of values.
#' @param frequencies Numeric vector of frequency counts corresponding to
#'   `values`.
#' @param p Quantile probability in `[0, 1]`.
#' @param interpolate Logical; whether to linearly interpolate within the
#'   selected frequency bin.
#' @param preSorted Logical; if TRUE, assumes `values` are already sorted in
#'   ascending order and `frequencies` aligned to that ordering.
#'
#' @return A single numeric quantile value.
#'
#' @keywords internal
.quantileFromFreq <- function(values, frequencies, p, interpolate, preSorted = FALSE) {
  if (!preSorted) {
    ordering <- order(values)
    values <- values[ordering]
    frequencies <- frequencies[ordering]
  }

  totalN <- sum(frequencies)
  target <- p * totalN
  if (target <= 0) return(values[1])
  if (target >= totalN) return(values[length(values)])

  cumulativeFrequencies <- cumsum(frequencies)
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

#' Weighted mean and standard deviation
#'
#' Computes weighted moments for value-frequency data.
#'
#' @param values Numeric vector of values.
#' @param frequencies Numeric vector of frequency counts corresponding to
#'   `values`.
#' @param sample Logical; if TRUE, use sample variance (`n - 1` denominator).
#'
#' @return A list with `mean` and `sd`.
#'
#' @keywords internal
.weightedMeanSd <- function(values, frequencies, sample = FALSE) {
  totalN <- sum(frequencies)
  meanValue <- sum(values * frequencies) / totalN
  denominator <- if (sample) totalN - 1 else totalN
  variance <- sum(frequencies * (values - meanValue)^2) / denominator
  list(mean = meanValue, sd = sqrt(variance))
}
