#' Scale feature columns for local outlier factor analysis
#'
#' Applies frequency-aware feature scaling without expanding rows by
#' `frequency`. This is useful before running local outlier factor when feature
#' columns are measured on different scales. The available methods are:
#'
#' - `weighted_zscore`: center by the frequency-weighted mean and divide by the
#'   frequency-weighted standard deviation.
#' - `weighted_mad`: center by the frequency-weighted median and divide by the
#'   frequency-weighted MAD.
#' - `weighted_iqr`: center by the frequency-weighted median and divide by the
#'   frequency-weighted IQR.
#'
#' Assumptions: feature columns are numeric, frequencies are non-negative
#' integer counts interpreted as duplicated observations, and each selected
#' feature has non-zero spread under the chosen scaling method.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param featureColumns Character vector of feature column names. If `NULL`,
#'   all columns except `valueColumn` and `frequencyColumn` are used.
#' @param method Scaling method. One of `"weighted_zscore"`,
#'   `"weighted_mad"`, or `"weighted_iqr"`.
#' @param interpolate Logical; whether to linearly interpolate quantiles when
#'   computing weighted medians, MAD, or IQR.
#'
#' @return A data.frame with the selected `featureColumns` replaced by their
#'   scaled values.
#' @export
scaleLofFeatures <- function(df,
                             valueColumn = "value",
                             frequencyColumn = "frequency",
                             featureColumns = NULL,
                             method = c("weighted_zscore", "weighted_mad", "weighted_iqr"),
                             interpolate = FALSE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  if (!is.logical(interpolate) || length(interpolate) != 1) {
    stop("`interpolate` must be TRUE or FALSE.")
  }

  method <- match.arg(method)
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
  if (anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  out <- df
  for (featureName in featureColumns) {
    out[[featureName]] <- .scaleWeightedFeature(
      values = out[[featureName]],
      frequencies = frequencies,
      positiveIdx = positiveIdx,
      method = method,
      interpolate = interpolate
    )
  }

  out
}

#' Identify outliers using local outlier factor
#'
#' Computes the local outlier factor (LOF) for multivariate data while treating
#' `frequency` as observation multiplicity rather than expanding rows. For each
#' row with positive frequency, the method builds a weighted `k`-nearest
#' neighborhood, computes local reachability density, and then compares that
#' density to the weighted average density of its neighbors. Rows with LOF
#' scores above `scoreCutoff` are flagged as outliers. This weighted,
#' non-expanded formulation is implemented here because common R LOF packages
#' do not natively support value/frequency inputs as case multiplicities.
#'
#' Assumptions: feature columns are numeric, frequencies are non-negative
#' integer counts interpreted as duplicated observations, and nearby points are
#' expected to form locally dense neighborhoods. This method is sensitive to
#' the choice of `k` and to feature scaling.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param featureColumns Character vector of feature column names. If `NULL`,
#'   all columns except `valueColumn` and `frequencyColumn` are used.
#' @param k A single positive integer giving the neighborhood size.
#' @param scoreCutoff A single positive number used to flag outliers by
#'   `localOutlierFactor > scoreCutoff`.
#'
#' @return A data.frame with `localOutlierFactor`, `kDistance`,
#'   `localReachabilityDensity`, `scoreCutoff`, `outlierProportion`, and
#'   `isOutlier` columns appended.
#' @export
localOutlierFactorOutliers <- function(df,
                                       valueColumn = "value",
                                       frequencyColumn = "frequency",
                                       featureColumns = NULL,
                                       k = 20,
                                       scoreCutoff = 1.5) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k %% 1 != 0) {
    stop("`k` must be a single positive integer.")
  }
  if (!is.numeric(scoreCutoff) || length(scoreCutoff) != 1 || scoreCutoff <= 0) {
    stop("`scoreCutoff` must be a single positive number.")
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
  if (length(positiveIdx) < 1) {
    stop("Need at least 1 row with positive frequency for LOF.")
  }
  if (anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  x <- as.matrix(featuresDf[positiveIdx, , drop = FALSE])
  w <- frequencies[positiveIdx]
  totalWeight <- sum(w)
  if (totalWeight <= 1) {
    stop("Total positive frequency must be > 1 for LOF.")
  }
  if (k >= totalWeight) {
    stop("`k` must be smaller than the total positive frequency.")
  }

  distanceMatrix <- as.matrix(stats::dist(x))
  if (max(distanceMatrix) == 0) {
    stop("All positive-frequency rows have identical feature values; LOF is undefined.")
  }

  weightedNeighborhood <- .computeWeightedLofNeighborhood(
    distanceMatrix = distanceMatrix,
    frequencies = w,
    k = k
  )

  kDistance <- weightedNeighborhood$kDistance
  neighborhoodCounts <- weightedNeighborhood$neighborCounts
  neighborTotals <- rowSums(neighborhoodCounts)

  reachabilityDistance <- pmax(
    distanceMatrix,
    matrix(kDistance, nrow = nrow(distanceMatrix), ncol = ncol(distanceMatrix), byrow = TRUE)
  )
  reachabilitySum <- rowSums(neighborhoodCounts * reachabilityDistance)
  localReachabilityDensity <- neighborTotals / reachabilitySum

  weightedNeighborLrd <- as.numeric(neighborhoodCounts %*% localReachabilityDensity) / neighborTotals
  localOutlierFactor <- weightedNeighborLrd / localReachabilityDensity
  localOutlierFactor[is.nan(localOutlierFactor)] <- 1

  outlierProportion <- rep(0, nrow(df))
  lofScore <- rep(NA_real_, nrow(df))
  kDistanceOut <- rep(NA_real_, nrow(df))
  lrdOut <- rep(NA_real_, nrow(df))

  lofScore[positiveIdx] <- localOutlierFactor
  kDistanceOut[positiveIdx] <- kDistance
  lrdOut[positiveIdx] <- localReachabilityDensity
  outlierProportion[positiveIdx] <- as.numeric(localOutlierFactor > scoreCutoff)

  out <- df
  out$localOutlierFactor <- lofScore
  out$kDistance <- kDistanceOut
  out$localReachabilityDensity <- lrdOut
  out$scoreCutoff <- scoreCutoff
  out$outlierProportion <- outlierProportion
  out$isOutlier <- outlierProportion > 0

  out
}

.scaleWeightedFeature <- function(values,
                                  frequencies,
                                  positiveIdx,
                                  method,
                                  interpolate) {
  scaledValues <- rep(NA_real_, length(values))
  positiveValues <- values[positiveIdx]
  positiveFrequencies <- frequencies[positiveIdx]

  if (method == "weighted_zscore") {
    stats <- .weightedMeanSd(positiveValues, positiveFrequencies)
    if (!is.finite(stats$sd) || stats$sd <= 0) {
      stop("Selected feature has 0 weighted standard deviation; weighted z-score scaling is undefined.")
    }
    scaledValues[positiveIdx] <- (positiveValues - stats$mean) / stats$sd
    return(scaledValues)
  }

  ordering <- order(positiveValues)
  sortedValues <- positiveValues[ordering]
  sortedFrequencies <- positiveFrequencies[ordering]
  medianValue <- .quantileFromFreq(sortedValues, sortedFrequencies, 0.5, interpolate, preSorted = TRUE)

  if (method == "weighted_mad") {
    absDeviations <- abs(sortedValues - medianValue)
    madOrdering <- order(absDeviations)
    madValue <- .quantileFromFreq(
      absDeviations[madOrdering],
      sortedFrequencies[madOrdering],
      0.5,
      interpolate,
      preSorted = TRUE
    )
    if (!is.finite(madValue) || madValue <= 0) {
      stop("Selected feature has 0 weighted MAD; weighted MAD scaling is undefined.")
    }
    scaledValues[positiveIdx] <- (positiveValues - medianValue) / madValue
    return(scaledValues)
  }

  q1 <- .quantileFromFreq(sortedValues, sortedFrequencies, 0.25, interpolate, preSorted = TRUE)
  q3 <- .quantileFromFreq(sortedValues, sortedFrequencies, 0.75, interpolate, preSorted = TRUE)
  iqrValue <- q3 - q1
  if (!is.finite(iqrValue) || iqrValue <= 0) {
    stop("Selected feature has 0 weighted IQR; weighted IQR scaling is undefined.")
  }
  scaledValues[positiveIdx] <- (positiveValues - medianValue) / iqrValue

  scaledValues
}

.computeWeightedLofNeighborhood <- function(distanceMatrix, frequencies, k) {
  n <- nrow(distanceMatrix)
  kDistance <- numeric(n)
  neighborCounts <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n)) {
    orderedDistances <- sort(unique(distanceMatrix[i, ]))
    cumulativeCount <- 0

    for (distanceValue in orderedDistances) {
      atDistance <- which(abs(distanceMatrix[i, ] - distanceValue) <= .Machine$double.eps^0.5)
      increment <- sum(frequencies[atDistance])
      if (distanceValue == 0) {
        increment <- increment - 1
      }
      cumulativeCount <- cumulativeCount + increment
      if (cumulativeCount >= k) {
        kDistance[i] <- distanceValue
        break
      }
    }

    included <- distanceMatrix[i, ] <= (kDistance[i] + .Machine$double.eps^0.5)
    neighborCounts[i, included] <- frequencies[included]
    if (included[i]) {
      neighborCounts[i, i] <- frequencies[i] - 1
    }
  }

  list(kDistance = kDistance, neighborCounts = neighborCounts)
}
