#' Compute simple diagnostics for a value-frequency stratum
#'
#' Summarizes stratum size, concentration, spread, ties, heaping, and method
#' readiness indicators for a value-frequency data frame.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the numeric value column.
#' @param frequencyColumn Name of the numeric frequency column.
#' @param kValues Optional positive integer vector used to evaluate LOF
#'   neighborhood readiness.
#' @param featureColumns Optional character vector of numeric feature columns
#'   used to evaluate Gaussian covariance readiness.
#'
#' @return A one-row data.frame of stratum diagnostics.
#' @export
summarizeStratumDiagnostics <- function(df,
                                        valueColumn = "value",
                                        frequencyColumn = "frequency",
                                        kValues = c(5, 10, 20),
                                        featureColumns = NULL) {
  .validateValueFrequencyDf(df, valueColumn = valueColumn, frequencyColumn = frequencyColumn)

  parts <- c(
    .computeStratumSizeDiagnostics(df, valueColumn, frequencyColumn),
    .computeStratumConcentrationDiagnostics(df, valueColumn, frequencyColumn),
    .computeStratumSpreadDiagnostics(df, valueColumn, frequencyColumn),
    .computeStratumTieDiagnostics(df, valueColumn, frequencyColumn),
    .computeStratumHeapingDiagnostics(df, valueColumn, frequencyColumn),
    .computeStratumNeighborReadiness(df, valueColumn, frequencyColumn, kValues),
    .computeStratumCovarianceReadiness(df, featureColumns, frequencyColumn)
  )

  as.data.frame(parts, stringsAsFactors = FALSE)
}

#' Summarize diagnostics across multiple strata
#'
#' Applies [summarizeStratumDiagnostics()] within each combination of the
#' supplied grouping columns.
#'
#' @param df A data.frame containing stratum identifiers plus value-frequency
#'   columns.
#' @param strataColumns Character vector of grouping columns that define strata.
#' @param valueColumn Name of the numeric value column.
#' @param frequencyColumn Name of the numeric frequency column.
#' @param kValues Optional positive integer vector used to evaluate LOF
#'   neighborhood readiness.
#' @param featureColumns Optional character vector of numeric feature columns
#'   used to evaluate Gaussian covariance readiness.
#'
#' @return A data.frame with one row per stratum.
#' @export
summarizeStrataDiagnostics <- function(df,
                                       strataColumns,
                                       valueColumn = "value",
                                       frequencyColumn = "frequency",
                                       kValues = c(5, 10, 20),
                                       featureColumns = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  if (!is.character(strataColumns) || length(strataColumns) < 1) {
    stop("`strataColumns` must be a non-empty character vector.")
  }
  if (!all(c(strataColumns, valueColumn, frequencyColumn) %in% names(df))) {
    stop("`df` is missing required stratum or value-frequency columns.")
  }

  splitKeys <- interaction(df[, strataColumns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, splitKeys)
  out <- lapply(pieces, function(piece) {
    diagRow <- summarizeStratumDiagnostics(
      df = piece,
      valueColumn = valueColumn,
      frequencyColumn = frequencyColumn,
      kValues = kValues,
      featureColumns = featureColumns
    )
    cbind(piece[1, strataColumns, drop = FALSE], diagRow, stringsAsFactors = FALSE)
  })

  do.call(rbind, out)
}

#' Flag strata that appear analyzable for LOF, IF, and GAD
#'
#' Applies simple readiness heuristics to the output of
#' [summarizeStratumDiagnostics()] or [summarizeStrataDiagnostics()].
#'
#' @param diagnosticsDf A data.frame of stratum diagnostics.
#' @param minDistinctValues Minimum number of distinct values required.
#' @param minTotalFrequency Minimum total frequency required.
#' @param lofK Target LOF neighborhood size for readiness checks.
#' @param maxTop1Mass Maximum allowed mass in the single most common value.
#' @param maxTop5Mass Maximum allowed mass in the five most common values.
#'
#' @return A data.frame with `usableForLof`, `usableForIf`, and `usableForGad`
#'   logical columns appended.
#' @export
flagAnalyzableStrata <- function(diagnosticsDf,
                                 minDistinctValues = 20,
                                 minTotalFrequency = 100,
                                 lofK = 10,
                                 maxTop1Mass = 0.5,
                                 maxTop5Mass = 1) {
  if (!is.data.frame(diagnosticsDf)) {
    stop("`diagnosticsDf` must be a data.frame.")
  }
  required <- c("distinctValues", "totalFrequency", "top1Mass", "top5Mass")
  if (!all(required %in% names(diagnosticsDf))) {
    stop("`diagnosticsDf` is missing required diagnostic columns.")
  }

  supportColumn <- paste0("supportsK_", lofK)
  if (!supportColumn %in% names(diagnosticsDf)) {
    diagnosticsDf[[supportColumn]] <- diagnosticsDf$distinctValues > lofK
  }

  out <- diagnosticsDf
  out$usableForLof <- (
    out$distinctValues >= max(minDistinctValues, lofK + 1) &
      out$totalFrequency >= minTotalFrequency &
      out[[supportColumn]] &
      out$top1Mass <= maxTop1Mass &
      out$top5Mass <= maxTop5Mass
  )
  out$usableForIf <- (
    out$distinctValues >= max(10, minDistinctValues %/% 2) &
      out$totalFrequency >= max(50, minTotalFrequency %/% 2) &
      out$top1Mass <= maxTop1Mass
  )
  gadVarianceOk <- is.na(out$hasZeroVarianceFeature) | !out$hasZeroVarianceFeature
  out$usableForGad <- (
    out$distinctValues >= max(10, minDistinctValues %/% 2) &
      out$totalFrequency >= minTotalFrequency &
      gadVarianceOk
  )

  out
}

#' Plot selected stratum diagnostics
#'
#' Creates a scatter plot of two diagnostic metrics, with optional color and
#' point size mappings.
#'
#' @param diagnosticsDf A data.frame of stratum diagnostics.
#' @param xColumn Name of the x-axis diagnostic column.
#' @param yColumn Name of the y-axis diagnostic column.
#' @param colorColumn Optional name of a column used for color.
#' @param sizeColumn Optional name of a numeric column used for point size.
#'
#' @return A ggplot object.
#' @export
plotStratumDiagnostics <- function(diagnosticsDf,
                                   xColumn,
                                   yColumn,
                                   colorColumn = NULL,
                                   sizeColumn = NULL) {
  if (!is.data.frame(diagnosticsDf)) {
    stop("`diagnosticsDf` must be a data.frame.")
  }
  required <- c(xColumn, yColumn)
  if (!is.null(colorColumn)) {
    required <- c(required, colorColumn)
  }
  if (!is.null(sizeColumn)) {
    required <- c(required, sizeColumn)
  }
  if (!all(required %in% names(diagnosticsDf))) {
    stop("`diagnosticsDf` is missing required plotting columns.")
  }

  aesArgs <- list(x = rlang::sym(xColumn), y = rlang::sym(yColumn))
  if (!is.null(colorColumn)) {
    aesArgs$color <- rlang::sym(colorColumn)
  }
  if (!is.null(sizeColumn)) {
    aesArgs$size <- rlang::sym(sizeColumn)
  }

  ggplot2::ggplot(diagnosticsDf, do.call(ggplot2::aes, aesArgs)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::labs(x = xColumn, y = yColumn, color = colorColumn, size = sizeColumn)
}

.computeStratumSizeDiagnostics <- function(df, valueColumn, frequencyColumn) {
  frequencies <- df[[frequencyColumn]]
  list(
    distinctValues = nrow(df),
    totalFrequency = sum(frequencies),
    meanFrequencyPerValue = mean(frequencies),
    medianFrequencyPerValue = stats::median(frequencies),
    maxFrequency = max(frequencies),
    singleCountValues = sum(frequencies == 1)
  )
}

.computeStratumConcentrationDiagnostics <- function(df, valueColumn, frequencyColumn) {
  frequencies <- sort(df[[frequencyColumn]], decreasing = TRUE)
  total <- sum(frequencies)
  probs <- frequencies / total
  entropy <- -sum(probs * log(probs))
  list(
    top1Mass = frequencies[1] / total,
    top5Mass = sum(utils::head(frequencies, 5)) / total,
    top10Mass = sum(utils::head(frequencies, 10)) / total,
    entropy = entropy,
    effectiveValueCount = exp(entropy)
  )
}

.computeStratumSpreadDiagnostics <- function(df, valueColumn, frequencyColumn) {
  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  stats <- .weightedMeanSd(values, frequencies)
  ordering <- order(values)
  sortedValues <- values[ordering]
  sortedFreq <- frequencies[ordering]
  medianValue <- .quantileFromFreq(sortedValues, sortedFreq, 0.5, FALSE, preSorted = TRUE)
  q1 <- .quantileFromFreq(sortedValues, sortedFreq, 0.25, FALSE, preSorted = TRUE)
  q3 <- .quantileFromFreq(sortedValues, sortedFreq, 0.75, FALSE, preSorted = TRUE)
  absDeviations <- abs(sortedValues - medianValue)
  madOrdering <- order(absDeviations)
  madValue <- .quantileFromFreq(absDeviations[madOrdering], sortedFreq[madOrdering], 0.5, FALSE, preSorted = TRUE)

  list(
    weightedMean = stats$mean,
    weightedSd = stats$sd,
    weightedMedian = medianValue,
    weightedIqr = q3 - q1,
    valueRange = max(values) - min(values),
    weightedMad = madValue
  )
}

.computeStratumTieDiagnostics <- function(df, valueColumn, frequencyColumn) {
  frequencies <- df[[frequencyColumn]]
  total <- sum(frequencies)
  list(
    propMassAtMode = max(frequencies) / total,
    propValuesWithFrequency1 = mean(frequencies == 1),
    propMassInRepeatedValues = sum(frequencies[frequencies > 1]) / total
  )
}

.computeStratumHeapingDiagnostics <- function(df, valueColumn, frequencyColumn) {
  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  total <- sum(frequencies)
  lastDigit <- abs(round(values)) %% 10
  list(
    propEndingIn0 = sum(frequencies[lastDigit == 0]) / total,
    propEndingIn5 = sum(frequencies[lastDigit == 5]) / total,
    propEndingIn0or5 = sum(frequencies[lastDigit %in% c(0, 5)]) / total
  )
}

.computeStratumNeighborReadiness <- function(df, valueColumn, frequencyColumn, kValues) {
  frequencies <- df[[frequencyColumn]]
  totalFrequency <- sum(frequencies)
  distinctValues <- nrow(df)

  readiness <- list(
    distinctToTotalRatio = distinctValues / totalFrequency
  )
  if (!is.numeric(kValues) || length(kValues) < 1 || any(kValues < 1) || any(kValues %% 1 != 0)) {
    return(readiness)
  }

  support <- setNames(
    as.list(vapply(kValues, function(k) distinctValues >= k + 1 && totalFrequency > k, logical(1))),
    paste0("supportsK_", kValues)
  )
  c(readiness, support)
}

.computeStratumCovarianceReadiness <- function(df, featureColumns, frequencyColumn) {
  if (is.null(featureColumns) || length(featureColumns) < 1) {
    return(list(
      featureCount = NA_integer_,
      hasZeroVarianceFeature = NA,
      maxAbsCorrelation = NA_real_
    ))
  }
  if (!all(featureColumns %in% names(df))) {
    stop("Some `featureColumns` are missing from `df`.")
  }
  featuresDf <- df[, featureColumns, drop = FALSE]
  if (!all(vapply(featuresDf, is.numeric, logical(1)))) {
    stop("All `featureColumns` must be numeric.")
  }

  positiveIdx <- which(df[[frequencyColumn]] > 0)
  featuresDf <- featuresDf[positiveIdx, , drop = FALSE]
  zeroVar <- vapply(featuresDf, function(x) stats::var(x) == 0, logical(1))
  corr <- if (ncol(featuresDf) < 2) {
    NA_real_
  } else {
    suppressWarnings(stats::cor(featuresDf, use = "pairwise.complete.obs"))
  }
  maxAbsCorrelation <- if (length(corr) == 1 && is.na(corr)) {
    NA_real_
  } else if (is.matrix(corr)) {
    max(abs(corr[upper.tri(corr)]), na.rm = TRUE)
  } else {
    NA_real_
  }

  list(
    featureCount = ncol(featuresDf),
    hasZeroVarianceFeature = any(zeroVar),
    maxAbsCorrelation = maxAbsCorrelation
  )
}
