#' Plot an anomaly score distribution
#'
#' Creates a weighted histogram of anomaly scores with an optional vertical
#' cutoff line.
#'
#' @param df A data.frame containing score results.
#' @param scoreColumn Name of the numeric anomaly score column in `df`.
#' @param frequencyColumn Name of the frequency column used as histogram
#'   weights.
#' @param cutoffColumn Optional name of a numeric column containing a score
#'   cutoff. If supplied, the first non-missing value is plotted as a reference
#'   line.
#' @param bins Number of histogram bins.
#'
#' @return A ggplot object.
#' @export
plotAnomalyScoreDistribution <- function(df,
                                         scoreColumn,
                                         frequencyColumn = "frequency",
                                         cutoffColumn = NULL,
                                         bins = 30) {
  .validateReportingDf(df, scoreColumn, frequencyColumn, cutoffColumn)
  if (!is.numeric(bins) || length(bins) != 1 || bins < 1 || bins %% 1 != 0) {
    stop("`bins` must be a single positive integer.")
  }

  plotData <- df[!is.na(df[[scoreColumn]]) & df[[frequencyColumn]] > 0, c(scoreColumn, frequencyColumn), drop = FALSE]
  names(plotData) <- c("score", "frequency")

  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$score, weight = .data$frequency)
  ) +
    ggplot2::geom_histogram(bins = bins, fill = "grey70", color = "grey35") +
    ggplot2::labs(x = scoreColumn, y = "Weighted Count")

  if (!is.null(cutoffColumn)) {
    cutoff <- .firstNonMissing(df[[cutoffColumn]])
    if (!is.na(cutoff)) {
      p <- p + ggplot2::geom_vline(xintercept = cutoff, linetype = 2, color = "red")
    }
  }

  p
}

#' Plot sorted anomaly scores by rank
#'
#' Sorts rows by anomaly score and plots the ordered scores, with point size
#' proportional to `frequency`.
#'
#' @param df A data.frame containing score results.
#' @param scoreColumn Name of the numeric anomaly score column in `df`.
#' @param frequencyColumn Name of the frequency column used for point size.
#' @param cutoffColumn Optional name of a numeric column containing a score
#'   cutoff.
#' @param decreasing Logical; whether to sort scores from largest to smallest.
#'
#' @return A ggplot object.
#' @export
plotAnomalyRank <- function(df,
                            scoreColumn,
                            frequencyColumn = "frequency",
                            cutoffColumn = NULL,
                            decreasing = TRUE) {
  .validateReportingDf(df, scoreColumn, frequencyColumn, cutoffColumn)
  if (!is.logical(decreasing) || length(decreasing) != 1) {
    stop("`decreasing` must be TRUE or FALSE.")
  }

  plotData <- df[!is.na(df[[scoreColumn]]) & df[[frequencyColumn]] > 0, c(scoreColumn, frequencyColumn), drop = FALSE]
  ordering <- order(plotData[[scoreColumn]], decreasing = decreasing)
  plotData <- plotData[ordering, , drop = FALSE]
  plotData$rank <- seq_len(nrow(plotData))
  names(plotData)[1:2] <- c("score", "frequency")

  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$rank, y = .data$score)
  ) +
    ggplot2::geom_line(color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(size = .data$frequency), color = "steelblue") +
    ggplot2::labs(x = "Rank", y = scoreColumn, size = "Frequency")

  if (!is.null(cutoffColumn)) {
    cutoff <- .firstNonMissing(df[[cutoffColumn]])
    if (!is.na(cutoff)) {
      p <- p + ggplot2::geom_hline(yintercept = cutoff, linetype = 2, color = "red")
    }
  }

  p
}

#' Plot anomaly results in feature space
#'
#' Creates a two-feature scatter plot with point size mapped to `frequency` and
#' optional color mapping to a score column or outlier flag.
#'
#' @param df A data.frame containing anomaly results.
#' @param xColumn Name of the feature column for the x-axis.
#' @param yColumn Name of the feature column for the y-axis.
#' @param frequencyColumn Name of the frequency column used for point size.
#' @param scoreColumn Optional name of a numeric score column used for color.
#' @param flagColumn Name of the logical outlier flag column used for color
#'   when `scoreColumn` is `NULL`.
#'
#' @return A ggplot object.
#' @export
plotAnomalyFeatureScatter <- function(df,
                                      xColumn,
                                      yColumn,
                                      frequencyColumn = "frequency",
                                      scoreColumn = NULL,
                                      flagColumn = "isOutlier") {
  .validateScatterDf(df, xColumn, yColumn, frequencyColumn, scoreColumn, flagColumn)

  plotData <- df[df[[frequencyColumn]] > 0, c(xColumn, yColumn, frequencyColumn), drop = FALSE]
  names(plotData)[1:3] <- c("x", "y", "frequency")

  if (!is.null(scoreColumn)) {
    plotData$colorValue <- df[df[[frequencyColumn]] > 0, scoreColumn]
    p <- ggplot2::ggplot(
      plotData,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$frequency, color = .data$colorValue)
    ) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_gradient(low = "grey70", high = "red") +
      ggplot2::labs(x = xColumn, y = yColumn, size = "Frequency", color = scoreColumn)
  } else {
    plotData$colorValue <- df[df[[frequencyColumn]] > 0, flagColumn]
    p <- ggplot2::ggplot(
      plotData,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$frequency, color = .data$colorValue)
    ) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey55", `TRUE` = "red")) +
      ggplot2::labs(x = xColumn, y = yColumn, size = "Frequency", color = flagColumn)
  }

  p
}

#' Plot anomaly results on principal components
#'
#' Computes a principal component projection from the selected feature columns
#' and plots the first two principal components, with point size proportional
#' to `frequency`.
#'
#' @param df A data.frame containing anomaly results and feature columns.
#' @param featureColumns Character vector of numeric feature column names.
#' @param frequencyColumn Name of the frequency column used for point size.
#' @param scoreColumn Optional name of a numeric score column used for color.
#' @param flagColumn Name of the logical outlier flag column used for color
#'   when `scoreColumn` is `NULL`.
#' @param center Logical; passed to `stats::prcomp()`.
#' @param scale. Logical; passed to `stats::prcomp()`.
#'
#' @return A ggplot object.
#' @export
plotAnomalyPcaProjection <- function(df,
                                     featureColumns,
                                     frequencyColumn = "frequency",
                                     scoreColumn = NULL,
                                     flagColumn = "isOutlier",
                                     center = TRUE,
                                     scale. = FALSE) {
  .validateFeatureColumnsDf(df, featureColumns, frequencyColumn, scoreColumn, flagColumn)
  if (!is.logical(center) || length(center) != 1) {
    stop("`center` must be TRUE or FALSE.")
  }
  if (!is.logical(scale.) || length(scale.) != 1) {
    stop("`scale.` must be TRUE or FALSE.")
  }

  plotIdx <- which(df[[frequencyColumn]] > 0)
  x <- as.matrix(df[plotIdx, featureColumns, drop = FALSE])
  if (ncol(x) < 2) {
    stop("Need at least 2 feature columns to compute a PCA projection.")
  }
  if (anyNA(x)) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  pca <- stats::prcomp(x, center = center, scale. = scale.)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  names(scores) <- c("PC1", "PC2")
  scores$frequency <- df[[frequencyColumn]][plotIdx]

  if (!is.null(scoreColumn)) {
    scores$colorValue <- df[[scoreColumn]][plotIdx]
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(x = .data$PC1, y = .data$PC2, size = .data$frequency, color = .data$colorValue)
    ) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_gradient(low = "grey70", high = "red") +
      ggplot2::labs(x = "PC1", y = "PC2", size = "Frequency", color = scoreColumn)
  } else {
    scores$colorValue <- df[[flagColumn]][plotIdx]
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(x = .data$PC1, y = .data$PC2, size = .data$frequency, color = .data$colorValue)
    ) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey55", `TRUE` = "red")) +
      ggplot2::labs(x = "PC1", y = "PC2", size = "Frequency", color = flagColumn)
  }

  p
}

#' Plot a chi-square Q-Q diagnostic for Gaussian anomaly detection
#'
#' Creates a Q-Q plot comparing observed squared Mahalanobis distances to the
#' theoretical chi-square quantiles implied by the number of feature columns.
#'
#' @param df A data.frame containing Gaussian anomaly results.
#' @param distanceColumn Name of the numeric Mahalanobis distance column.
#' @param frequencyColumn Name of the frequency column.
#' @param featureColumns Character vector of feature column names used to fit
#'   the Gaussian model.
#' @param addLine Logical; whether to add a 45-degree reference line.
#'
#' @return A ggplot object.
#' @export
plotGadDistanceQQ <- function(df,
                              distanceColumn = "mahalanobisDistance",
                              frequencyColumn = "frequency",
                              featureColumns) {
  .validateFeatureColumnsDf(df, featureColumns, frequencyColumn, distanceColumn, NULL)

  distances <- df[[distanceColumn]]
  frequencies <- df[[frequencyColumn]]
  positiveIdx <- which(frequencies > 0 & !is.na(distances))
  squaredDistances <- distances[positiveIdx] ^ 2
  weights <- frequencies[positiveIdx]
  ordering <- order(squaredDistances)
  squaredDistances <- squaredDistances[ordering]
  weights <- weights[ordering]
  totalN <- sum(weights)
  cumulative <- cumsum(weights)
  probs <- (cumulative - 0.5 * weights) / totalN
  theoretical <- stats::qchisq(probs, df = length(featureColumns))

  plotData <- data.frame(
    theoretical = theoretical,
    observed = squaredDistances,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$theoretical, y = .data$observed)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
    ggplot2::labs(x = "Theoretical Chi-square Quantiles", y = "Observed Squared Distances")
}

#' Summarize LOF sensitivity across neighborhood sizes
#'
#' Runs local outlier factor across multiple `k` values and summarizes each
#' row's score variability and flag frequency.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param featureColumns Character vector of numeric feature column names.
#' @param kValues Positive integer vector of neighborhood sizes.
#' @param valueColumn Name of the value column.
#' @param frequencyColumn Name of the frequency column.
#' @param scoreCutoff LOF score cutoff used for flagging.
#'
#' @return A data.frame with score means, score standard deviations, and flag
#'   rates across the supplied `kValues`.
#' @export
summarizeLofSensitivity <- function(df,
                                    featureColumns,
                                    kValues,
                                    valueColumn = "value",
                                    frequencyColumn = "frequency",
                                    scoreCutoff = 1.5) {
  if (!is.numeric(kValues) || length(kValues) < 1 || any(kValues < 1) || any(kValues %% 1 != 0)) {
    stop("`kValues` must be a non-empty vector of positive integers.")
  }
  results <- lapply(kValues, function(k) {
    localOutlierFactorOutliers(
      df = df,
      valueColumn = valueColumn,
      frequencyColumn = frequencyColumn,
      featureColumns = featureColumns,
      k = k,
      scoreCutoff = scoreCutoff
    )
  })

  scoreMatrix <- do.call(cbind, lapply(results, function(x) x$localOutlierFactor))
  flagMatrix <- do.call(cbind, lapply(results, function(x) x$isOutlier))
  out <- df
  out$meanLofScore <- rowMeans(scoreMatrix, na.rm = TRUE)
  out$sdLofScore <- apply(scoreMatrix, 1, stats::sd, na.rm = TRUE)
  out$flagRate <- rowMeans(flagMatrix, na.rm = TRUE)
  out$kValues <- paste(kValues, collapse = ", ")
  out
}

#' Summarize Isolation Forest stability across random seeds
#'
#' Re-fits Isolation Forest across multiple random seeds and summarizes score
#' variability and flag frequency for each row.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param featureColumns Character vector of numeric feature column names.
#' @param seeds Integer vector of seeds to evaluate.
#' @param ... Additional arguments passed to [isolationForestOutliers()].
#'
#' @return A data.frame with score means, score standard deviations, and flag
#'   rates across the supplied seeds.
#' @export
summarizeIsolationForestStability <- function(df,
                                              featureColumns,
                                              seeds,
                                              ...) {
  if (!is.numeric(seeds) || length(seeds) < 1 || any(seeds %% 1 != 0)) {
    stop("`seeds` must be a non-empty vector of integers.")
  }

  results <- lapply(seeds, function(seed) {
    isolationForestOutliers(
      df = df,
      featureColumns = featureColumns,
      seed = seed,
      ...
    )
  })

  scoreMatrix <- do.call(cbind, lapply(results, function(x) x$isolationForestScore))
  flagMatrix <- do.call(cbind, lapply(results, function(x) x$isOutlier))
  out <- df
  out$meanIsolationForestScore <- rowMeans(scoreMatrix, na.rm = TRUE)
  out$sdIsolationForestScore <- apply(scoreMatrix, 1, stats::sd, na.rm = TRUE)
  out$flagRate <- rowMeans(flagMatrix, na.rm = TRUE)
  out$seeds <- paste(seeds, collapse = ", ")
  out
}

#' Compare anomaly results across multiple detectors
#'
#' Combines several detector result data frames into a single comparison table
#' with one score column, one flag column, and one consensus count per method.
#'
#' @param results Named list of result data.frames.
#' @param idColumns Character vector of columns used to align rows across
#'   methods.
#' @param scoreColumns Optional named character vector mapping method names to
#'   score column names. If `NULL`, default score columns are inferred.
#' @param flagColumn Name of the logical flag column shared across result data
#'   frames.
#'
#' @return A data.frame with per-method score and flag columns plus a
#'   `consensusCount` column.
#' @export
compareAnomalyResults <- function(results,
                                  idColumns = c("value", "frequency"),
                                  scoreColumns = NULL,
                                  flagColumn = "isOutlier") {
  if (!is.list(results) || length(results) < 1 || is.null(names(results)) || any(names(results) == "")) {
    stop("`results` must be a named non-empty list of data.frames.")
  }
  if (!is.character(idColumns) || length(idColumns) < 1) {
    stop("`idColumns` must be a non-empty character vector.")
  }
  if (!is.null(scoreColumns) && (!is.character(scoreColumns) || is.null(names(scoreColumns)))) {
    stop("`scoreColumns` must be NULL or a named character vector.")
  }

  inferredScores <- c(
    lof = "localOutlierFactor",
    isolation_forest = "isolationForestScore",
    gaussian = "mahalanobisDistance"
  )

  merged <- NULL
  flagColumns <- character(0)
  for (methodName in names(results)) {
    resultDf <- results[[methodName]]
    if (!is.data.frame(resultDf)) {
      stop("Each element of `results` must be a data.frame.")
    }
    if (!all(idColumns %in% names(resultDf))) {
      stop(paste0("Result `", methodName, "` is missing required id columns."))
    }
    if (!flagColumn %in% names(resultDf)) {
      stop(paste0("Result `", methodName, "` is missing flag column `", flagColumn, "`."))
    }

    scoreColumn <- if (is.null(scoreColumns)) {
      inferred <- unname(inferredScores[methodName])
      if (is.na(inferred) || !(inferred %in% names(resultDf))) {
        numericColumns <- names(resultDf)[vapply(resultDf, is.numeric, logical(1))]
        numericColumns <- setdiff(numericColumns, idColumns)
        numericColumns <- setdiff(numericColumns, c(flagColumn, "outlierProportion"))
        if (length(numericColumns) < 1) {
          stop(paste0("Could not infer score column for method `", methodName, "`."))
        }
        numericColumns[1]
      } else {
        inferred
      }
    } else {
      scoreColumns[[methodName]]
    }
    if (!scoreColumn %in% names(resultDf)) {
      stop(paste0("Score column `", scoreColumn, "` not found for method `", methodName, "`."))
    }

    methodDf <- resultDf[, c(idColumns, scoreColumn, flagColumn), drop = FALSE]
    names(methodDf)[(length(idColumns) + 1):(length(idColumns) + 2)] <- c(
      paste0(methodName, "_score"),
      paste0(methodName, "_flag")
    )
    flagColumns <- c(flagColumns, paste0(methodName, "_flag"))

    merged <- if (is.null(merged)) {
      methodDf
    } else {
      merge(merged, methodDf, by = idColumns, all = FALSE, sort = FALSE)
    }
  }

  merged$consensusCount <- rowSums(as.data.frame(merged[, flagColumns, drop = FALSE]))
  merged
}

#' Compare univariate outlier detection results across methods
#'
#' Combines several univariate detector result data frames into a single
#' comparison table. For each method, the output includes a standardized
#' severity score, the outlier flag, and any available threshold columns.
#'
#' For threshold-based methods, the severity score is the distance beyond the
#' nearest threshold and is zero for in-threshold values. For score-based
#' methods, the severity score is inferred from the absolute score or method-
#' specific numeric output.
#'
#' @param results Named list of result data.frames.
#' @param idColumns Character vector of columns used to align rows across
#'   methods.
#' @param flagColumn Name of the logical flag column shared across result data
#'   frames.
#'
#' @return A data.frame with per-method score, flag, threshold columns, and a
#'   `consensusCount` column.
#' @export
compareUnivariateOutlierResults <- function(results,
                                            idColumns = c("value", "frequency"),
                                            flagColumn = "isOutlier") {
  if (!is.list(results) || length(results) < 1 || is.null(names(results)) || any(names(results) == "")) {
    stop("`results` must be a named non-empty list of data.frames.")
  }
  if (!is.character(idColumns) || length(idColumns) < 1) {
    stop("`idColumns` must be a non-empty character vector.")
  }

  merged <- NULL
  flagColumns <- character(0)
  for (methodName in names(results)) {
    resultDf <- results[[methodName]]
    if (!is.data.frame(resultDf)) {
      stop("Each element of `results` must be a data.frame.")
    }
    if (!all(idColumns %in% names(resultDf))) {
      stop(paste0("Result `", methodName, "` is missing required id columns."))
    }
    if (!flagColumn %in% names(resultDf)) {
      stop(paste0("Result `", methodName, "` is missing flag column `", flagColumn, "`."))
    }

    components <- .extractUnivariateComparisonComponents(
      resultDf = resultDf,
      methodName = methodName,
      idColumns = idColumns,
      flagColumn = flagColumn
    )

    methodDf <- cbind(
      resultDf[, idColumns, drop = FALSE],
      components,
      stringsAsFactors = FALSE
    )
    flagColumns <- c(flagColumns, paste0(methodName, "_flag"))
    merged <- if (is.null(merged)) {
      methodDf
    } else {
      merge(merged, methodDf, by = idColumns, all = FALSE, sort = FALSE)
    }
  }

  merged$consensusCount <- rowSums(as.data.frame(merged[, flagColumns, drop = FALSE]))
  merged
}

#' Learn plausible value ranges from univariate consensus voting
#'
#' Computes consensus outlier scores from the binary method flags returned by
#' [compareUnivariateOutlierResults()] and derives plausible value ranges from
#' the retained values. The consensus score can be based on a simple vote count,
#' frequency-weighted votes, or support-weighted votes. Plausible ranges can be
#' learned using weighted quantiles or raw min/max values.
#'
#' @param comparisonDf Optional output from [compareUnivariateOutlierResults()].
#' @param results Optional named list of univariate result data.frames. Used only
#'   when `comparisonDf` is `NULL`.
#' @param valueColumn Name of the numeric value column.
#' @param frequencyColumn Name of the numeric frequency column.
#' @param consensusMethod One or more consensus scoring methods. Supported
#'   values are `"count"`, `"frequency_weighted"`, and `"support_weighted"`.
#' @param consensusThreshold Numeric threshold used to flag consensus outliers.
#'   For `"count"`, values with score `>= consensusThreshold` are excluded. For
#'   weighted methods, the threshold is applied to the weighted score scale.
#' @param rangeMethod One or more range derivation methods. Supported values are
#'   `"weighted_quantile"` and `"raw_min_max"`.
#' @param lowerProb Lower quantile probability used when
#'   `rangeMethod = "weighted_quantile"`.
#' @param upperProb Upper quantile probability used when
#'   `rangeMethod = "weighted_quantile"`.
#'
#' @return A list with `valueSummary` and `rangeSummary` data.frames.
#' @export
learnPlausibleRangeConsensus <- function(comparisonDf = NULL,
                                         results = NULL,
                                         valueColumn = "value",
                                         frequencyColumn = "frequency",
                                         consensusMethod = c("count", "frequency_weighted", "support_weighted"),
                                         consensusThreshold = 3,
                                         rangeMethod = c("weighted_quantile", "raw_min_max"),
                                         lowerProb = 0.001,
                                         upperProb = 0.999) {
  if (is.null(comparisonDf)) {
    if (is.null(results)) {
      stop("Provide either `comparisonDf` or `results`.")
    }
    comparisonDf <- compareUnivariateOutlierResults(
      results = results,
      idColumns = c(valueColumn, frequencyColumn)
    )
  }
  if (!is.data.frame(comparisonDf)) {
    stop("`comparisonDf` must be a data.frame.")
  }
  if (!all(c(valueColumn, frequencyColumn) %in% names(comparisonDf))) {
    stop("`comparisonDf` must contain `valueColumn` and `frequencyColumn`.")
  }
  if (!is.numeric(comparisonDf[[valueColumn]]) || !is.numeric(comparisonDf[[frequencyColumn]])) {
    stop("`valueColumn` and `frequencyColumn` must be numeric.")
  }
  if (!is.numeric(consensusThreshold) || length(consensusThreshold) != 1 || consensusThreshold < 0) {
    stop("`consensusThreshold` must be a single non-negative number.")
  }
  if (!is.numeric(lowerProb) || length(lowerProb) != 1 || lowerProb < 0 || lowerProb > 1) {
    stop("`lowerProb` must be a single number in [0, 1].")
  }
  if (!is.numeric(upperProb) || length(upperProb) != 1 || upperProb < 0 || upperProb > 1) {
    stop("`upperProb` must be a single number in [0, 1].")
  }
  if (lowerProb >= upperProb) {
    stop("`lowerProb` must be < `upperProb`.")
  }

  consensusMethod <- match.arg(
    consensusMethod,
    choices = c("count", "frequency_weighted", "support_weighted"),
    several.ok = TRUE
  )
  rangeMethod <- match.arg(
    rangeMethod,
    choices = c("weighted_quantile", "raw_min_max"),
    several.ok = TRUE
  )

  flagColumns <- grep("_flag$", names(comparisonDf), value = TRUE)
  if (length(flagColumns) < 1) {
    stop("`comparisonDf` must contain at least one `*_flag` column.")
  }
  if (!all(vapply(comparisonDf[, flagColumns, drop = FALSE], is.logical, logical(1)))) {
    stop("All `*_flag` columns in `comparisonDf` must be logical.")
  }

  values <- comparisonDf[[valueColumn]]
  frequencies <- comparisonDf[[frequencyColumn]]
  totalFrequency <- sum(frequencies)
  flagMatrix <- as.data.frame(comparisonDf[, flagColumns, drop = FALSE])
  flagCount <- rowSums(flagMatrix)
  support <- frequencies / totalFrequency

  valueSummary <- comparisonDf[, c(valueColumn, frequencyColumn), drop = FALSE]
  valueSummary$consensusCount <- flagCount
  valueSummary$support <- support
  valueSummary$consensusScore_count <- flagCount
  valueSummary$consensusScore_frequency_weighted <- flagCount * frequencies
  valueSummary$consensusScore_support_weighted <- flagCount * support

  rangeRows <- list()
  idx <- 1L
  for (scoring in consensusMethod) {
    scoreColumn <- paste0("consensusScore_", scoring)
    consensusScore <- valueSummary[[scoreColumn]]
    isConsensusOutlier <- consensusScore >= consensusThreshold
    valueSummary[[paste0("isConsensusOutlier_", scoring)]] <- isConsensusOutlier

    retainedIdx <- which(!isConsensusOutlier)
    if (length(retainedIdx) < 1) {
      stop(paste0(
        "No retained values remain under `consensusMethod = \"", scoring,
        "\"` and `consensusThreshold = ", consensusThreshold, "`."
      ))
    }

    retainedValues <- values[retainedIdx]
    retainedFrequencies <- frequencies[retainedIdx]
    for (method in rangeMethod) {
      bounds <- .deriveConsensusBounds(
        values = retainedValues,
        frequencies = retainedFrequencies,
        rangeMethod = method,
        lowerProb = lowerProb,
        upperProb = upperProb
      )
      rangeRows[[idx]] <- data.frame(
        consensusMethod = scoring,
        consensusThreshold = consensusThreshold,
        rangeMethod = method,
        lowerProb = if (method == "weighted_quantile") lowerProb else NA_real_,
        upperProb = if (method == "weighted_quantile") upperProb else NA_real_,
        minPlausible = bounds[1],
        maxPlausible = bounds[2],
        retainedDistinctValues = length(retainedValues),
        retainedTotalFrequency = sum(retainedFrequencies),
        excludedDistinctValues = sum(isConsensusOutlier),
        excludedTotalFrequency = sum(frequencies[isConsensusOutlier]),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  list(
    valueSummary = valueSummary,
    rangeSummary = do.call(rbind, rangeRows)
  )
}

.validateReportingDf <- function(df, scoreColumn, frequencyColumn, cutoffColumn = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  required <- c(scoreColumn, frequencyColumn)
  if (!is.null(cutoffColumn)) {
    required <- c(required, cutoffColumn)
  }
  if (!all(required %in% names(df))) {
    stop("`df` is missing required reporting columns.")
  }
  if (!is.numeric(df[[scoreColumn]]) || !is.numeric(df[[frequencyColumn]])) {
    stop("`scoreColumn` and `frequencyColumn` must be numeric.")
  }
  if (!is.null(cutoffColumn) && !is.numeric(df[[cutoffColumn]])) {
    stop("`cutoffColumn` must be numeric.")
  }
}

.validateScatterDf <- function(df, xColumn, yColumn, frequencyColumn, scoreColumn, flagColumn) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  required <- c(xColumn, yColumn, frequencyColumn)
  if (!is.null(scoreColumn)) {
    required <- c(required, scoreColumn)
  } else {
    required <- c(required, flagColumn)
  }
  if (!all(required %in% names(df))) {
    stop("`df` is missing required plotting columns.")
  }
  if (!is.numeric(df[[xColumn]]) || !is.numeric(df[[yColumn]]) || !is.numeric(df[[frequencyColumn]])) {
    stop("`xColumn`, `yColumn`, and `frequencyColumn` must be numeric.")
  }
  if (!is.null(scoreColumn) && !is.numeric(df[[scoreColumn]])) {
    stop("`scoreColumn` must be numeric.")
  }
  if (is.null(scoreColumn) && !is.logical(df[[flagColumn]])) {
    stop("`flagColumn` must be logical when `scoreColumn` is NULL.")
  }
}

.validateFeatureColumnsDf <- function(df, featureColumns, frequencyColumn, scoreColumn = NULL, flagColumn = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  if (!is.character(featureColumns) || length(featureColumns) < 1) {
    stop("`featureColumns` must be a non-empty character vector.")
  }
  required <- c(featureColumns, frequencyColumn)
  if (!is.null(scoreColumn)) {
    required <- c(required, scoreColumn)
  }
  if (!is.null(flagColumn)) {
    required <- c(required, flagColumn)
  }
  if (!all(required %in% names(df))) {
    stop("`df` is missing required feature or result columns.")
  }
  if (!all(vapply(df[, featureColumns, drop = FALSE], is.numeric, logical(1)))) {
    stop("All `featureColumns` must be numeric.")
  }
}

.firstNonMissing <- function(x) {
  idx <- which(!is.na(x))[1]
  if (is.na(idx)) {
    return(NA_real_)
  }
  x[idx]
}

.extractUnivariateComparisonComponents <- function(resultDf,
                                                   methodName,
                                                   idColumns,
                                                   flagColumn) {
  valueColumn <- idColumns[1]
  values <- resultDf[[valueColumn]]
  scoreName <- paste0(methodName, "_score")
  flagName <- paste0(methodName, "_flag")
  components <- list()

  if (all(c("lowerFence", "upperFence") %in% names(resultDf))) {
    severity <- .thresholdSeverity(values, resultDf$lowerFence, resultDf$upperFence)
    components[[scoreName]] <- severity
    components[[paste0(methodName, "_lowerThreshold")]] <- resultDf$lowerFence
    components[[paste0(methodName, "_upperThreshold")]] <- resultDf$upperFence
  } else if (all(c("lowerThreshold", "upperThreshold") %in% names(resultDf))) {
    severity <- .thresholdSeverity(values, resultDf$lowerThreshold, resultDf$upperThreshold)
    components[[scoreName]] <- severity
    components[[paste0(methodName, "_lowerThreshold")]] <- resultDf$lowerThreshold
    components[[paste0(methodName, "_upperThreshold")]] <- resultDf$upperThreshold
  } else if ("zScore" %in% names(resultDf)) {
    components[[scoreName]] <- abs(resultDf$zScore)
  } else if ("modifiedZScore" %in% names(resultDf)) {
    components[[scoreName]] <- abs(resultDf$modifiedZScore)
  } else if ("outlierCount" %in% names(resultDf)) {
    components[[scoreName]] <- resultDf$outlierCount
  } else {
    numericColumns <- names(resultDf)[vapply(resultDf, is.numeric, logical(1))]
    numericColumns <- setdiff(numericColumns, c(idColumns, "outlierProportion"))
    numericColumns <- setdiff(numericColumns, c("lowerFence", "upperFence", "lowerThreshold", "upperThreshold"))
    if (length(numericColumns) < 1) {
      stop(paste0("Could not infer a score-like column for method `", methodName, "`."))
    }
    components[[scoreName]] <- resultDf[[numericColumns[1]]]
  }

  components[[flagName]] <- resultDf[[flagColumn]]
  as.data.frame(components, stringsAsFactors = FALSE)
}

.thresholdSeverity <- function(values, lower, upper) {
  below <- pmax(lower - values, 0)
  above <- pmax(values - upper, 0)
  pmax(below, above)
}

.deriveConsensusBounds <- function(values, frequencies, rangeMethod, lowerProb, upperProb) {
  if (rangeMethod == "raw_min_max") {
    return(c(min(values), max(values)))
  }

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]
  c(
    .quantileFromFreq(values, frequencies, lowerProb, interpolate = FALSE, preSorted = TRUE),
    .quantileFromFreq(values, frequencies, upperProb, interpolate = FALSE, preSorted = TRUE)
  )
}
