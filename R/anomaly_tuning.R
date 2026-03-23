#' Evaluate LOF settings across a parameter grid
#'
#' Runs local outlier factor across combinations of neighborhood sizes and
#' score cutoffs, then summarizes score ranking and flagging behavior.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param featureColumns Character vector of numeric feature column names.
#' @param kValues Positive integer vector of neighborhood sizes.
#' @param scoreCutoffs Positive numeric vector of LOF cutoffs.
#' @param topN Number of highest-scoring rows used for top-set stability.
#' @param valueColumn Name of the value column.
#' @param frequencyColumn Name of the frequency column.
#'
#' @return A list with `summary` and `stability` data.frames.
#' @export
evaluateLofGrid <- function(df,
                            featureColumns,
                            kValues,
                            scoreCutoffs,
                            topN = 5,
                            valueColumn = "value",
                            frequencyColumn = "frequency") {
  .validateGridInputs(df, featureColumns, valueColumn, frequencyColumn, topN)
  if (!is.numeric(kValues) || length(kValues) < 1 || any(kValues < 1) || any(kValues %% 1 != 0)) {
    stop("`kValues` must be a non-empty vector of positive integers.")
  }
  if (!is.numeric(scoreCutoffs) || length(scoreCutoffs) < 1 || any(scoreCutoffs <= 0)) {
    stop("`scoreCutoffs` must be a non-empty vector of positive numbers.")
  }

  baseResults <- lapply(kValues, function(k) {
    localOutlierFactorOutliers(
      df = df,
      valueColumn = valueColumn,
      frequencyColumn = frequencyColumn,
      featureColumns = featureColumns,
      k = k,
      scoreCutoff = 1
    )
  })
  names(baseResults) <- paste0("k=", kValues)

  summaryRows <- lapply(seq_along(kValues), function(i) {
    baseOut <- baseResults[[i]]
    lapply(scoreCutoffs, function(cutoff) {
      .summarizeGridRun(
        df = baseOut,
        methodLabel = "lof",
        config = c(k = kValues[i], scoreCutoff = cutoff),
        score = baseOut$localOutlierFactor,
        flag = baseOut$localOutlierFactor > cutoff,
        frequency = df[[frequencyColumn]],
        values = df[[valueColumn]],
        topN = topN
      )
    })
  })
  summaryDf <- do.call(rbind, unlist(summaryRows, recursive = FALSE))

  stabilityDf <- .pairwiseGridStability(
    scoreList = lapply(baseResults, function(x) x$localOutlierFactor),
    configLabels = names(baseResults),
    values = df[[valueColumn]],
    topN = topN
  )

  list(summary = summaryDf, stability = stabilityDf)
}

#' Evaluate Isolation Forest settings across a parameter grid
#'
#' Runs Isolation Forest across combinations of tree counts, sample sizes, and
#' random seeds, then summarizes score stability and cutoff sensitivity.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param featureColumns Character vector of numeric feature column names.
#' @param ntreesValues Positive integer vector of tree counts.
#' @param sampleSizeValues Positive integer vector of sample sizes.
#' @param scoreCutoffs Positive numeric vector of score cutoffs.
#' @param seeds Integer vector of random seeds.
#' @param topN Number of highest-scoring rows used for top-set stability.
#' @param ... Additional arguments passed to [isolationForestOutliers()].
#'
#' @return A list with `summary` and `stability` data.frames.
#' @export
evaluateIsolationForestGrid <- function(df,
                                        featureColumns,
                                        ntreesValues,
                                        sampleSizeValues,
                                        scoreCutoffs,
                                        seeds = c(1, 11, 23),
                                        topN = 5,
                                        ...) {
  .validateGridInputs(df, featureColumns, "value", "frequency", topN)
  if (!is.numeric(ntreesValues) || length(ntreesValues) < 1 || any(ntreesValues < 1) || any(ntreesValues %% 1 != 0)) {
    stop("`ntreesValues` must be a non-empty vector of positive integers.")
  }
  if (!is.numeric(sampleSizeValues) || length(sampleSizeValues) < 1 || any(sampleSizeValues < 1) || any(sampleSizeValues %% 1 != 0)) {
    stop("`sampleSizeValues` must be a non-empty vector of positive integers.")
  }
  if (!is.numeric(scoreCutoffs) || length(scoreCutoffs) < 1 || any(scoreCutoffs <= 0)) {
    stop("`scoreCutoffs` must be a non-empty vector of positive numbers.")
  }
  if (!is.numeric(seeds) || length(seeds) < 1 || any(seeds %% 1 != 0)) {
    stop("`seeds` must be a non-empty vector of integers.")
  }

  configGrid <- expand.grid(
    ntrees = ntreesValues,
    sampleSize = sampleSizeValues,
    stringsAsFactors = FALSE
  )

  summaryList <- vector("list", nrow(configGrid) * length(scoreCutoffs))
  stabilityList <- vector("list", nrow(configGrid))
  idxSummary <- 1L

  for (i in seq_len(nrow(configGrid))) {
    config <- configGrid[i, , drop = FALSE]
    seedResults <- lapply(seeds, function(seed) {
      isolationForestOutliers(
        df = df,
        featureColumns = featureColumns,
        ntrees = config$ntrees,
        sampleSize = config$sampleSize,
        scoreCutoff = 0.5,
        seed = seed,
        ...
      )
    })
    names(seedResults) <- paste0("seed=", seeds)
    scoreMatrix <- do.call(cbind, lapply(seedResults, function(x) x$isolationForestScore))
    meanScores <- rowMeans(scoreMatrix, na.rm = TRUE)
    meanScoreSd <- mean(apply(scoreMatrix, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
    meanPairwiseSpearman <- .meanPairwiseScoreCorrelation(scoreMatrix)
    meanPairwiseJaccard <- .meanPairwiseTopJaccard(scoreMatrix, df$value, topN)

    for (cutoff in scoreCutoffs) {
      flagMatrix <- scoreMatrix > cutoff
      meanFlags <- rowMeans(flagMatrix, na.rm = TRUE)
      summaryList[[idxSummary]] <- data.frame(
        method = "isolation_forest",
        ntrees = config$ntrees,
        sampleSize = config$sampleSize,
        scoreCutoff = cutoff,
        meanFlagCount = mean(colSums(flagMatrix)),
        meanWeightedFlagCount = mean(colSums(flagMatrix * df$frequency)),
        meanRowScoreSd = meanScoreSd,
        meanPairwiseSpearman = meanPairwiseSpearman,
        meanPairwiseTopJaccard = meanPairwiseJaccard,
        topNValues = paste(df$value[order(meanScores, decreasing = TRUE)][seq_len(min(topN, length(meanScores)))], collapse = ", "),
        stringsAsFactors = FALSE
      )
      idxSummary <- idxSummary + 1L
    }

    stabilityList[[i]] <- .pairwiseGridStability(
      scoreList = lapply(seedResults, function(x) x$isolationForestScore),
      configLabels = paste0("ntrees=", config$ntrees, ",sampleSize=", config$sampleSize, ",", names(seedResults)),
      values = df$value,
      topN = topN
    )
  }

  list(
    summary = do.call(rbind, summaryList),
    stability = do.call(rbind, stabilityList)
  )
}

#' Evaluate Gaussian anomaly detection settings across a parameter grid
#'
#' Runs Gaussian anomaly detection across combinations of preprocessing,
#' covariance estimation, and tail cutoffs, then summarizes outlier flags and
#' Gaussian fit diagnostics.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param featureColumns Character vector of numeric feature column names.
#' @param transformMethods Character vector of transform methods. Use `NA` or
#'   `"none"` to indicate no preprocessing.
#' @param covarianceMethods Character vector of covariance methods.
#' @param tailProbValues Numeric vector of upper-tail probabilities in `(0, 1)`.
#' @param topN Number of highest-distance rows used for top-set stability.
#' @param ... Additional arguments passed to [gaussianAnomalyOutliers()].
#'
#' @return A list with `summary` and `stability` data.frames.
#' @export
evaluateGadGrid <- function(df,
                            featureColumns,
                            transformMethods = c(NA, "weighted_zscore", "log1p"),
                            covarianceMethods = c("classical", "mcd"),
                            tailProbValues = c(0.95, 0.99),
                            topN = 5,
                            ...) {
  .validateGridInputs(df, featureColumns, "value", "frequency", topN)
  if (!is.character(covarianceMethods) || length(covarianceMethods) < 1) {
    stop("`covarianceMethods` must be a non-empty character vector.")
  }
  if (!is.numeric(tailProbValues) || length(tailProbValues) < 1 || any(tailProbValues <= 0) || any(tailProbValues >= 1)) {
    stop("`tailProbValues` must be a non-empty numeric vector in (0, 1).")
  }

  transformLabels <- ifelse(is.na(transformMethods) | transformMethods == "none", NA_character_, transformMethods)
  configGrid <- expand.grid(
    transformMethod = seq_along(transformLabels),
    covarianceMethod = covarianceMethods,
    tailProb = tailProbValues,
    stringsAsFactors = FALSE
  )

  resultCache <- vector("list", nrow(configGrid))
  summaryRows <- vector("list", nrow(configGrid))
  for (i in seq_len(nrow(configGrid))) {
    transformLabel <- transformLabels[configGrid$transformMethod[i]]
    transformMethod <- transformLabel
    if (is.na(transformMethod)) {
      transformMethod <- NULL
    }
    out <- gaussianAnomalyOutliers(
      df = df,
      featureColumns = featureColumns,
      transformMethod = transformMethod,
      covarianceMethod = configGrid$covarianceMethod[i],
      tailProb = configGrid$tailProb[i],
      ...
    )
    resultCache[[i]] <- out
    qqCorrelation <- .gadQqCorrelation(
      distances = out$mahalanobisDistance,
      frequencies = df$frequency,
      dfDegrees = length(featureColumns)
    )
    summaryRows[[i]] <- .summarizeGridRun(
      df = out,
      methodLabel = "gaussian",
      config = c(
        transformMethod = ifelse(is.na(transformLabel), "none", transformLabel),
        covarianceMethod = configGrid$covarianceMethod[i],
        tailProb = configGrid$tailProb[i]
      ),
      score = out$mahalanobisDistance,
      flag = out$isOutlier,
      frequency = df$frequency,
      values = df$value,
      topN = topN,
      extra = c(qqCorrelation = qqCorrelation)
    )
  }

  stabilityDf <- .pairwiseGridStability(
    scoreList = lapply(resultCache, function(x) x$mahalanobisDistance),
    configLabels = vapply(seq_len(nrow(configGrid)), function(i) {
      paste0(
        "transform=", ifelse(is.na(transformLabels[configGrid$transformMethod[i]]), "none", transformLabels[configGrid$transformMethod[i]]),
        ",cov=", configGrid$covarianceMethod[i],
        ",tail=", configGrid$tailProb[i]
      )
    }, character(1)),
    values = df$value,
    topN = topN
  )

  list(
    summary = do.call(rbind, summaryRows),
    stability = stabilityDf
  )
}

.validateGridInputs <- function(df, featureColumns, valueColumn, frequencyColumn, topN) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }
  if (!all(c(valueColumn, frequencyColumn, featureColumns) %in% names(df))) {
    stop("`df` is missing required columns.")
  }
  if (!is.numeric(topN) || length(topN) != 1 || topN < 1 || topN %% 1 != 0) {
    stop("`topN` must be a single positive integer.")
  }
}

.summarizeGridRun <- function(df,
                              methodLabel,
                              config,
                              score,
                              flag,
                              frequency,
                              values,
                              topN,
                              extra = NULL) {
  topIdx <- order(score, decreasing = TRUE)[seq_len(min(topN, sum(!is.na(score))))]
  row <- c(
    list(
      method = methodLabel
    ),
    as.list(config),
    list(
      flagCount = sum(flag, na.rm = TRUE),
      weightedFlagCount = sum(frequency[flag], na.rm = TRUE),
      topNValues = paste(values[topIdx], collapse = ", "),
      scoreMean = mean(score, na.rm = TRUE),
      scoreSd = stats::sd(score, na.rm = TRUE)
    ),
    as.list(extra)
  )
  as.data.frame(row, stringsAsFactors = FALSE)
}

.pairwiseGridStability <- function(scoreList, configLabels, values, topN) {
  if (length(scoreList) < 2) {
    return(data.frame(
      configA = character(0),
      configB = character(0),
      spearman = numeric(0),
      topNJaccard = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  combos <- utils::combn(seq_along(scoreList), 2)
  rows <- lapply(seq_len(ncol(combos)), function(i) {
    a <- combos[1, i]
    b <- combos[2, i]
    topA <- values[order(scoreList[[a]], decreasing = TRUE)[seq_len(min(topN, length(scoreList[[a]])))]]
    topB <- values[order(scoreList[[b]], decreasing = TRUE)[seq_len(min(topN, length(scoreList[[b]])))]]
    data.frame(
      configA = configLabels[a],
      configB = configLabels[b],
      spearman = suppressWarnings(stats::cor(scoreList[[a]], scoreList[[b]], method = "spearman", use = "pairwise.complete.obs")),
      topNJaccard = .jaccard(topA, topB),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.meanPairwiseScoreCorrelation <- function(scoreMatrix) {
  if (ncol(scoreMatrix) < 2) {
    return(NA_real_)
  }
  combos <- utils::combn(seq_len(ncol(scoreMatrix)), 2)
  cors <- apply(combos, 2, function(idx) {
    suppressWarnings(stats::cor(scoreMatrix[, idx[1]], scoreMatrix[, idx[2]], method = "spearman", use = "pairwise.complete.obs"))
  })
  mean(cors, na.rm = TRUE)
}

.meanPairwiseTopJaccard <- function(scoreMatrix, values, topN) {
  if (ncol(scoreMatrix) < 2) {
    return(NA_real_)
  }
  combos <- utils::combn(seq_len(ncol(scoreMatrix)), 2)
  jaccards <- apply(combos, 2, function(idx) {
    topA <- values[order(scoreMatrix[, idx[1]], decreasing = TRUE)[seq_len(min(topN, nrow(scoreMatrix)))]]
    topB <- values[order(scoreMatrix[, idx[2]], decreasing = TRUE)[seq_len(min(topN, nrow(scoreMatrix)))]]
    .jaccard(topA, topB)
  })
  mean(jaccards, na.rm = TRUE)
}

.gadQqCorrelation <- function(distances, frequencies, dfDegrees) {
  positiveIdx <- which(frequencies > 0 & !is.na(distances))
  squaredDistances <- distances[positiveIdx] ^ 2
  weights <- frequencies[positiveIdx]
  ordering <- order(squaredDistances)
  squaredDistances <- squaredDistances[ordering]
  weights <- weights[ordering]
  probs <- (cumsum(weights) - 0.5 * weights) / sum(weights)
  theoretical <- stats::qchisq(probs, df = dfDegrees)
  suppressWarnings(stats::cor(theoretical, squaredDistances, method = "spearman"))
}

.jaccard <- function(x, y) {
  x <- unique(x)
  y <- unique(y)
  unionLength <- length(union(x, y))
  if (unionLength == 0) {
    return(NA_real_)
  }
  length(intersect(x, y)) / unionLength
}
