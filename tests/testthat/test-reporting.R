test_that("reporting plots return ggplot objects", {
  skip_if_not_installed("isotree")

  df <- data.frame(
    value = c(1, 2, 3, 4, 5, 6),
    frequency = c(3, 3, 3, 3, 3, 1),
    featureA = c(0.0, 0.1, -0.1, 4.0, 4.1, 8.0),
    featureB = c(0.0, -0.1, 0.1, 4.2, 3.9, 8.1),
    featureC = c(1.0, 1.1, 0.9, 1.0, 1.2, 2.5)
  )

  lofOut <- localOutlierFactorOutliers(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    k = 3,
    scoreCutoff = 1.2
  )
  ifOut <- isolationForestOutliers(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    ntrees = 100,
    sampleSize = 6,
    scoreCutoff = 0.6,
    seed = 123,
    nthreads = 1
  )
  gadOut <- gaussianAnomalyOutliers(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    tailProb = 0.99
  )

  expect_s3_class(
    plotAnomalyScoreDistribution(lofOut, "localOutlierFactor", cutoffColumn = "scoreCutoff"),
    "ggplot"
  )
  expect_s3_class(
    plotAnomalyRank(ifOut, "isolationForestScore", cutoffColumn = "scoreCutoff"),
    "ggplot"
  )
  expect_s3_class(
    plotAnomalyFeatureScatter(lofOut, "featureA", "featureB", scoreColumn = "localOutlierFactor"),
    "ggplot"
  )
  expect_s3_class(
    plotAnomalyPcaProjection(gadOut, c("featureA", "featureB", "featureC"), scoreColumn = "mahalanobisDistance"),
    "ggplot"
  )
  expect_s3_class(
    plotGadDistanceQQ(gadOut, featureColumns = c("featureA", "featureB", "featureC")),
    "ggplot"
  )
})

test_that("reporting summaries and comparison tables return expected fields", {
  skip_if_not_installed("isotree")

  df <- data.frame(
    value = c(1, 2, 3, 4, 5, 6),
    frequency = c(3, 3, 3, 3, 3, 1),
    featureA = c(0.0, 0.1, -0.1, 4.0, 4.1, 8.0),
    featureB = c(0.0, -0.1, 0.1, 4.2, 3.9, 8.1),
    featureC = c(1.0, 1.1, 0.9, 1.0, 1.2, 2.5)
  )

  lofOut <- localOutlierFactorOutliers(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    k = 3,
    scoreCutoff = 1.2
  )
  ifOut <- isolationForestOutliers(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    ntrees = 100,
    sampleSize = 6,
    scoreCutoff = 0.6,
    seed = 123,
    nthreads = 1
  )
  gadOut <- gaussianAnomalyOutliers(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    tailProb = 0.99
  )

  lofSensitivity <- summarizeLofSensitivity(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    kValues = c(2, 3, 4),
    scoreCutoff = 1.2
  )
  ifStability <- summarizeIsolationForestStability(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    seeds = c(1, 2, 3),
    ntrees = 100,
    sampleSize = 6,
    scoreCutoff = 0.6,
    nthreads = 1
  )
  comparison <- compareAnomalyResults(
    list(
      lof = lofOut,
      isolation_forest = ifOut,
      gaussian = gadOut
    )
  )

  expect_true(all(c("meanLofScore", "sdLofScore", "flagRate") %in% names(lofSensitivity)))
  expect_true(all(c("meanIsolationForestScore", "sdIsolationForestScore", "flagRate") %in% names(ifStability)))
  expect_true(all(c(
    "lof_score",
    "lof_flag",
    "isolation_forest_score",
    "isolation_forest_flag",
    "gaussian_score",
    "gaussian_flag",
    "consensusCount"
  ) %in% names(comparison)))
  expect_equal(nrow(comparison), nrow(df))
})

test_that("compareUnivariateOutlierResults returns score, flag, and threshold fields", {
  df <- data.frame(
    value = c(10, 11, 12, 13, 100),
    frequency = c(8, 10, 10, 8, 1)
  )

  tukeyOut <- tukeyFences(df)
  quantileOut <- quantileThresholds(df, lowerProb = 0.05, upperProb = 0.95)
  zOut <- zScoreOutliers(df, zCutoff = 2)
  modifiedOut <- modifiedZScoreOutliers(df, zCutoff = 2.5)
  esdOut <- generalizedESDOutliers(df, maxOutliers = 3, alpha = 0.05)

  comparison <- compareUnivariateOutlierResults(
    list(
      tukey = tukeyOut,
      quantile = quantileOut,
      zscore = zOut,
      modified = modifiedOut,
      esd = esdOut
    )
  )

  expect_true(all(c(
    "tukey_score",
    "tukey_flag",
    "tukey_lowerThreshold",
    "tukey_upperThreshold",
    "quantile_score",
    "quantile_flag",
    "quantile_lowerThreshold",
    "quantile_upperThreshold",
    "zscore_score",
    "zscore_flag",
    "modified_score",
    "modified_flag",
    "esd_score",
    "esd_flag",
    "consensusCount"
  ) %in% names(comparison)))
  expect_equal(nrow(comparison), nrow(df))
  expect_true(comparison$tukey_score[comparison$value == 100] > 0)
  expect_equal(comparison$zscore_score, abs(zOut$zScore))
  expect_equal(comparison$modified_score, abs(modifiedOut$modifiedZScore))
  expect_equal(comparison$esd_score, esdOut$outlierCount)
})

test_that("learnPlausibleRangeConsensus returns value and range summaries", {
  df <- data.frame(
    value = c(10, 11, 12, 13, 100),
    frequency = c(8, 10, 10, 8, 1)
  )

  comparison <- compareUnivariateOutlierResults(
    list(
      tukey = tukeyFences(df),
      quantile = quantileThresholds(df, lowerProb = 0.05, upperProb = 0.95),
      zscore = zScoreOutliers(df, zCutoff = 2),
      modified = modifiedZScoreOutliers(df, zCutoff = 2.5),
      esd = generalizedESDOutliers(df, maxOutliers = 3, alpha = 0.05)
    )
  )

  out <- learnPlausibleRangeConsensus(
    comparisonDf = comparison,
    consensusMethod = c("count", "support_weighted"),
    consensusThreshold = 3,
    rangeMethod = c("weighted_quantile", "raw_min_max"),
    lowerProb = 0.05,
    upperProb = 0.95
  )

  expect_true(all(c("valueSummary", "rangeSummary") %in% names(out)))
  expect_true(all(c(
    "consensusCount",
    "support",
    "consensusScore_count",
    "consensusScore_frequency_weighted",
    "consensusScore_support_weighted",
    "isConsensusOutlier_count",
    "isConsensusOutlier_support_weighted"
  ) %in% names(out$valueSummary)))
  expect_true(all(c(
    "consensusMethod",
    "rangeMethod",
    "minPlausible",
    "maxPlausible",
    "retainedDistinctValues",
    "retainedTotalFrequency"
  ) %in% names(out$rangeSummary)))
  expect_true(out$valueSummary$isConsensusOutlier_count[out$valueSummary$value == 100])
  countRanges <- subset(out$rangeSummary, consensusMethod == "count")
  expect_true(all(countRanges$maxPlausible <= 13))
})
