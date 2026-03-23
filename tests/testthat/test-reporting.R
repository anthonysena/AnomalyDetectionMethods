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
