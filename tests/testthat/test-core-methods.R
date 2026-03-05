test_that("summarizeValueFrequency returns expected columns", {
  df <- data.frame(
    value = c(1, 2, 3, 4, 100),
    frequency = c(10, 8, 6, 4, 1)
  )

  out <- summarizeValueFrequency(df)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("n", "mean", "median", "mode", "q1", "q3", "iqr"))
  expect_equal(out$n, sum(df$frequency))
  expect_true(out$q3 >= out$q1)
})

test_that("tukeyFences flags obvious high outlier", {
  df <- data.frame(
    value = c(10, 11, 12, 13, 100),
    frequency = c(8, 10, 10, 8, 1)
  )

  out <- tukeyFences(df)

  expect_true(is.logical(out$isOutlier))
  expect_true(out$isOutlier[df$value == 100])
})

test_that("quantileThresholds validates probability arguments", {
  df <- data.frame(value = c(1, 2, 3), frequency = c(1, 1, 1))

  expect_error(quantileThresholds(df, lowerProb = 0.8, upperProb = 0.2), "`lowerProb` must be < `upperProb`.")
  expect_error(quantileThresholds(df, lowerProb = -0.1, upperProb = 0.9), "`lowerProb` must be a single number in \\[0, 1\\].")
})

test_that("zScoreOutliers handles zero variance", {
  df <- data.frame(value = c(10), frequency = c(5))
  expect_error(zScoreOutliers(df), "Standard deviation is 0; z-scores are undefined.")
})

test_that("modifiedZScoreOutliers and generalizedESDOutliers return expected fields", {
  df <- data.frame(
    value = c(5, 6, 7, 8, 30),
    frequency = c(12, 10, 10, 8, 1)
  )

  mz <- modifiedZScoreOutliers(df)
  esd <- generalizedESDOutliers(df, maxOutliers = 5, alpha = 0.05)

  expect_true(all(c("modifiedZScore", "isOutlier") %in% names(mz)))
  expect_true(all(c("outlierCount", "isOutlier") %in% names(esd)))
})

test_that("input validator rejects malformed data", {
  badDf <- data.frame(value = c(1, 2, NA), frequency = c(1, 2, 3))
  expect_error(summarizeValueFrequency(badDf), "`value` and `frequency` columns cannot contain NA.")

  badFreq <- data.frame(value = c(1, 2, 3), frequency = c(1, -1, 2))
  expect_error(summarizeValueFrequency(badFreq), "`frequency` values must be >= 0.")
})

test_that("mcdOutliersRrcovHD flags known multivariate outlier from value-frequency data", {
  skip_if_not_installed("rrcovHD")

  path <- system.file("extdata", "mcd_test_data.csv", package = "AnomalyDetectionMethods")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata", "mcd_test_data.csv")
  }
  expect_true(file.exists(path))

  df <- read.csv(path, stringsAsFactors = FALSE)
  out <- mcdOutliersRrcovHD(
    df = df,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("featureA", "featureB", "featureC"),
    control = "mcd"
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c("robustDistance", "distanceCutoff", "outlierProportion", "isOutlier") %in% names(out)))
  expect_true(is.logical(out$isOutlier))
  expect_true(out$isOutlier[df$value == 40])
  expect_gt(out$outlierProportion[df$value == 40], 0)
})

test_that("mcdOutliersRrcovHD errors for default value-frequency data without feature columns", {
  skip_if_not_installed("rrcovHD")

  df <- data.frame(
    value = c(5, 6, 7, 8, 30),
    frequency = c(12, 10, 10, 8, 1)
  )

  expect_error(
    mcdOutliersRrcovHD(df),
    "No feature columns available. Provide `featureColumns` explicitly."
  )
})

test_that("studentizedResidualOutliers flags known regression outlier", {
  df <- data.frame(
    value = c(11, 14, 17, 20, 23, 40),
    frequency = c(6, 7, 8, 8, 7, 2),
    featureA = c(1, 2, 3, 4, 5, 6),
    featureB = c(1.4, 2.1, 3.9, 3.8, 5.7, 7.5)
  )

  out <- studentizedResidualOutliers(
    df = df,
    valueColumn = "value",
    frequencyColumn = "frequency",
    responseColumn = "value",
    featureColumns = c("featureA", "featureB"),
    residualCutoff = 1.5,
    deleted = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c("studentizedResidual", "residualCutoff", "outlierProportion", "isOutlier") %in% names(out)))
  expect_true(is.logical(out$isOutlier))
  expect_true(out$isOutlier[df$value == 40])
  expect_gt(out$outlierProportion[df$value == 40], 0)
})
