test_that("scaleLofFeatures performs weighted z-score scaling", {
  df <- data.frame(
    value = c(1, 2, 3, 4),
    frequency = c(1, 2, 1, 0),
    featureA = c(10, 20, 30, 40),
    featureB = c(1, 3, 5, 7)
  )

  out <- scaleLofFeatures(
    df = df,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("featureA", "featureB"),
    method = "weighted_zscore"
  )

  positiveIdx <- which(df$frequency > 0)
  weightedMeanA <- sum(out$featureA[positiveIdx] * df$frequency[positiveIdx]) / sum(df$frequency[positiveIdx])
  weightedMeanB <- sum(out$featureB[positiveIdx] * df$frequency[positiveIdx]) / sum(df$frequency[positiveIdx])
  weightedSdA <- sqrt(sum(df$frequency[positiveIdx] * out$featureA[positiveIdx]^2) / sum(df$frequency[positiveIdx]))
  weightedSdB <- sqrt(sum(df$frequency[positiveIdx] * out$featureB[positiveIdx]^2) / sum(df$frequency[positiveIdx]))

  expect_equal(weightedMeanA, 0, tolerance = 1e-10)
  expect_equal(weightedMeanB, 0, tolerance = 1e-10)
  expect_equal(weightedSdA, 1, tolerance = 1e-10)
  expect_equal(weightedSdB, 1, tolerance = 1e-10)
  expect_true(is.na(out$featureA[df$frequency == 0]))
})

test_that("scaleLofFeatures performs weighted MAD and IQR scaling", {
  df <- data.frame(
    value = c(1, 2, 3, 4, 5),
    frequency = c(1, 1, 1, 1, 0),
    featureA = c(1, 2, 3, 4, 5)
  )

  madOut <- scaleLofFeatures(
    df = df,
    featureColumns = "featureA",
    method = "weighted_mad"
  )
  iqrOut <- scaleLofFeatures(
    df = df,
    featureColumns = "featureA",
    method = "weighted_iqr"
  )

  expect_equal(madOut$featureA, c(-1, 0, 1, 2, NA))
  expect_equal(iqrOut$featureA, c(-0.5, 0, 0.5, 1, NA))
})

test_that("scaleLofFeatures rejects zero-spread features under the chosen method", {
  df <- data.frame(
    value = c(1, 2, 3),
    frequency = c(1, 1, 1),
    featureA = c(5, 5, 5)
  )

  expect_error(
    scaleLofFeatures(df, featureColumns = "featureA", method = "weighted_zscore"),
    "Selected feature has 0 weighted standard deviation; weighted z-score scaling is undefined."
  )
  expect_error(
    scaleLofFeatures(df, featureColumns = "featureA", method = "weighted_mad"),
    "Selected feature has 0 weighted MAD; weighted MAD scaling is undefined."
  )
  expect_error(
    scaleLofFeatures(df, featureColumns = "featureA", method = "weighted_iqr"),
    "Selected feature has 0 weighted IQR; weighted IQR scaling is undefined."
  )
})

test_that("localOutlierFactorOutliers flags an isolated multivariate point", {
  df <- data.frame(
    value = c(10, 11, 12, 13, 14, 100),
    frequency = c(2, 2, 2, 2, 2, 1),
    featureA = c(0.00, 0.10, -0.10, 0.05, -0.05, 4.00),
    featureB = c(0.00, 0.05, -0.05, -0.10, 0.10, 4.20)
  )

  out <- localOutlierFactorOutliers(
    df = df,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("featureA", "featureB"),
    k = 4,
    scoreCutoff = 1.2
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c(
    "localOutlierFactor",
    "kDistance",
    "localReachabilityDensity",
    "scoreCutoff",
    "outlierProportion",
    "isOutlier"
  ) %in% names(out)))
  expect_true(is.logical(out$isOutlier))
  expect_true(out$isOutlier[df$value == 100])
  expect_gt(out$localOutlierFactor[df$value == 100], out$scoreCutoff[df$value == 100])
  expect_false(any(out$isOutlier[df$value != 100]))
})

test_that("localOutlierFactorOutliers validates neighborhood size and feature variation", {
  df <- data.frame(
    value = c(10, 11, 12),
    frequency = c(1, 1, 1),
    featureA = c(0, 1, 2),
    featureB = c(0, 1, 2)
  )

  expect_error(
    localOutlierFactorOutliers(df, featureColumns = c("featureA", "featureB"), k = 3),
    "`k` must be smaller than the total positive frequency."
  )

  identicalDf <- data.frame(
    value = c(10, 11),
    frequency = c(2, 1),
    featureA = c(0, 0),
    featureB = c(0, 0)
  )

  expect_error(
    localOutlierFactorOutliers(identicalDf, featureColumns = c("featureA", "featureB"), k = 1),
    "All positive-frequency rows have identical feature values; LOF is undefined."
  )
})
