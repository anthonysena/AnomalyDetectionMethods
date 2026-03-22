test_that("transformGadFeatures performs weighted z-score scaling", {
  df <- data.frame(
    value = c(1, 2, 3, 4),
    frequency = c(1, 2, 1, 0),
    featureA = c(10, 20, 30, 40),
    featureB = c(1, 3, 5, 7)
  )

  out <- transformGadFeatures(
    df = df,
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

test_that("transformGadFeatures performs logarithmic transformations", {
  df <- data.frame(
    value = c(1, 2, 3, 4),
    frequency = c(1, 1, 1, 0),
    featureA = c(1, 3, 8, 10),
    featureB = c(-2, 0, 3, 5)
  )

  logOut <- transformGadFeatures(df, featureColumns = "featureA", method = "log")
  log1pOut <- transformGadFeatures(df, featureColumns = "featureA", method = "log1p")
  signedOut <- transformGadFeatures(df, featureColumns = "featureB", method = "signed_log1p")

  expect_equal(logOut$featureA, c(log(1), log(3), log(8), NA_real_))
  expect_equal(log1pOut$featureA, c(log1p(1), log1p(3), log1p(8), NA_real_))
  expect_equal(
    signedOut$featureB,
    c(-log1p(2), 0, log1p(3), NA_real_)
  )
})

test_that("transformGadFeatures validates transformation domains", {
  badLogDf <- data.frame(
    value = c(1, 2),
    frequency = c(1, 1),
    featureA = c(0, 2)
  )
  badLog1pDf <- data.frame(
    value = c(1, 2),
    frequency = c(1, 1),
    featureA = c(-1, 2)
  )

  expect_error(
    transformGadFeatures(badLogDf, featureColumns = "featureA", method = "log"),
    "Feature `featureA` contains non-positive values; `log` transformation is undefined."
  )
  expect_error(
    transformGadFeatures(badLog1pDf, featureColumns = "featureA", method = "log1p"),
    "Feature `featureA` contains values <= -1; `log1p` transformation is undefined."
  )
})

test_that("gaussianAnomalyOutliers flags an isolated multivariate point", {
  df <- data.frame(
    value = c(10, 11, 12, 13, 100),
    frequency = c(10, 8, 8, 8, 1),
    featureA = c(0.00, 0.10, -0.10, 0.05, 4.00),
    featureB = c(0.00, 0.05, -0.05, -0.10, 4.20)
  )

  out <- gaussianAnomalyOutliers(
    df = df,
    featureColumns = c("featureA", "featureB"),
    tailProb = 0.99
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c(
    "gaussianLogDensity",
    "gaussianDensity",
    "mahalanobisDistance",
    "distanceCutoff",
    "outlierProportion",
    "isOutlier"
  ) %in% names(out)))
  expect_true(out$isOutlier[df$value == 100])
  expect_gt(out$mahalanobisDistance[df$value == 100], out$distanceCutoff[df$value == 100])
})

test_that("gaussianAnomalyOutliers supports preprocessing transforms", {
  df <- data.frame(
    value = c(10, 11, 12, 13, 100),
    frequency = c(10, 8, 8, 8, 1),
    featureA = c(1, 2, 3, 4, 100),
    featureB = c(1, 2, 4, 8, 256)
  )

  out <- gaussianAnomalyOutliers(
    df = df,
    featureColumns = c("featureA", "featureB"),
    transformMethod = "log1p",
    tailProb = 0.99
  )

  expect_true(out$isOutlier[df$value == 100])
  expect_true(all(is.finite(out$gaussianLogDensity[df$frequency > 0])))
  expect_true(all(out$gaussianDensity[df$frequency > 0] > 0))
})

test_that("gaussianAnomalyOutliers supports robust covariance sensitivity analysis", {
  skip_if_not_installed("robustbase")

  df <- data.frame(
    value = 1:11,
    frequency = c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1),
    featureA = c(-1.2, -0.8, -0.5, -0.1, 0.0, 0.2, 0.5, 0.8, 1.0, 1.3, 5.0),
    featureB = c(-0.9, -0.6, -0.4, -0.2, 0.1, 0.2, 0.4, 0.7, 0.9, 1.1, 4.7)
  )

  classicalOut <- gaussianAnomalyOutliers(
    df = df,
    featureColumns = c("featureA", "featureB"),
    covarianceMethod = "classical",
    tailProb = 0.99
  )
  robustOut <- gaussianAnomalyOutliers(
    df = df,
    featureColumns = c("featureA", "featureB"),
    covarianceMethod = "mcd",
    mcdAlpha = 0.9,
    tailProb = 0.99
  )

  expect_true(robustOut$isOutlier[df$value == 11])
  expect_gt(
    robustOut$mahalanobisDistance[df$value == 11],
    classicalOut$mahalanobisDistance[df$value == 11]
  )
  expect_false(isTRUE(all.equal(
    classicalOut$mahalanobisDistance,
    robustOut$mahalanobisDistance
  )))
})
