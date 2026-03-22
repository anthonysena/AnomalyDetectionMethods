test_that("isolationForestOutliers flags an isolated point", {
  skip_if_not_installed("isotree")

  df <- data.frame(
    value = c(10, 11, 12, 100),
    frequency = c(3, 2, 2, 1),
    featureA = c(0.0, 0.1, -0.1, 5.0),
    featureB = c(0.0, -0.1, 0.1, 5.2)
  )

  out <- isolationForestOutliers(
    df = df,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("featureA", "featureB"),
    ntrees = 200,
    sampleSize = 4,
    scoreCutoff = 0.6,
    weightsAsSampleProb = FALSE,
    standardizeData = FALSE,
    seed = 123,
    nthreads = 1
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c("isolationForestScore", "scoreCutoff", "outlierProportion", "isOutlier") %in% names(out)))
  expect_true(is.logical(out$isOutlier))
  expect_true(out$isOutlier[df$value == 100])
  expect_false(any(out$isOutlier[df$value != 100]))
})

test_that("isolationForestOutliers tracks expanded-data ranking under multiplicity weights", {
  skip_if_not_installed("isotree")

  weightedDf <- data.frame(
    value = c(10, 11, 12, 100),
    frequency = c(3, 2, 2, 1),
    featureA = c(0.0, 0.1, -0.1, 5.0),
    featureB = c(0.0, -0.1, 0.1, 5.2)
  )

  expandedIdx <- rep(seq_len(nrow(weightedDf)), weightedDf$frequency)
  expandedDf <- weightedDf[expandedIdx, c("value", "featureA", "featureB")]
  expandedDf$frequency <- 1

  weightedOut <- isolationForestOutliers(
    df = weightedDf,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("featureA", "featureB"),
    ntrees = 2000,
    sampleSize = 4,
    weightsAsSampleProb = FALSE,
    standardizeData = FALSE,
    seed = 123,
    nthreads = 1
  )
  expandedOut <- isolationForestOutliers(
    df = expandedDf,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("featureA", "featureB"),
    ntrees = 2000,
    sampleSize = 4,
    weightsAsSampleProb = FALSE,
    standardizeData = FALSE,
    seed = 123,
    nthreads = 1
  )

  aggregatedExpandedScores <- tapply(expandedOut$isolationForestScore, expandedIdx, mean)

  expect_equal(order(weightedOut$isolationForestScore), order(aggregatedExpandedScores))
  expect_equal(unname(which.max(weightedOut$isolationForestScore)), unname(which.max(aggregatedExpandedScores)))
  expect_gt(
    stats::cor(weightedOut$isolationForestScore, as.numeric(aggregatedExpandedScores)),
    0.95
  )
})
