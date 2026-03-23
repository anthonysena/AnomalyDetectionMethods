test_that("evaluateLofGrid returns expected summaries", {
  df <- data.frame(
    value = c(1, 2, 3, 4, 5, 6),
    frequency = c(3, 3, 3, 3, 3, 1),
    featureA = c(0.0, 0.1, -0.1, 4.0, 4.1, 8.0),
    featureB = c(0.0, -0.1, 0.1, 4.2, 3.9, 8.1),
    featureC = c(1.0, 1.1, 0.9, 1.0, 1.2, 2.5)
  )

  out <- evaluateLofGrid(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    kValues = c(2, 3),
    scoreCutoffs = c(1.1, 1.3),
    topN = 3
  )

  expect_true(all(c("summary", "stability") %in% names(out)))
  expect_true(all(c("k", "scoreCutoff", "flagCount", "weightedFlagCount", "topNValues") %in% names(out$summary)))
  expect_true(all(c("configA", "configB", "spearman", "topNJaccard") %in% names(out$stability)))
})

test_that("evaluateIsolationForestGrid returns expected summaries", {
  skip_if_not_installed("isotree")

  df <- data.frame(
    value = c(1, 2, 3, 4, 5, 6),
    frequency = c(3, 3, 3, 3, 3, 1),
    featureA = c(0.0, 0.1, -0.1, 4.0, 4.1, 8.0),
    featureB = c(0.0, -0.1, 0.1, 4.2, 3.9, 8.1),
    featureC = c(1.0, 1.1, 0.9, 1.0, 1.2, 2.5)
  )

  out <- evaluateIsolationForestGrid(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    ntreesValues = c(50, 100),
    sampleSizeValues = c(4, 6),
    scoreCutoffs = c(0.55, 0.65),
    seeds = c(1, 2),
    topN = 3,
    weightsAsSampleProb = FALSE,
    nthreads = 1
  )

  expect_true(all(c("summary", "stability") %in% names(out)))
  expect_true(all(c(
    "ntrees",
    "sampleSize",
    "scoreCutoff",
    "meanFlagCount",
    "meanWeightedFlagCount",
    "meanPairwiseSpearman",
    "meanPairwiseTopJaccard"
  ) %in% names(out$summary)))
  expect_true(all(c("configA", "configB", "spearman", "topNJaccard") %in% names(out$stability)))
})

test_that("evaluateGadGrid returns expected summaries", {
  skip_if_not_installed("robustbase")

  df <- data.frame(
    value = c(1, 2, 3, 4, 5, 6),
    frequency = c(3, 3, 3, 3, 3, 1),
    featureA = c(0.0, 0.1, -0.1, 4.0, 4.1, 8.0),
    featureB = c(0.0, -0.1, 0.1, 4.2, 3.9, 8.1),
    featureC = c(1.0, 1.1, 0.9, 1.0, 1.2, 2.5)
  )

  out <- evaluateGadGrid(
    df,
    featureColumns = c("featureA", "featureB", "featureC"),
    transformMethods = c(NA, "weighted_zscore"),
    covarianceMethods = c("classical", "mcd"),
    tailProbValues = c(0.95, 0.99),
    topN = 3,
    mcdAlpha = 0.9
  )

  expect_true(all(c("summary", "stability") %in% names(out)))
  expect_true(all(c(
    "transformMethod",
    "covarianceMethod",
    "tailProb",
    "flagCount",
    "weightedFlagCount",
    "qqCorrelation"
  ) %in% names(out$summary)))
  expect_true(all(c("configA", "configB", "spearman", "topNJaccard") %in% names(out$stability)))
})
