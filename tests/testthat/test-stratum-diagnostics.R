test_that("summarizeStratumDiagnostics returns expected fields", {
  df <- data.frame(
    value = c(110, 115, 120, 125, 130),
    frequency = c(5, 10, 20, 10, 5),
    ageMid = c(42, 42, 42, 42, 42)
  )

  out <- summarizeStratumDiagnostics(
    df,
    kValues = c(3, 5),
    featureColumns = "ageMid"
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c(
    "distinctValues",
    "totalFrequency",
    "top1Mass",
    "weightedSd",
    "propEndingIn0or5",
    "supportsK_3",
    "supportsK_5",
    "featureCount"
  ) %in% names(out)))
  expect_equal(out$distinctValues, 5)
  expect_equal(out$totalFrequency, 50)
})

test_that("summarizeStrataDiagnostics and flagAnalyzableStrata return expected outputs", {
  df <- data.frame(
    site = c("A", "A", "A", "B", "B", "B"),
    gender = c("F", "F", "F", "M", "M", "M"),
    value = c(110, 120, 130, 115, 120, 125),
    frequency = c(10, 20, 10, 2, 2, 1),
    ageMid = c(40, 40, 40, 65, 65, 65)
  )

  diagOut <- summarizeStrataDiagnostics(
    df,
    strataColumns = c("site", "gender"),
    kValues = c(2, 4),
    featureColumns = "ageMid"
  )
  flagged <- flagAnalyzableStrata(
    diagOut,
    minDistinctValues = 3,
    minTotalFrequency = 10,
    lofK = 2
  )

  expect_equal(nrow(diagOut), 2)
  expect_true(all(c("site", "gender", "distinctValues", "totalFrequency") %in% names(diagOut)))
  expect_true(all(c("usableForLof", "usableForIf", "usableForGad") %in% names(flagged)))
  expect_true(flagged$usableForLof[flagged$site == "A"])
  expect_false(flagged$usableForLof[flagged$site == "B"])
})

test_that("plotStratumDiagnostics returns a ggplot object", {
  diagnosticsDf <- data.frame(
    totalFrequency = c(50, 20, 100),
    distinctValues = c(10, 5, 25),
    usableForLof = c(TRUE, FALSE, TRUE)
  )

  expect_s3_class(
    plotStratumDiagnostics(
      diagnosticsDf,
      xColumn = "totalFrequency",
      yColumn = "distinctValues",
      colorColumn = "usableForLof"
    ),
    "ggplot"
  )
})
