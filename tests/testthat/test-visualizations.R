test_that("plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    value = c(1, 2, 3, 4, 5),
    frequency = c(2, 3, 4, 3, 2)
  )

  p1 <- plotWeightedQQ(df)
  p2 <- plotWeightedHistogram(df)
  p3 <- plotWeightedDensity(df)
  p4 <- plotWeightedECDF(df)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})
