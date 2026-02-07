source("CrossSectionalMethods.R")

# ---- Example data (same 40 observations as before, aggregated) ----
df <- data.frame(
  value = c(
    12, 15, 18, 20, 22, 25, 27, 28, 30, 32,
    35, 38, 40, 42, 45, 48, 50, 52, 55, 58,
    60, 62, 65, 68, 70, 75
  ),
  frequency = c(
    1, 2, 1, 3, 1, 2, 1, 1, 4, 1,
    2, 1, 2, 1, 2, 1, 3, 1, 2, 1,
    2, 1, 1, 1, 1, 1
  )
)

valueFrequencySummary <- summarizeValueFrequency(df)

# TukeyFences -------------
dfWithFences <- tukeyFences(df)
dfWithFences[dfWithFences$isOutlier, ]

# QuantileThreshold -------------
dfWithQuantileThresholds <- quantileThresholds(df)
dfWithQuantileThresholds[dfWithQuantileThresholds$isOutlier,]

# zScore ---------------
dfWithZscoreOutliers <- zScoreOutliers(df)
dfWithZscoreOutliers[dfWithZscoreOutliers$isOutlier, ]