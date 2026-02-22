library(AnomalyDetectionMethods)

# ---- Baseline example data (aggregated value-frequency) ----
baselineDf <- data.frame(
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

baselineSummary <- summarizeValueFrequency(baselineDf)
baselineSummary

baselineWithFences <- tukeyFences(baselineDf)
baselineWithFences[baselineWithFences$isOutlier, ]

baselineWithQuantileThresholds <- quantileThresholds(baselineDf)
baselineWithQuantileThresholds[baselineWithQuantileThresholds$isOutlier, ]

baselineWithZscoreOutliers <- zScoreOutliers(baselineDf)
baselineWithZscoreOutliers[baselineWithZscoreOutliers$isOutlier, ]

baselineWithModZscoreOutliers <- modifiedZScoreOutliers(baselineDf)
baselineWithModZscoreOutliers[baselineWithModZscoreOutliers$isOutlier, ]

baselineWithGeneralizedEsd <- generalizedESDOutliers(baselineDf)
baselineWithGeneralizedEsd[baselineWithGeneralizedEsd$isOutlier, ]

plotWeightedQQ(baselineDf)
plotWeightedHistogram(baselineDf)
plotWeightedDensity(baselineDf)
plotWeightedECDF(baselineDf)

# ---- Example data with 5 intentional outliers (aggregated value-frequency) ----
outlierDf <- data.frame(
  value = c(
    35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
    85, 90, 95, 100, 8, 12, 140, 145, 150
  ),
  frequency = c(
    1, 2, 4, 8, 12, 18, 22, 24, 22, 18,
    12, 8, 4, 2, 1, 1, 1, 1, 1
  )
)

# The 5 rare extreme values we expect to be outliers in this distribution.
expectedOutlierValues <- c(8, 12, 140, 145, 150)

outlierSummary <- summarizeValueFrequency(outlierDf)
outlierSummary

outlierWithFences <- tukeyFences(outlierDf)
outlierWithFences[outlierWithFences$isOutlier, c("value", "frequency")]

outlierWithQuantileThresholds <- quantileThresholds(
  outlierDf,
  lowerProb = 0.015,
  upperProb = 0.985,
  interpolate = TRUE
)
outlierWithQuantileThresholds[
  outlierWithQuantileThresholds$isOutlier,
  c("value", "frequency")
]

outlierWithZscore <- zScoreOutliers(outlierDf, zCutoff = 2.5)
outlierWithZscore[
  outlierWithZscore$isOutlier,
  c("value", "frequency", "zScore")
]

outlierWithModZscore <- modifiedZScoreOutliers(outlierDf, zCutoff = 3.5)
outlierWithModZscore[
  outlierWithModZscore$isOutlier,
  c("value", "frequency", "modifiedZScore")
]

# Generalized ESD can be more conservative on weighted/aggregated data.
outlierWithGeneralizedEsd <- generalizedESDOutliers(
  outlierDf,
  maxOutliers = 5,
  alpha = 0.05
)
outlierWithGeneralizedEsd[
  outlierWithGeneralizedEsd$isOutlier,
  c("value", "frequency", "outlierCount")
]

# Quick checks for methods expected to return the 5 designed outliers.
identical(
  sort(outlierWithFences$value[outlierWithFences$isOutlier]),
  expectedOutlierValues
)
identical(
  sort(
    outlierWithQuantileThresholds$value[
      outlierWithQuantileThresholds$isOutlier
    ]
  ),
  expectedOutlierValues
)
identical(
  sort(outlierWithZscore$value[outlierWithZscore$isOutlier]),
  expectedOutlierValues
)
identical(
  sort(outlierWithModZscore$value[outlierWithModZscore$isOutlier]),
  expectedOutlierValues
)

plotWeightedQQ(outlierDf)
plotWeightedHistogram(outlierDf)
plotWeightedDensity(outlierDf)
plotWeightedECDF(outlierDf)
