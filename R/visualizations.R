#' Plot a normal Q-Q plot for value-frequency data
#'
#' Creates a normal Q-Q plot using frequency-weighted plotting positions. Each
#' unique value is plotted once using its cumulative frequency to determine the
#' corresponding theoretical normal quantile.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param addLine Logical; whether to add a reference line based on the
#'   weighted mean and standard deviation.
#'
#' @return A ggplot object with theoretical and observed quantiles.
#' @export
plotWeightedQQ <- function(df,
                           valueColumn = "value",
                           frequencyColumn = "frequency",
                           addLine = TRUE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  totalN <- sum(frequencies)

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]
  cumulativeFrequencies <- cumsum(frequencies)
  probs <- (cumulativeFrequencies - 0.5 * frequencies) / totalN

  stats <- .weightedMeanSd(values, frequencies)
  theoretical <- stats::qnorm(probs, mean = stats$mean, sd = stats$sd)

  plotData <- data.frame(
    theoretical = theoretical,
    observed = values,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$theoretical, y = .data$observed)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Theoretical Quantiles", y = "Observed Values")

  if (addLine) {
    p <- p + ggplot2::geom_abline(intercept = stats$mean, slope = stats$sd, linetype = 2, color = "red")
  }

  p
}

#' Plot a weighted histogram with normal overlay
#'
#' Plots a histogram based on value-frequency data using weighted bin counts.
#' Optionally overlays the fitted normal curve using the weighted mean and
#' standard deviation.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param breaks Histogram breaks passed to `hist()` to define bin edges.
#' @param addNormal Logical; whether to overlay the fitted normal curve.
#'
#' @return A ggplot object with weighted histogram and optional normal curve.
#' @export
plotWeightedHistogram <- function(df,
                                  valueColumn = "value",
                                  frequencyColumn = "frequency",
                                  breaks = "Sturges",
                                  addNormal = TRUE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  totalN <- sum(frequencies)

  baseHist <- graphics::hist(values, breaks = breaks, plot = FALSE)
  bins <- cut(values, breaks = baseHist$breaks, include.lowest = TRUE, right = TRUE)
  counts <- tapply(frequencies, bins, sum)
  counts <- as.numeric(counts)
  counts[is.na(counts)] <- 0
  mids <- baseHist$mids
  binWidth <- diff(baseHist$breaks)[1]

  plotData <- data.frame(
    mids = mids,
    counts = counts,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$mids, y = .data$counts)
  ) +
    ggplot2::geom_col(width = binWidth, fill = "grey70", color = "grey40") +
    ggplot2::labs(x = "Value", y = "Weighted Count")

  if (addNormal) {
    stats <- .weightedMeanSd(values, frequencies)
    xGrid <- seq(min(baseHist$breaks), max(baseHist$breaks), length.out = 512)
    yNormal <- stats::dnorm(xGrid, mean = stats$mean, sd = stats$sd) * totalN * binWidth
    normalData <- data.frame(x = xGrid, y = yNormal, stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_line(
      data = normalData,
      ggplot2::aes(x = .data$x, y = .data$y),
      linetype = 2,
      color = "red"
    )
  }

  p
}

#' Plot a weighted density with normal overlay
#'
#' Computes a kernel density estimate for value-frequency data. By default, the
#' data are expanded by frequency up to `maxExpand` to use `density()`; if the
#' total count exceeds `maxExpand`, a weighted approximation is used.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param n Number of points used to estimate the density.
#' @param maxExpand Maximum total frequency to expand directly for `density()`.
#' @param addNormal Logical; whether to overlay the fitted normal density.
#'
#' @return A ggplot object with weighted density and optional normal curve.
#' @export
plotWeightedDensity <- function(df,
                                valueColumn = "value",
                                frequencyColumn = "frequency",
                                n = 512,
                                maxExpand = 50000,
                                addNormal = TRUE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  totalN <- sum(frequencies)

  if (totalN <= maxExpand) {
    expanded <- rep(values, frequencies)
    dens <- stats::density(expanded, n = n)
  } else {
    bw <- stats::bw.nrd0(values)
    xGrid <- seq(min(values), max(values), length.out = n)
    weights <- frequencies / totalN
    kernel <- exp(-0.5 * ((outer(xGrid, values, "-") / bw) ^ 2)) / (bw * sqrt(2 * pi))
    yGrid <- kernel %*% weights
    dens <- list(x = xGrid, y = as.numeric(yGrid))
  }

  plotData <- data.frame(x = dens$x, y = dens$y, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Value", y = "Density")

  if (addNormal) {
    stats <- .weightedMeanSd(values, frequencies)
    yNormal <- stats::dnorm(dens$x, mean = stats$mean, sd = stats$sd)
    normalData <- data.frame(x = dens$x, y = yNormal, stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_line(
      data = normalData,
      ggplot2::aes(x = .data$x, y = .data$y),
      linetype = 2,
      color = "red"
    )
  }

  p
}

#' Plot a weighted ECDF with normal CDF overlay
#'
#' Plots the empirical CDF derived from value-frequency data. Optionally
#' overlays the fitted normal CDF using the weighted mean and standard
#' deviation.
#'
#' @param df A data.frame containing a value column and a frequency column.
#' @param valueColumn Name of the column in `df` that holds the numeric values.
#' @param frequencyColumn Name of the column in `df` that holds the numeric
#'   frequency counts.
#' @param addNormal Logical; whether to overlay the fitted normal CDF.
#'
#' @return A ggplot object with weighted ECDF and optional normal CDF.
#' @export
plotWeightedECDF <- function(df,
                             valueColumn = "value",
                             frequencyColumn = "frequency",
                             addNormal = TRUE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )

  values <- df[[valueColumn]]
  frequencies <- df[[frequencyColumn]]
  totalN <- sum(frequencies)

  ordering <- order(values)
  values <- values[ordering]
  frequencies <- frequencies[ordering]
  ecdfValues <- cumsum(frequencies) / totalN

  plotData <- data.frame(values = values, ecdf = ecdfValues, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data$values, y = .data$ecdf)
  ) +
    ggplot2::geom_step() +
    ggplot2::labs(x = "Value", y = "ECDF")

  if (addNormal) {
    stats <- .weightedMeanSd(values, frequencies)
    normalData <- data.frame(values = values, cdf = stats::pnorm(values, mean = stats$mean, sd = stats$sd), stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_line(
      data = normalData,
      ggplot2::aes(x = .data$values, y = .data$cdf),
      linetype = 2,
      color = "red"
    )
  }

  p
}
