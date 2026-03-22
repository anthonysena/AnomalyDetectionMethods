#' Transform feature columns for Gaussian anomaly detection
#'
#' Applies frequency-aware preprocessing without expanding rows by `frequency`.
#' This is useful before Gaussian anomaly detection when features are on
#' different scales or deviate from approximate normality. The available
#' methods are:
#'
#' - `weighted_zscore`: center by the frequency-weighted mean and divide by the
#'   frequency-weighted standard deviation.
#' - `weighted_mad`: center by the frequency-weighted median and divide by the
#'   frequency-weighted MAD.
#' - `weighted_iqr`: center by the frequency-weighted median and divide by the
#'   frequency-weighted IQR.
#' - `log`: apply the natural logarithm to strictly positive values.
#' - `log1p`: apply `log(1 + x)` to values greater than `-1`.
#' - `signed_log1p`: apply `sign(x) * log(1 + abs(x))`.
#'
#' Assumptions: feature columns are numeric, frequencies are non-negative
#' integer counts interpreted as duplicated observations, and each selected
#' feature satisfies the domain requirements of the chosen transformation.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param featureColumns Character vector of feature column names. If `NULL`,
#'   all columns except `valueColumn` and `frequencyColumn` are used.
#' @param method Transformation method. One of `"weighted_zscore"`,
#'   `"weighted_mad"`, `"weighted_iqr"`, `"log"`, `"log1p"`, or
#'   `"signed_log1p"`.
#' @param interpolate Logical; whether to linearly interpolate quantiles when
#'   computing weighted medians, MAD, or IQR.
#'
#' @return A data.frame with the selected `featureColumns` replaced by their
#'   transformed values.
#' @export
transformGadFeatures <- function(df,
                                 valueColumn = "value",
                                 frequencyColumn = "frequency",
                                 featureColumns = NULL,
                                 method = c(
                                   "weighted_zscore",
                                   "weighted_mad",
                                   "weighted_iqr",
                                   "log",
                                   "log1p",
                                   "signed_log1p"
                                 ),
                                 interpolate = FALSE) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  if (!is.logical(interpolate) || length(interpolate) != 1) {
    stop("`interpolate` must be TRUE or FALSE.")
  }

  method <- match.arg(method)
  if (is.null(featureColumns)) {
    featureColumns <- setdiff(names(df), c(valueColumn, frequencyColumn))
  }
  if (length(featureColumns) < 1) {
    stop("No feature columns available. Provide `featureColumns` explicitly.")
  }
  if (!all(featureColumns %in% names(df))) {
    missing <- setdiff(featureColumns, names(df))
    stop(paste0("Missing feature columns: ", paste(missing, collapse = ", "), "."))
  }

  featuresDf <- df[, featureColumns, drop = FALSE]
  if (!all(vapply(featuresDf, is.numeric, logical(1)))) {
    stop("All `featureColumns` must be numeric.")
  }

  frequencies <- as.numeric(df[[frequencyColumn]])
  positiveIdx <- which(frequencies > 0)
  if (anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  out <- df
  for (featureName in featureColumns) {
    out[[featureName]] <- .transformGadFeature(
      values = out[[featureName]],
      frequencies = frequencies,
      positiveIdx = positiveIdx,
      method = method,
      interpolate = interpolate,
      featureName = featureName
    )
  }

  out
}

#' Identify outliers using Gaussian anomaly detection
#'
#' Fits a frequency-weighted multivariate Gaussian model over the selected
#' feature columns, without expanding rows by `frequency`. The method estimates
#' a weighted mean vector and covariance matrix, computes Mahalanobis distance
#' and Gaussian density for each positive-frequency row, and flags rows whose
#' squared Mahalanobis distance exceeds a chi-square cutoff.
#'
#' Optional preprocessing can be applied through `transformMethod` before
#' fitting the Gaussian model. This is often useful when features are strongly
#' skewed or on very different scales. Robust covariance estimation can also be
#' enabled as a sensitivity analysis through `covarianceMethod = "mcd"`. The
#' robust MCD path estimates center and covariance on an expanded version of
#' the positive-frequency rows because the underlying estimator does not accept
#' frequency weights directly.
#'
#' Assumptions: feature columns are numeric, frequencies are non-negative
#' integer counts interpreted as case weights, and the inlier distribution is
#' approximately Gaussian after any chosen preprocessing step.
#'
#' @param df A data.frame containing value/frequency columns and feature columns.
#' @param valueColumn Name of the column in `df` that holds numeric values.
#' @param frequencyColumn Name of the column in `df` that holds numeric
#'   frequency counts.
#' @param featureColumns Character vector of feature column names. If `NULL`,
#'   all columns except `valueColumn` and `frequencyColumn` are used.
#' @param transformMethod Optional preprocessing method passed to
#'   [transformGadFeatures()]. If `NULL`, no preprocessing is applied.
#' @param covarianceMethod Covariance estimation method. One of `"classical"`
#'   for the standard frequency-weighted covariance or `"mcd"` for a robust
#'   minimum covariance determinant estimate.
#' @param mcdAlpha Optional trimming parameter passed to
#'   [robustbase::covMcd()] when `covarianceMethod = "mcd"`. Must be in
#'   `[0.5, 1]`. Larger values use a larger subset for the robust fit.
#' @param interpolate Logical; whether to linearly interpolate quantiles when
#'   preprocessing uses weighted medians, MAD, or IQR.
#' @param tailProb Upper-tail probability used for chi-square cutoff in `(0, 1)`.
#' @param regularization Non-negative diagonal ridge term added to the weighted
#'   covariance matrix before inversion and density calculation.
#'
#' @return A data.frame with `gaussianLogDensity`, `gaussianDensity`,
#'   `mahalanobisDistance`, `distanceCutoff`, `outlierProportion`, and
#'   `isOutlier` columns appended.
#' @export
gaussianAnomalyOutliers <- function(df,
                                    valueColumn = "value",
                                    frequencyColumn = "frequency",
                                    featureColumns = NULL,
                                    transformMethod = NULL,
                                    covarianceMethod = c("classical", "mcd"),
                                    mcdAlpha = 0.75,
                                    interpolate = FALSE,
                                    tailProb = 0.99,
                                    regularization = 1e-8) {
  .validateValueFrequencyDf(
    df = df,
    valueColumn = valueColumn,
    frequencyColumn = frequencyColumn
  )
  if (!is.null(featureColumns) && (!is.character(featureColumns) || length(featureColumns) < 1)) {
    stop("`featureColumns` must be NULL or a non-empty character vector.")
  }
  covarianceMethod <- match.arg(covarianceMethod)
  validTransformMethods <- c(
    "weighted_zscore",
    "weighted_mad",
    "weighted_iqr",
    "log",
    "log1p",
    "signed_log1p"
  )
  if (!is.null(transformMethod) && (!is.character(transformMethod) || length(transformMethod) != 1 || !(transformMethod %in% validTransformMethods))) {
    stop(paste0(
      "`transformMethod` must be NULL or one of: ",
      paste(validTransformMethods, collapse = ", "),
      "."
    ))
  }
  if (!is.numeric(mcdAlpha) || length(mcdAlpha) != 1 || mcdAlpha < 0.5 || mcdAlpha > 1) {
    stop("`mcdAlpha` must be a single number in [0.5, 1].")
  }
  if (!is.logical(interpolate) || length(interpolate) != 1) {
    stop("`interpolate` must be TRUE or FALSE.")
  }
  if (!is.numeric(tailProb) || length(tailProb) != 1 || tailProb <= 0 || tailProb >= 1) {
    stop("`tailProb` must be a single number in (0, 1).")
  }
  if (!is.numeric(regularization) || length(regularization) != 1 || regularization < 0) {
    stop("`regularization` must be a single non-negative number.")
  }

  if (is.null(featureColumns)) {
    featureColumns <- setdiff(names(df), c(valueColumn, frequencyColumn))
  }
  if (length(featureColumns) < 1) {
    stop("No feature columns available. Provide `featureColumns` explicitly.")
  }
  if (!all(featureColumns %in% names(df))) {
    missing <- setdiff(featureColumns, names(df))
    stop(paste0("Missing feature columns: ", paste(missing, collapse = ", "), "."))
  }

  workingDf <- if (is.null(transformMethod)) {
    df
  } else {
    transformGadFeatures(
      df = df,
      valueColumn = valueColumn,
      frequencyColumn = frequencyColumn,
      featureColumns = featureColumns,
      method = transformMethod,
      interpolate = interpolate
    )
  }

  featuresDf <- workingDf[, featureColumns, drop = FALSE]
  if (!all(vapply(featuresDf, is.numeric, logical(1)))) {
    stop("All `featureColumns` must be numeric.")
  }

  frequencies <- as.numeric(df[[frequencyColumn]])
  positiveIdx <- which(frequencies > 0)
  if (length(positiveIdx) < 2) {
    stop("Need at least 2 rows with positive frequency for Gaussian fitting.")
  }
  if (anyNA(featuresDf[positiveIdx, , drop = FALSE])) {
    stop("Feature columns cannot contain NA for rows with positive frequency.")
  }

  x <- as.matrix(featuresDf[positiveIdx, , drop = FALSE])
  w <- frequencies[positiveIdx]
  totalWeight <- sum(w)
  if (totalWeight <= 1) {
    stop("Total positive frequency must be > 1 for covariance estimation.")
  }

  covarianceFit <- .fitGadCovariance(
    x = x,
    frequencies = w,
    covarianceMethod = covarianceMethod,
    mcdAlpha = mcdAlpha
  )
  center <- covarianceFit$center
  covariance <- covarianceFit$covariance
  xCentered <- sweep(x, 2, center, "-")

  p <- ncol(x)
  diagScale <- max(1, mean(diag(covariance)))
  covarianceReg <- covariance + diag(regularization * diagScale, p)

  covarianceInv <- tryCatch(
    solve(covarianceReg),
    error = function(e) {
      sv <- svd(covarianceReg)
      tol <- max(dim(covarianceReg)) * max(sv$d) * .Machine$double.eps
      dInv <- ifelse(sv$d > tol, 1 / sv$d, 0)
      sv$v %*% diag(dInv, length(dInv)) %*% t(sv$u)
    }
  )

  determinantInfo <- determinant(covarianceReg, logarithm = TRUE)
  if (determinantInfo$sign <= 0 || !is.finite(as.numeric(determinantInfo$modulus))) {
    stop("Covariance matrix is singular; Gaussian density is undefined. Increase `regularization`.")
  }

  distanceSqPositive <- rowSums((xCentered %*% covarianceInv) * xCentered)
  distanceSqPositive <- pmax(distanceSqPositive, 0)
  logDet <- as.numeric(determinantInfo$modulus)
  gaussianLogDensityPositive <- -0.5 * (p * log(2 * pi) + logDet + distanceSqPositive)
  gaussianDensityPositive <- exp(gaussianLogDensityPositive)
  cutoffSq <- stats::qchisq(tailProb, df = p)

  gaussianLogDensity <- rep(NA_real_, nrow(df))
  gaussianDensity <- rep(NA_real_, nrow(df))
  mahalanobisDistance <- rep(NA_real_, nrow(df))
  outlierProportion <- rep(0, nrow(df))

  gaussianLogDensity[positiveIdx] <- gaussianLogDensityPositive
  gaussianDensity[positiveIdx] <- gaussianDensityPositive
  mahalanobisDistance[positiveIdx] <- sqrt(distanceSqPositive)
  outlierProportion[positiveIdx] <- as.numeric(distanceSqPositive > cutoffSq)

  out <- df
  out$gaussianLogDensity <- gaussianLogDensity
  out$gaussianDensity <- gaussianDensity
  out$mahalanobisDistance <- mahalanobisDistance
  out$distanceCutoff <- sqrt(cutoffSq)
  out$outlierProportion <- outlierProportion
  out$isOutlier <- outlierProportion > 0

  out
}

.transformGadFeature <- function(values,
                                 frequencies,
                                 positiveIdx,
                                 method,
                                 interpolate,
                                 featureName) {
  if (method %in% c("weighted_zscore", "weighted_mad", "weighted_iqr")) {
    return(.scaleWeightedFeature(
      values = values,
      frequencies = frequencies,
      positiveIdx = positiveIdx,
      method = method,
      interpolate = interpolate
    ))
  }

  transformedValues <- rep(NA_real_, length(values))
  positiveValues <- values[positiveIdx]

  if (method == "log") {
    if (any(positiveValues <= 0)) {
      stop(paste0(
        "Feature `", featureName, "` contains non-positive values; `log` transformation is undefined."
      ))
    }
    transformedValues[positiveIdx] <- log(positiveValues)
    return(transformedValues)
  }

  if (method == "log1p") {
    if (any(positiveValues <= -1)) {
      stop(paste0(
        "Feature `", featureName, "` contains values <= -1; `log1p` transformation is undefined."
      ))
    }
    transformedValues[positiveIdx] <- log1p(positiveValues)
    return(transformedValues)
  }

  transformedValues[positiveIdx] <- sign(positiveValues) * log1p(abs(positiveValues))
  transformedValues
}

.fitGadCovariance <- function(x, frequencies, covarianceMethod, mcdAlpha) {
  if (covarianceMethod == "classical") {
    totalWeight <- sum(frequencies)
    center <- colSums(x * frequencies) / totalWeight
    xCentered <- sweep(x, 2, center, "-")
    weightedXCentered <- xCentered * sqrt(frequencies)
    covariance <- crossprod(weightedXCentered) / (totalWeight - 1)
    return(list(center = center, covariance = covariance))
  }

  if (!requireNamespace("robustbase", quietly = TRUE)) {
    stop("Package `robustbase` is required for `covarianceMethod = \"mcd\"`. Install it first.")
  }

  expandedIdx <- rep(seq_len(nrow(x)), frequencies)
  if (length(expandedIdx) < ncol(x) + 1) {
    stop("Need more expanded observations than feature columns for robust covariance estimation.")
  }
  xExpanded <- x[expandedIdx, , drop = FALSE]
  mcdFit <- robustbase::covMcd(xExpanded, alpha = mcdAlpha)
  list(center = as.numeric(mcdFit$center), covariance = as.matrix(mcdFit$cov))
}
