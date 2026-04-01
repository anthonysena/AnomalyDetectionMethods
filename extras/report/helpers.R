report_check_required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required package(s): ",
        paste(missing, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
}

report_has_andromeda <- function() {
  requireNamespace("Andromeda", quietly = TRUE)
}

report_has_quarto_package <- function() {
  requireNamespace("quarto", quietly = TRUE)
}

report_has_quarto_cli <- function() {
  nzchar(Sys.which("quarto"))
}

report_age_band <- function(age,
                            breaks = c(18, 40, 65, 80, Inf),
                            labels = c("18-39", "40-64", "65-79", "80+")) {
  cut(
    age,
    breaks = breaks,
    labels = labels,
    right = FALSE,
    include.lowest = TRUE
  )
}

report_gender_map <- function(genderConceptIds) {
  genderIds <- unique(as.numeric(genderConceptIds))
  labels <- as.character(genderIds)
  knownLabels <- c(
    `8507` = "Male",
    `8532` = "Female",
    `8551` = "Unknown"
  )
  matched <- names(knownLabels) %in% as.character(genderIds)
  labels[match(names(knownLabels)[matched], as.character(genderIds))] <- knownLabels[matched]

  data.frame(
    gender_concept_id = genderIds,
    gender_label = labels,
    stringsAsFactors = FALSE
  )
}

report_encode_gender <- function(genderConceptIds) {
  uniqueIds <- sort(unique(as.numeric(genderConceptIds)))
  encoded <- seq_along(uniqueIds) - 1
  stats::setNames(encoded, uniqueIds)[as.character(genderConceptIds)]
}

report_repo_root <- function(baseDir = getwd()) {
  normalizePath(baseDir, winslash = "/", mustWork = TRUE)
}

report_paths <- function(baseDir = getwd()) {
  root <- report_repo_root(baseDir)
  reportDir <- file.path(root, "extras", "report")
  list(
    root = root,
    report = reportDir,
    data = file.path(reportDir, "data"),
    output = file.path(reportDir, "output"),
    outputTables = file.path(reportDir, "output", "tables"),
    outputFigures = file.path(reportDir, "output", "figures"),
    outputMetadata = file.path(reportDir, "output", "metadata"),
    outputReport = file.path(reportDir, "output", "report"),
    quarto = file.path(reportDir, "quarto")
  )
}

report_ensure_directories <- function(baseDir = getwd()) {
  paths <- report_paths(baseDir)
  for (path in unname(paths[c(
    "report",
    "data",
    "output",
    "outputTables",
    "outputFigures",
    "outputMetadata",
    "outputReport",
    "quarto"
  )])) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(paths)
}

report_write_table_bundle <- function(df, fileStem, outputDir) {
  csvPath <- file.path(outputDir, paste0(fileStem, ".csv"))
  rdsPath <- file.path(outputDir, paste0(fileStem, ".rds"))
  utils::write.csv(df, csvPath, row.names = FALSE)
  saveRDS(df, rdsPath)
  invisible(list(csv = csvPath, rds = rdsPath))
}

report_skip_reason <- function(distinctValues,
                               totalFrequency,
                               minDistinctValues = 20,
                               minTotalFrequency = 100) {
  if (distinctValues < minDistinctValues) {
    return("insufficient_distinct_values")
  }
  if (totalFrequency < minTotalFrequency) {
    return("insufficient_total_frequency")
  }
  NA_character_
}

report_extract_table <- function(andromeda, tableName) {
  if (inherits(andromeda, "report_mock_andromeda")) {
    return(andromeda[[tableName]])
  }

  if (is.list(andromeda) && !is.null(andromeda[[tableName]])) {
    return(andromeda[[tableName]])
  }

  if (!report_has_andromeda()) {
    stop(
      "Package `Andromeda` is required to read a real Andromeda object in this environment.",
      call. = FALSE
    )
  }

  table <- tryCatch(andromeda[[tableName]], error = function(e) NULL)
  if (is.null(table)) {
    stop(
      paste0("Table `", tableName, "` was not found in the supplied Andromeda object."),
      call. = FALSE
    )
  }

  as.data.frame(table, stringsAsFactors = FALSE)
}

report_first_non_missing <- function(x) {
  idx <- which(!is.na(x))[1]
  if (length(idx) == 0 || is.na(idx)) {
    return(NA)
  }
  x[idx]
}

report_get_bounds <- function(resultDf) {
  if (all(c("lowerFence", "upperFence") %in% names(resultDf))) {
    return(c(lower = report_first_non_missing(resultDf$lowerFence), upper = report_first_non_missing(resultDf$upperFence)))
  }
  if (all(c("lowerThreshold", "upperThreshold") %in% names(resultDf))) {
    return(c(lower = report_first_non_missing(resultDf$lowerThreshold), upper = report_first_non_missing(resultDf$upperThreshold)))
  }
  c(lower = NA_real_, upper = NA_real_)
}

report_source_package_functions <- function(baseDir = getwd()) {
  root <- report_repo_root(baseDir)
  rFiles <- list.files(
    file.path(root, "R"),
    pattern = "\\.[Rr]$",
    full.names = TRUE
  )
  for (path in rFiles) {
    sys.source(path, envir = .GlobalEnv)
  }
  invisible(rFiles)
}

report_require_analysis_functions <- function() {
  required <- c(
    "tukeyFences",
    "quantileThresholds",
    "zScoreOutliers",
    "modifiedZScoreOutliers",
    "generalizedESDOutliers",
    "compareUnivariateOutlierResults",
    "learnPlausibleRangeConsensus",
    "summarizeStrataDiagnostics",
    "flagAnalyzableStrata",
    "gaussianAnomalyOutliers",
    "isolationForestOutliers",
    "plotWeightedHistogram",
    "plotAnomalyScoreDistribution",
    "plotAnomalyRank",
    "plotStratumDiagnostics",
    "plotGadDistanceQQ"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing) > 0) {
    report_source_package_functions()
    missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
  }
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing analysis function(s): ",
        paste(missing, collapse = ", "),
        ". Source the package R scripts before running the report pipeline."
      ),
      call. = FALSE
    )
  }
}
