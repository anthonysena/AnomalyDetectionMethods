prepare_measurement_concept_values <- function(andromeda,
                                               tableName = "measurement_concept_values",
                                               ageColumn = "age",
                                               genderColumn = "gender_concept_id",
                                               careSiteColumn = "care_site_id",
                                               valueColumn = "value_as_number",
                                               frequencyColumn = "frequency",
                                               ageBreaks = c(18, 40, 65, 80, Inf),
                                               ageLabels = c("18-39", "40-64", "65-79", "80+")) {
  report_check_required_packages(c("ggplot2"))
  rawDf <- report_extract_table(andromeda, tableName)
  requiredColumns <- c(ageColumn, genderColumn, careSiteColumn, valueColumn, frequencyColumn)
  missingColumns <- setdiff(requiredColumns, names(rawDf))
  if (length(missingColumns) > 0) {
    stop(
      paste0("Missing required column(s): ", paste(missingColumns, collapse = ", ")),
      call. = FALSE
    )
  }

  prepared <- rawDf[, requiredColumns, drop = FALSE]
  names(prepared) <- c(
    "age",
    "gender_concept_id",
    "care_site_id",
    "value",
    "frequency"
  )

  prepared$age <- as.numeric(prepared$age)
  prepared$gender_concept_id <- as.numeric(prepared$gender_concept_id)
  prepared$care_site_id <- as.numeric(prepared$care_site_id)
  prepared$value <- as.numeric(prepared$value)
  prepared$frequency <- as.numeric(prepared$frequency)

  initialRows <- nrow(prepared)
  keep <- !is.na(prepared$value) &
    !is.na(prepared$age) &
    !is.na(prepared$gender_concept_id) &
    !is.na(prepared$care_site_id) &
    !is.na(prepared$frequency) &
    prepared$frequency > 0
  excludedRows <- sum(!keep)
  prepared <- prepared[keep, , drop = FALSE]

  if (nrow(prepared) < 1) {
    stop("No valid rows remain after filtering invalid measurements.", call. = FALSE)
  }

  prepared$age_band <- as.character(
    report_age_band(
      age = prepared$age,
      breaks = ageBreaks,
      labels = ageLabels
    )
  )
  prepared <- prepared[!is.na(prepared$age_band), , drop = FALSE]

  genderDf <- report_gender_map(prepared$gender_concept_id)
  prepared <- merge(
    prepared,
    genderDf,
    by = "gender_concept_id",
    all.x = TRUE,
    sort = FALSE
  )
  prepared$gender_numeric <- as.numeric(report_encode_gender(prepared$gender_concept_id))

  prepared <- stats::aggregate(
    frequency ~ age + age_band + gender_concept_id + gender_label + gender_numeric + care_site_id + value,
    data = prepared,
    FUN = sum
  )

  prepared$row_id <- seq_len(nrow(prepared))
  attr(prepared, "excluded_rows") <- excludedRows
  attr(prepared, "input_rows") <- initialRows
  prepared[order(prepared$row_id), ]
}

report_run_univariate_methods <- function(valueDf) {
  list(
    tukey = tukeyFences(valueDf),
    quantile = quantileThresholds(valueDf, lowerProb = 0.01, upperProb = 0.99),
    zscore = zScoreOutliers(valueDf, zCutoff = 3),
    modified = modifiedZScoreOutliers(valueDf, zCutoff = 3.5),
    esd = generalizedESDOutliers(valueDf, maxOutliers = 10, alpha = 0.05)
  )
}

run_global_univariate_analysis <- function(analysisDf) {
  globalDf <- stats::aggregate(
    frequency ~ value,
    data = analysisDf[, c("value", "frequency")],
    FUN = sum
  )
  globalDf <- globalDf[order(globalDf$value), , drop = FALSE]

  methods <- report_run_univariate_methods(globalDf)
  comparison <- compareUnivariateOutlierResults(methods)
  consensus <- learnPlausibleRangeConsensus(
    comparisonDf = comparison,
    consensusMethod = "count",
    consensusThreshold = 3,
    rangeMethod = "weighted_quantile",
    lowerProb = 0.001,
    upperProb = 0.999
  )

  rangeRows <- lapply(names(methods), function(methodName) {
    bounds <- report_get_bounds(methods[[methodName]])
    data.frame(
      method = methodName,
      lower_bound = bounds[["lower"]],
      upper_bound = bounds[["upper"]],
      flagged_distinct_values = sum(methods[[methodName]]$isOutlier),
      flagged_total_frequency = sum(methods[[methodName]]$frequency[methods[[methodName]]$isOutlier]),
      stringsAsFactors = FALSE
    )
  })

  list(
    aggregated_values = globalDf,
    method_results = methods,
    comparison = comparison,
    consensus = consensus,
    method_ranges = do.call(rbind, rangeRows)
  )
}

run_stratified_univariate_analysis <- function(analysisDf,
                                               minDistinctValues = 20,
                                               minTotalFrequency = 100) {
  diagnostics <- suppressWarnings(
    summarizeStrataDiagnostics(
      df = analysisDf,
      strataColumns = c("age_band", "gender_concept_id", "gender_label"),
      valueColumn = "value",
      frequencyColumn = "frequency",
      featureColumns = c("value", "age", "gender_numeric")
    )
  )

  diagnostics$stratumAnalyzable <- diagnostics$distinctValues >= minDistinctValues &
    diagnostics$totalFrequency >= minTotalFrequency
  diagnostics$skipReason <- mapply(
    report_skip_reason,
    diagnostics$distinctValues,
    diagnostics$totalFrequency,
    MoreArgs = list(
      minDistinctValues = minDistinctValues,
      minTotalFrequency = minTotalFrequency
    )
  )

  splitKeys <- interaction(
    analysisDf[, c("age_band", "gender_concept_id", "gender_label"), drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  pieces <- split(analysisDf, splitKeys)

  comparisonRows <- list()
  rangeRows <- list()
  idxComparison <- 1L
  idxRanges <- 1L

  for (piece in pieces) {
    key <- piece[1, c("age_band", "gender_concept_id", "gender_label"), drop = FALSE]
    diagRow <- diagnostics[
      diagnostics$age_band == key$age_band &
        diagnostics$gender_concept_id == key$gender_concept_id &
        diagnostics$gender_label == key$gender_label,
      ,
      drop = FALSE
    ]

    if (!isTRUE(diagRow$stratumAnalyzable[1])) {
      rangeRows[[idxRanges]] <- cbind(
        key,
        data.frame(
          consensusMethod = "count",
          consensusThreshold = 3,
          rangeMethod = "weighted_quantile",
          lowerProb = 0.001,
          upperProb = 0.999,
          minPlausible = NA_real_,
          maxPlausible = NA_real_,
          retainedDistinctValues = NA_real_,
          retainedTotalFrequency = NA_real_,
          excludedDistinctValues = NA_real_,
          excludedTotalFrequency = NA_real_,
          stratumAnalyzable = FALSE,
          skipReason = diagRow$skipReason[1],
          stringsAsFactors = FALSE
        )
      )
      idxRanges <- idxRanges + 1L
      next
    }

    stratumDf <- stats::aggregate(
      frequency ~ value,
      data = piece[, c("value", "frequency")],
      FUN = sum
    )
    stratumDf <- stratumDf[order(stratumDf$value), , drop = FALSE]

    methods <- report_run_univariate_methods(stratumDf)
    comparison <- compareUnivariateOutlierResults(methods)
    consensus <- learnPlausibleRangeConsensus(
      comparisonDf = comparison,
      consensusMethod = "count",
      consensusThreshold = 3,
      rangeMethod = "weighted_quantile",
      lowerProb = 0.001,
      upperProb = 0.999
    )

    comparisonRows[[idxComparison]] <- cbind(
      key[rep(1, nrow(consensus$valueSummary)), , drop = FALSE],
      consensus$valueSummary,
      stratumAnalyzable = TRUE,
      skipReason = NA_character_,
      stringsAsFactors = FALSE
    )
    idxComparison <- idxComparison + 1L

    rangeRows[[idxRanges]] <- cbind(
      key,
      consensus$rangeSummary,
      stratumAnalyzable = TRUE,
      skipReason = NA_character_,
      stringsAsFactors = FALSE
    )
    idxRanges <- idxRanges + 1L
  }

  list(
    diagnostics = diagnostics,
    comparison = if (length(comparisonRows) > 0) do.call(rbind, comparisonRows) else data.frame(),
    ranges = do.call(rbind, rangeRows)
  )
}

run_contextual_models <- function(analysisDf) {
  gadOut <- gaussianAnomalyOutliers(
    df = analysisDf,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("value", "age", "gender_numeric"),
    transformMethod = "weighted_zscore",
    covarianceMethod = "classical",
    tailProb = 0.99
  )

  ifOut <- isolationForestOutliers(
    df = analysisDf,
    valueColumn = "value",
    frequencyColumn = "frequency",
    featureColumns = c("value", "age", "gender_numeric", "frequency"),
    ntrees = 500,
    tailProb = 0.99,
    weightsAsSampleProb = FALSE,
    standardizeData = FALSE,
    seed = 1,
    nthreads = 1
  )

  contextual <- merge(
    gadOut[, c(
      "row_id",
      "gaussianLogDensity",
      "gaussianDensity",
      "mahalanobisDistance",
      "distanceCutoff",
      "outlierProportion",
      "isOutlier"
    )],
    ifOut[, c(
      "row_id",
      "isolationForestScore",
      "scoreCutoff",
      "outlierProportion",
      "isOutlier"
    )],
    by = "row_id",
    suffixes = c("_gad", "_if"),
    all = TRUE,
    sort = FALSE
  )

  names(contextual)[names(contextual) == "isOutlier_gad"] <- "gaussianFlag"
  names(contextual)[names(contextual) == "isOutlier_if"] <- "isolationForestFlag"
  names(contextual)[names(contextual) == "outlierProportion_gad"] <- "gaussianOutlierProportion"
  names(contextual)[names(contextual) == "outlierProportion_if"] <- "isolationForestOutlierProportion"

  contextual
}

classify_measurement_rows <- function(analysisDf,
                                      globalAnalysis,
                                      stratifiedAnalysis,
                                      contextualAnalysis) {
  globalFlags <- globalAnalysis$consensus$valueSummary[, c(
    "value",
    "consensusCount",
    "isConsensusOutlier_count"
  )]
  names(globalFlags) <- c(
    "value",
    "globalConsensusCount",
    "globalConsensusFlag"
  )

  stratifiedFlags <- stratifiedAnalysis$comparison
  if (nrow(stratifiedFlags) > 0) {
    stratifiedFlags <- stratifiedFlags[, c(
      "age_band",
      "gender_concept_id",
      "gender_label",
      "value",
      "consensusCount",
      "isConsensusOutlier_count",
      "stratumAnalyzable",
      "skipReason"
    )]
    names(stratifiedFlags) <- c(
      "age_band",
      "gender_concept_id",
      "gender_label",
      "value",
      "stratumConsensusCount",
      "stratumConsensusFlag",
      "stratumAnalyzable",
      "skipReason"
    )
  } else {
    stratifiedFlags <- data.frame(
      age_band = character(0),
      gender_concept_id = numeric(0),
      gender_label = character(0),
      value = numeric(0),
      stratumConsensusCount = numeric(0),
      stratumConsensusFlag = logical(0),
      stratumAnalyzable = logical(0),
      skipReason = character(0),
      stringsAsFactors = FALSE
    )
  }

  analyzability <- stratifiedAnalysis$diagnostics[, c(
    "age_band",
    "gender_concept_id",
    "gender_label",
    "stratumAnalyzable",
    "skipReason"
  )]

  merged <- merge(analysisDf, globalFlags, by = "value", all.x = TRUE, sort = FALSE)
  merged <- merge(
    merged,
    stratifiedFlags,
    by = c("age_band", "gender_concept_id", "gender_label", "value"),
    all.x = TRUE,
    sort = FALSE
  )
  merged <- merge(
    merged,
    analyzability,
    by = c("age_band", "gender_concept_id", "gender_label"),
    all.x = TRUE,
    suffixes = c("", "_diagnostic"),
    sort = FALSE
  )
  merged <- merge(merged, contextualAnalysis, by = "row_id", all.x = TRUE, sort = FALSE)

  if ("stratumAnalyzable_diagnostic" %in% names(merged)) {
    merged$stratumAnalyzable <- ifelse(
      is.na(merged$stratumAnalyzable),
      merged$stratumAnalyzable_diagnostic,
      merged$stratumAnalyzable
    )
    merged$skipReason <- ifelse(
      is.na(merged$skipReason),
      merged$skipReason_diagnostic,
      merged$skipReason
    )
    merged$stratumAnalyzable_diagnostic <- NULL
    merged$skipReason_diagnostic <- NULL
  }

  merged$globalConsensusFlag[is.na(merged$globalConsensusFlag)] <- FALSE
  merged$stratumConsensusFlag[is.na(merged$stratumConsensusFlag)] <- FALSE
  merged$stratumAnalyzable[is.na(merged$stratumAnalyzable)] <- FALSE
  merged$gaussianFlag[is.na(merged$gaussianFlag)] <- FALSE
  merged$isolationForestFlag[is.na(merged$isolationForestFlag)] <- FALSE

  merged$finalClassification <- "plausible"
  merged$finalClassification[!merged$stratumAnalyzable &
    !merged$globalConsensusFlag &
    !merged$gaussianFlag &
    !merged$isolationForestFlag] <- "not_evaluable_stratum"
  merged$finalClassification[!merged$globalConsensusFlag &
    (merged$stratumConsensusFlag | (merged$gaussianFlag & merged$isolationForestFlag))] <- "context_dependent_anomaly"
  merged$finalClassification[!merged$globalConsensusFlag &
    merged$finalClassification == "plausible" &
    xor(merged$gaussianFlag, merged$isolationForestFlag)] <- "borderline_disagreement"
  merged$finalClassification[merged$globalConsensusFlag] <- "globally_implausible"

  merged
}

build_report_tables <- function(analysisDf,
                                globalAnalysis,
                                stratifiedAnalysis,
                                contextualAnalysis,
                                rowLevelResults) {
  inputSummary <- data.frame(
    input_rows = attr(analysisDf, "input_rows"),
    excluded_rows = attr(analysisDf, "excluded_rows"),
    analyzed_rows = nrow(analysisDf),
    distinct_values = length(unique(analysisDf$value)),
    total_weighted_frequency = sum(analysisDf$frequency),
    age_min = min(analysisDf$age),
    age_max = max(analysisDf$age),
    stringsAsFactors = FALSE
  )

  ageBandDistribution <- stats::aggregate(
    frequency ~ age_band,
    data = analysisDf[, c("age_band", "frequency")],
    FUN = sum
  )
  genderDistribution <- stats::aggregate(
    frequency ~ gender_concept_id + gender_label,
    data = analysisDf[, c("gender_concept_id", "gender_label", "frequency")],
    FUN = sum
  )

  contextualSummary <- data.frame(
    gaussian_cutoff = report_first_non_missing(contextualAnalysis$distanceCutoff),
    gaussian_flagged_rows = sum(contextualAnalysis$gaussianFlag, na.rm = TRUE),
    gaussian_flagged_frequency = sum(analysisDf$frequency[contextualAnalysis$gaussianFlag], na.rm = TRUE),
    isolation_forest_cutoff = report_first_non_missing(contextualAnalysis$scoreCutoff),
    isolation_forest_flagged_rows = sum(contextualAnalysis$isolationForestFlag, na.rm = TRUE),
    isolation_forest_flagged_frequency = sum(analysisDf$frequency[contextualAnalysis$isolationForestFlag], na.rm = TRUE),
    overlap_rows = sum(contextualAnalysis$gaussianFlag & contextualAnalysis$isolationForestFlag, na.rm = TRUE),
    overlap_frequency = sum(analysisDf$frequency[contextualAnalysis$gaussianFlag & contextualAnalysis$isolationForestFlag], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  classificationSummary <- stats::aggregate(
    frequency ~ finalClassification,
    data = rowLevelResults[, c("finalClassification", "frequency")],
    FUN = sum
  )
  classificationSummary <- classificationSummary[order(-classificationSummary$frequency), , drop = FALSE]

  list(
    input_summary = inputSummary,
    age_band_distribution = ageBandDistribution,
    gender_distribution = genderDistribution,
    global_method_ranges = globalAnalysis$method_ranges,
    global_consensus_range = globalAnalysis$consensus$rangeSummary,
    global_value_summary = globalAnalysis$consensus$valueSummary,
    strata_diagnostics = stratifiedAnalysis$diagnostics,
    stratified_ranges = stratifiedAnalysis$ranges,
    contextual_model_summary = contextualSummary,
    row_level_results = rowLevelResults,
    classification_summary = classificationSummary
  )
}

build_report_figures <- function(analysisDf,
                                 globalAnalysis,
                                 stratifiedAnalysis,
                                 contextualAnalysis,
                                 rowLevelResults) {
  globalRange <- globalAnalysis$consensus$rangeSummary[1, , drop = FALSE]
  mergedContext <- merge(analysisDf, contextualAnalysis, by = "row_id", all.x = TRUE, sort = FALSE)

  globalHistogram <- plotWeightedHistogram(globalAnalysis$aggregated_values, addNormal = TRUE) +
    ggplot2::geom_vline(
      xintercept = c(globalRange$minPlausible[1], globalRange$maxPlausible[1]),
      linetype = 2,
      color = "red"
    ) +
    ggplot2::labs(title = "Global Weighted Distribution")

  methodRangesDf <- globalAnalysis$method_ranges[
    !is.na(globalAnalysis$method_ranges$lower_bound) &
      !is.na(globalAnalysis$method_ranges$upper_bound),
    ,
    drop = FALSE
  ]
  methodRangesPlot <- ggplot2::ggplot(
    methodRangesDf,
    ggplot2::aes(
      y = stats::reorder(.data$method, .data$lower_bound),
      x = .data$lower_bound,
      xend = .data$upper_bound,
      yend = .data$method
    )
  ) +
    ggplot2::geom_segment(linewidth = 1.2, color = "steelblue") +
    ggplot2::geom_point(size = 2.4, color = "steelblue") +
    ggplot2::geom_point(ggplot2::aes(x = .data$upper_bound), size = 2.4, color = "steelblue") +
    ggplot2::labs(x = "Value", y = "Method", title = "Global Method Bounds")

  consensusPlot <- ggplot2::ggplot(
    globalAnalysis$consensus$valueSummary,
    ggplot2::aes(x = .data$value, y = .data$frequency, color = .data$consensusCount)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_gradient(low = "grey70", high = "firebrick") +
    ggplot2::labs(
      x = "Value",
      y = "Frequency",
      color = "Consensus Count",
      title = "Consensus Count by Value"
    )

  diagnosticsPlot <- plotStratumDiagnostics(
    stratifiedAnalysis$diagnostics,
    xColumn = "distinctValues",
    yColumn = "totalFrequency",
    colorColumn = "stratumAnalyzable"
  ) + ggplot2::labs(title = "Stratum Diagnostics")

  analyzableRanges <- stratifiedAnalysis$ranges[stratifiedAnalysis$ranges$stratumAnalyzable, , drop = FALSE]
  if (nrow(analyzableRanges) > 0) {
    analyzableRanges$stratumLabel <- paste(
      analyzableRanges$age_band,
      analyzableRanges$gender_label,
      sep = " / "
    )
    stratifiedRangesPlot <- ggplot2::ggplot(
      analyzableRanges,
      ggplot2::aes(
        y = stats::reorder(.data$stratumLabel, .data$minPlausible),
        x = .data$minPlausible,
        xend = .data$maxPlausible,
        yend = .data$stratumLabel,
        color = .data$gender_label
      )
    ) +
      ggplot2::geom_segment(linewidth = 1.1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = .data$maxPlausible), size = 2) +
      ggplot2::labs(x = "Value", y = "Stratum", title = "Stratified Plausible Ranges")
  } else {
    stratifiedRangesPlot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "No analyzable strata")
  }

  gadScorePlot <- plotAnomalyScoreDistribution(
    mergedContext,
    scoreColumn = "mahalanobisDistance",
    frequencyColumn = "frequency",
    cutoffColumn = "distanceCutoff"
  ) + ggplot2::labs(title = "Gaussian Anomaly Score Distribution")

  gadQqPlot <- plotGadDistanceQQ(
    mergedContext,
    distanceColumn = "mahalanobisDistance",
    frequencyColumn = "frequency",
    featureColumns = c("value", "age", "gender_numeric")
  ) + ggplot2::labs(title = "Gaussian Distance Q-Q Plot")

  ifScorePlot <- plotAnomalyRank(
    mergedContext,
    scoreColumn = "isolationForestScore",
    frequencyColumn = "frequency",
    cutoffColumn = "scoreCutoff"
  ) + ggplot2::labs(title = "Isolation Forest Ranked Scores")

  contextualScatter <- ggplot2::ggplot(
    rowLevelResults,
    ggplot2::aes(x = .data$age, y = .data$value, size = .data$frequency, color = .data$finalClassification)
  ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::labs(x = "Age", y = "Value", title = "Contextual Value Scatter")

  classificationSummary <- stats::aggregate(
    frequency ~ finalClassification,
    data = rowLevelResults[, c("finalClassification", "frequency")],
    FUN = sum
  )
  classificationPlot <- ggplot2::ggplot(
    classificationSummary,
    ggplot2::aes(
      x = stats::reorder(.data$finalClassification, .data$frequency),
      y = .data$frequency,
      fill = .data$finalClassification
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Classification", y = "Weighted Frequency", title = "Final Classification Summary")

  list(
    global_weighted_distribution = globalHistogram,
    global_method_bounds = methodRangesPlot,
    consensus_count_by_value = consensusPlot,
    stratum_diagnostics = diagnosticsPlot,
    stratified_ranges = stratifiedRangesPlot,
    gad_score_distribution = gadScorePlot,
    gad_distance_qq = gadQqPlot,
    isolation_forest_rank = ifScorePlot,
    contextual_value_scatter = contextualScatter,
    final_classification_bar = classificationPlot
  )
}

write_report_outputs <- function(pipelineResults, outputDir) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  tableDir <- file.path(outputDir, "tables")
  figureDir <- file.path(outputDir, "figures")
  metadataDir <- file.path(outputDir, "metadata")
  reportDir <- file.path(outputDir, "report")
  dir.create(tableDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figureDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadataDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reportDir, recursive = TRUE, showWarnings = FALSE)

  for (name in names(pipelineResults$tables)) {
    report_write_table_bundle(pipelineResults$tables[[name]], name, tableDir)
  }

  for (name in names(pipelineResults$figures)) {
    ggplot2::ggsave(
      filename = file.path(figureDir, paste0(name, ".png")),
      plot = pipelineResults$figures[[name]],
      width = 10,
      height = 6,
      dpi = 300
    )
  }

  metadata <- list(
    generated_at = as.character(Sys.time()),
    global_consensus_threshold = 3,
    range_method = "weighted_quantile",
    lower_prob = 0.001,
    upper_prob = 0.999
  )
  saveRDS(metadata, file.path(metadataDir, "run_metadata.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(metadataDir, "session_info.txt"))

  invisible(
    list(
      tables = tableDir,
      figures = figureDir,
      metadata = metadataDir,
      report = reportDir
    )
  )
}

render_report <- function(outputDir,
                          reportTitle = "Measurement Plausibility Summary Report",
                          baseDir = getwd()) {
  paths <- report_paths(baseDir)
  qmdPath <- file.path(paths$quarto, "measurement_plausibility_report.qmd")
  outputFile <- file.path(outputDir, "report", "measurement_plausibility_report.html")
  finalizeOutput <- function() {
    if (file.exists(outputFile)) {
      return(invisible(outputFile))
    }

    candidates <- list.files(
      path = paths$root,
      pattern = paste0("^", basename(outputFile), "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    candidates <- setdiff(normalizePath(candidates, winslash = "/", mustWork = FALSE),
      normalizePath(outputFile, winslash = "/", mustWork = FALSE)
    )
    if (length(candidates) > 0) {
      newest <- candidates[which.max(file.info(candidates)$mtime)]
      dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)
      ok <- file.rename(newest, outputFile)
      if (!ok) {
        file.copy(newest, outputFile, overwrite = TRUE)
        unlink(newest, force = TRUE)
      }
    }

    if (!file.exists(outputFile)) {
      stop("Rendered report HTML was not found after Quarto completed.", call. = FALSE)
    }

    invisible(outputFile)
  }

  if (report_has_quarto_cli()) {
    oldOutputDir <- Sys.getenv("REPORT_OUTPUT_DIR", unset = NA_character_)
    oldTitle <- Sys.getenv("REPORT_TITLE", unset = NA_character_)
    on.exit({
      if (is.na(oldOutputDir)) {
        Sys.unsetenv("REPORT_OUTPUT_DIR")
      } else {
        Sys.setenv(REPORT_OUTPUT_DIR = oldOutputDir)
      }
      if (is.na(oldTitle)) {
        Sys.unsetenv("REPORT_TITLE")
      } else {
        Sys.setenv(REPORT_TITLE = oldTitle)
      }
    }, add = TRUE)
    Sys.setenv(
      REPORT_OUTPUT_DIR = normalizePath(outputDir, winslash = "/", mustWork = TRUE),
      REPORT_TITLE = reportTitle
    )
    cmd <- c(
      "render",
      qmdPath,
      "--to",
      "html",
      "--output",
      basename(outputFile),
      "--output-dir",
      dirname(outputFile)
    )
    status <- system2(Sys.which("quarto"), cmd)
    if (status != 0) {
      stop("Quarto CLI failed to render the report.", call. = FALSE)
    }
    return(finalizeOutput())
  }

  if (report_has_quarto_package()) {
    oldWd <- getwd()
    oldOutputDir <- Sys.getenv("REPORT_OUTPUT_DIR", unset = NA_character_)
    oldTitle <- Sys.getenv("REPORT_TITLE", unset = NA_character_)
    on.exit(setwd(oldWd), add = TRUE)
    on.exit({
      if (is.na(oldOutputDir)) {
        Sys.unsetenv("REPORT_OUTPUT_DIR")
      } else {
        Sys.setenv(REPORT_OUTPUT_DIR = oldOutputDir)
      }
      if (is.na(oldTitle)) {
        Sys.unsetenv("REPORT_TITLE")
      } else {
        Sys.setenv(REPORT_TITLE = oldTitle)
      }
    }, add = TRUE)
    setwd(dirname(outputFile))
    Sys.setenv(
      REPORT_OUTPUT_DIR = normalizePath(outputDir, winslash = "/", mustWork = TRUE),
      REPORT_TITLE = reportTitle
    )
    quarto::quarto_render(
      input = qmdPath,
      output_file = basename(outputFile),
      quiet = TRUE
    )
    return(finalizeOutput())
  }

  stop(
    "Quarto is not available. Install the `quarto` R package or the Quarto CLI to render the report.",
    call. = FALSE
  )
}

run_measurement_plausibility_pipeline <- function(andromeda,
                                                  outputDir = NULL,
                                                  renderQuarto = FALSE,
                                                  baseDir = getwd()) {
  report_check_required_packages(c("ggplot2"))
  report_require_analysis_functions()

  analysisDf <- prepare_measurement_concept_values(andromeda)
  globalAnalysis <- run_global_univariate_analysis(analysisDf)
  stratifiedAnalysis <- run_stratified_univariate_analysis(analysisDf)
  contextualAnalysis <- run_contextual_models(analysisDf)
  rowLevelResults <- classify_measurement_rows(
    analysisDf,
    globalAnalysis,
    stratifiedAnalysis,
    contextualAnalysis
  )
  tables <- build_report_tables(
    analysisDf,
    globalAnalysis,
    stratifiedAnalysis,
    contextualAnalysis,
    rowLevelResults
  )
  figures <- build_report_figures(
    analysisDf,
    globalAnalysis,
    stratifiedAnalysis,
    contextualAnalysis,
    rowLevelResults
  )

  out <- list(
    analysis_data = analysisDf,
    global = globalAnalysis,
    stratified = stratifiedAnalysis,
    contextual = contextualAnalysis,
    row_level_results = rowLevelResults,
    tables = tables,
    figures = figures
  )

  if (!is.null(outputDir)) {
    write_report_outputs(out, outputDir)
    if (renderQuarto) {
      render_report(outputDir = outputDir, baseDir = baseDir)
    }
  }

  out
}
