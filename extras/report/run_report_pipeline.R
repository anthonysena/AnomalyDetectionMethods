run_report_pipeline <- function(baseDir = getwd(),
                                createSample = TRUE,
                                renderQuarto = FALSE,
                                useMockIfUnavailable = TRUE) {
  paths <- report_ensure_directories(baseDir)
  report_source_package_functions(baseDir)

  sys.source(file.path(paths$report, "helpers.R"), envir = .GlobalEnv)
  sys.source(file.path(paths$report, "sample_andromeda.R"), envir = .GlobalEnv)
  sys.source(file.path(paths$report, "pipeline_functions.R"), envir = .GlobalEnv)

  if (createSample) {
    andromeda <- write_sample_measurement_andromeda(
      baseDir = baseDir,
      useMockIfUnavailable = useMockIfUnavailable
    )
  } else {
    mockPath <- file.path(paths$data, "sample_andromeda_mock.rds")
    if (file.exists(mockPath)) {
      andromeda <- readRDS(mockPath)
    } else {
      stop(
        "No sample object is available under extras/report/data. Set `createSample = TRUE` first.",
        call. = FALSE
      )
    }
  }

  run_measurement_plausibility_pipeline(
    andromeda = andromeda,
    outputDir = paths$output,
    renderQuarto = renderQuarto,
    baseDir = baseDir
  )
}

if (sys.nframe() == 0) {
  source(file.path("extras", "report", "helpers.R"))
  source(file.path("extras", "report", "sample_andromeda.R"))
  source(file.path("extras", "report", "pipeline_functions.R"))
  run_report_pipeline()
}
