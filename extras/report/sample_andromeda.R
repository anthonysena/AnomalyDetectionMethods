create_sample_measurement_concept_values <- function() {
  genders <- c(8507, 8532)
  careSites <- c(101, 102, 103)
  ageGrid <- list(
    `18-39` = c(22, 28, 34),
    `40-64` = c(45, 54, 63),
    `65-79` = c(67, 73, 78),
    `80+` = c(82, 88)
  )
  offsets <- c(-10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 10)
  weightMap <- c(3, 5, 7, 10, 13, 16, 20, 16, 13, 10, 7, 5, 3)
  names(weightMap) <- as.character(offsets)

  rows <- list()
  idx <- 1L
  for (bandName in names(ageGrid)) {
    for (age in ageGrid[[bandName]]) {
      for (gender in genders) {
        for (careSite in careSites) {
          baseMean <- 104 + floor((age - 18) / 8) + if (gender == 8507) 3 else 0
          siteShift <- if (careSite == 103) 1 else 0
          for (offset in offsets) {
            value <- baseMean + siteShift + offset
            frequency <- weightMap[[as.character(offset)]] + if (careSite == 102) 1 else 0
            rows[[idx]] <- data.frame(
              age = age,
              gender_concept_id = gender,
              care_site_id = careSite,
              value_as_number = value,
              frequency = frequency,
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
      }
    }
  }

  weakStratum <- data.frame(
    age = c(83, 83, 83, 83),
    gender_concept_id = c(8551, 8551, 8551, 8551),
    care_site_id = c(201, 201, 202, 202),
    value_as_number = c(96, 97, 98, 130),
    frequency = c(2, 2, 1, 1),
    stringsAsFactors = FALSE
  )

  globalAnomalies <- data.frame(
    age = c(25, 88),
    gender_concept_id = c(8507, 8532),
    care_site_id = c(103, 102),
    value_as_number = c(260, 35),
    frequency = c(1, 1),
    stringsAsFactors = FALSE
  )

  contextualAnomaly <- data.frame(
    age = c(24, 24),
    gender_concept_id = c(8532, 8532),
    care_site_id = c(101, 102),
    value_as_number = c(142, 144),
    frequency = c(2, 1),
    stringsAsFactors = FALSE
  )

  sampleDf <- do.call(
    rbind,
    c(rows, list(weakStratum, globalAnomalies, contextualAnomaly))
  )

  sampleDf[order(
    sampleDf$age,
    sampleDf$gender_concept_id,
    sampleDf$care_site_id,
    sampleDf$value_as_number
  ), ]
}

create_sample_measurement_andromeda <- function(path = NULL,
                                                useMockIfUnavailable = TRUE) {
  measurementDf <- create_sample_measurement_concept_values()

  if (report_has_andromeda()) {
    if (is.null(path)) {
      path <- tempfile("sample_andromeda_")
    }
    if (dir.exists(path)) {
      unlink(path, recursive = TRUE, force = TRUE)
    }
    andromeda <- Andromeda::andromeda(path)
    andromeda$measurement_concept_values <- measurementDf
    return(andromeda)
  }

  if (!useMockIfUnavailable) {
    stop(
      "Package `Andromeda` is not installed, so a real Andromeda object cannot be created.",
      call. = FALSE
    )
  }

  structure(
    list(
      measurement_concept_values = measurementDf,
      dataPath = path
    ),
    class = c("report_mock_andromeda", "list")
  )
}

write_sample_measurement_andromeda <- function(path = NULL,
                                               baseDir = getwd(),
                                               useMockIfUnavailable = TRUE) {
  paths <- report_ensure_directories(baseDir)
  if (is.null(path)) {
    path <- file.path(paths$data, "sample_andromeda")
  }

  andromeda <- create_sample_measurement_andromeda(
    path = path,
    useMockIfUnavailable = useMockIfUnavailable
  )

  measurementDf <- report_extract_table(andromeda, "measurement_concept_values")
  utils::write.csv(
    measurementDf,
    file.path(paths$data, "sample_measurement_concept_values.csv"),
    row.names = FALSE
  )

  if (inherits(andromeda, "report_mock_andromeda")) {
    saveRDS(andromeda, file.path(paths$data, "sample_andromeda_mock.rds"))
  }

  invisible(andromeda)
}
