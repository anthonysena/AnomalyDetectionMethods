library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(AnomalyDetectionMethods)

# Weighted example datasets from inst/examples/main.R
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

makeReviewData <- function(df, qLower, qUpper, qInterp, zCutoff, modZCutoff, esdMax, esdAlpha) {
  withFences <- tukeyFences(df)
  withQuantile <- quantileThresholds(
    df,
    lowerProb = qLower,
    upperProb = qUpper,
    interpolate = qInterp
  )
  withZscore <- zScoreOutliers(df, zCutoff = zCutoff)
  withModZscore <- modifiedZScoreOutliers(df, zCutoff = modZCutoff)
  withEsd <- generalizedESDOutliers(df, maxOutliers = esdMax, alpha = esdAlpha)

  dplyr::bind_rows(
    withFences %>%
      transmute(
        method = "Tukey Fences",
        value,
        frequency,
        isOutlier,
        score = NA_real_
      ),
    withQuantile %>%
      transmute(
        method = "Quantile Thresholds",
        value,
        frequency,
        isOutlier,
        score = NA_real_
      ),
    withZscore %>%
      transmute(
        method = "Z-Score",
        value,
        frequency,
        isOutlier,
        score = zScore
      ),
    withModZscore %>%
      transmute(
        method = "Modified Z-Score",
        value,
        frequency,
        isOutlier,
        score = modifiedZScore
      ),
    withEsd %>%
      transmute(
        method = "Generalized ESD",
        value,
        frequency,
        isOutlier,
        score = outlierCount
      )
  )
}

ui <- fluidPage(
  titlePanel("Anomaly Detection Interactive Review"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "dataset",
        "Dataset",
        choices = c("Baseline" = "baseline", "Intentional Outliers" = "outlier"),
        selected = "outlier"
      ),
      checkboxGroupInput(
        "methods",
        "Methods to display",
        choices = c(
          "Tukey Fences",
          "Quantile Thresholds",
          "Z-Score",
          "Modified Z-Score",
          "Generalized ESD"
        ),
        selected = c(
          "Tukey Fences",
          "Quantile Thresholds",
          "Z-Score",
          "Modified Z-Score",
          "Generalized ESD"
        )
      ),
      checkboxInput("showOutliersOnly", "Show outliers only", value = FALSE),
      hr(),
      h4("Method Parameters"),
      sliderInput("qLower", "Quantile lower", min = 0.001, max = 0.2, value = 0.015, step = 0.001),
      sliderInput("qUpper", "Quantile upper", min = 0.8, max = 0.999, value = 0.985, step = 0.001),
      checkboxInput("qInterp", "Quantile interpolate", value = TRUE),
      sliderInput("zCutoff", "Z-score cutoff", min = 1.5, max = 5, value = 2.5, step = 0.1),
      sliderInput("modZCutoff", "Modified Z-score cutoff", min = 2, max = 8, value = 3.5, step = 0.1),
      sliderInput("esdMax", "ESD max outliers", min = 1, max = 10, value = 5, step = 1),
      sliderInput("esdAlpha", "ESD alpha", min = 0.001, max = 0.2, value = 0.05, step = 0.001)
    ),
    mainPanel(
      plotlyOutput("outlierPlot", height = "650px"),
      hr(),
      DTOutput("summaryTable")
    )
  )
)

server <- function(input, output, session) {
  selectedData <- reactive({
    if (identical(input$dataset, "baseline")) {
      return(baselineDf)
    }
    outlierDf
  })

  reviewData <- reactive({
    makeReviewData(
      df = selectedData(),
      qLower = input$qLower,
      qUpper = input$qUpper,
      qInterp = input$qInterp,
      zCutoff = input$zCutoff,
      modZCutoff = input$modZCutoff,
      esdMax = input$esdMax,
      esdAlpha = input$esdAlpha
    )
  })

  filteredReviewData <- reactive({
    df <- reviewData()
    if (length(input$methods) > 0) {
      df <- df %>% filter(method %in% input$methods)
    } else {
      df <- df[0, ]
    }
    if (isTRUE(input$showOutliersOnly)) {
      df <- df %>% filter(isOutlier)
    }
    df
  })

  output$outlierPlot <- renderPlotly({
    df <- filteredReviewData()
    validate(need(nrow(df) > 0, "No rows to display for current filters."))

    p <- ggplot(
      df,
      aes(
        x = value,
        y = frequency,
        color = isOutlier,
        size = frequency,
        text = paste0(
          "Method: ", method,
          "<br>Value: ", value,
          "<br>Frequency: ", frequency,
          "<br>Outlier: ", isOutlier,
          "<br>Score: ", ifelse(is.na(score), "NA", round(score, 3))
        )
      )
    ) +
      geom_point(alpha = 0.85) +
      facet_wrap(~method, scales = "free_y") +
      scale_color_manual(values = c("FALSE" = "#2c7fb8", "TRUE" = "#d7191c")) +
      theme_minimal(base_size = 12) +
      labs(
        x = "Value",
        y = "Frequency",
        color = "Outlier",
        size = "Frequency"
      )

    ggplotly(p, tooltip = "text")
  })

  output$summaryTable <- renderDT({
    filteredReviewData() %>%
      arrange(method, desc(isOutlier), value) %>%
      mutate(score = ifelse(is.na(score), NA, round(score, 3))) %>%
      datatable(
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 15, autoWidth = TRUE)
      )
  })
}

shinyApp(ui = ui, server = server)
