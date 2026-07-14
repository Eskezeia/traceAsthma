# TRACE-Asthma point-of-care Shiny app
#
# A minimal clinician-facing interface: load a deployed traceAsthma_model
# (.rds file, produced by build_deployable_model() / run_trace_asthma_pipeline()),
# upload a single patient's molecular data, and view the predicted risk
# together with the score breakdown so the result is interpretable rather
# than a black-box number.
#
# Launch with: shiny::runApp(system.file("shiny", package = "traceAsthma"))

library(shiny)
library(traceAsthma)

ui <- fluidPage(
  titlePanel("TRACE-Asthma risk calculator (research use only)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("model_file", "Deployed model (.rds)", accept = ".rds"),
      fileInput("methylation_file", "Patient methylation (CSV: CpG x patient)", accept = ".csv"),
      fileInput("expression_file", "Patient expression (CSV: gene x patient)", accept = ".csv"),
      numericInput("age", "Age", value = NA, min = 0, max = 120),
      selectInput("sex", "Sex", choices = c("", "F", "M")),
      numericInput("bmi", "BMI", value = NA, min = 0),
      selectInput("smoking", "Smoking status", choices = c("", "never", "former", "current")),
      actionButton("run", "Compute risk", class = "btn-primary"),
      tags$hr(),
      tags$p(tags$strong("Disclaimer: "),
             "Research tool only. Not a diagnostic device. Predictions require ",
             "site-specific validation before any clinical use. See package vignette ",
             "'clinical-deployment-notes'.")
    ),
    mainPanel(
      h4("Predicted risk"),
      tableOutput("risk_table"),
      h4("Score breakdown"),
      plotOutput("score_plot", height = "260px"),
      h4("Model provenance"),
      verbatimTextOutput("model_info")
    )
  )
)

server <- function(input, output, session) {
  model <- eventReactive(input$model_file, {
    readRDS(input$model_file$datapath)
  })

  output$model_info <- renderPrint({
    req(model())
    print(model())
  })

  result <- eventReactive(input$run, {
    req(model(), input$methylation_file, input$expression_file)

    meth <- as.matrix(read.csv(input$methylation_file$datapath, row.names = 1, check.names = FALSE))
    expr <- as.matrix(read.csv(input$expression_file$datapath, row.names = 1, check.names = FALSE))

    clinical <- data.frame(
      age = input$age,
      sex = input$sex,
      bmi = input$bmi,
      smoking = input$smoking
    )
    clinical <- clinical[, colnames(clinical) %in% names(stats::coef(model()$risk_model)), drop = FALSE]

    predict_asthma_risk(model(), methylation = meth, expression = expr,
                         clinical = if (ncol(clinical) > 0) clinical else NULL)
  })

  output$risk_table <- renderTable({
    req(result())
    result()[, c("patient_id", "predicted_probability", "risk_category")]
  })

  output$score_plot <- renderPlot({
    req(result())
    r <- result()
    scores <- c(eQTM = r$eqtm_score[1], TRACE = r$trace_score[1])
    if ("mprs_score" %in% names(r)) scores <- c(MPRS = r$mprs_score[1], scores)
    barplot(scores, col = ifelse(scores > 0, "#D85A30", "#378ADD"),
            ylab = "Standardized score (Z)", main = "Patient molecular score profile")
    abline(h = 0, lty = 2)
  })
}

shinyApp(ui, server)
