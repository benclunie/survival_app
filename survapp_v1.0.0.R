# Install / load packages

library(shiny)
library(survival)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(survminer)
library(AICcmodavg)

# User Interface
ui <- fluidPage(
  titlePanel("Survival Analysis and Kaplan-Meier Plots"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      tags$hr(),
      checkboxInput("header", "Header", TRUE),
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                   selected = ","),
      radioButtons("quote", "Quote",
                   choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'"),
                   selected = '"'),
      tags$hr(),
      textInput("timeVar", "Time Variable", "Time"),
      textInput("eventVar", "Event Variable", "Status"),
      textInput("treatmentVar", "Treatment Variable", "Treatment"),
      textInput("doseVar", "Dose Variable (optional)", ""),
      checkboxInput("facetByTreatment", "Facet by Treatment", value = FALSE),
      actionButton("analyze", "Analyze")
    ),
    mainPanel(
      plotOutput("kmPlot"),
      tableOutput("summaryTable")
    )
  )
)

server <- function(input, output) {
  data <- reactive({
    req(input$file)
    inFile <- input$file
    df <- read.csv(inFile$datapath, header = input$header, sep = input$sep, quote = input$quote)
    # Debugging: Print the first few rows of the data frame
    print("Data preview:")
    print(head(df))
    df
  })
  analysisResults <- eventReactive(input$analyze, {
    df <- data()
    timeVar <- input$timeVar
    eventVar <- input$eventVar
    treatmentVar <- input$treatmentVar
    doseVar <- input$doseVar
    # Ensure the specified columns exist in the data frame
    if (!(timeVar %in% names(df))) stop("Time variable not found in the data")
    if (!(eventVar %in% names(df))) stop("Event variable not found in the data")
    if (!(treatmentVar %in% names(df))) stop("Treatment variable not found in the data")
    if (doseVar != "" && !(doseVar %in% names(df))) stop("Dose variable not found in the data")
    # Convert variables to appropriate types
    df[[timeVar]] <- as.numeric(df[[timeVar]])
    df[[eventVar]] <- as.numeric(df[[eventVar]])
    df[[treatmentVar]] <- as.factor(df[[treatmentVar]])
    if (doseVar != "") {
      df[[doseVar]] <- as.factor(df[[doseVar]])
    }
    # Fit the survival models
    models <- list()
    model_names <- c()
    surv_T <- survreg(as.formula(paste("Surv(", timeVar, ",", eventVar, ") ~", treatmentVar)), data = df)
    models <- append(models, list(surv_T))
    model_names <- append(model_names, "Treatment")
    if (doseVar != "") {
      surv_T_x_D <- survreg(as.formula(paste("Surv(", timeVar, ",", eventVar, ") ~", treatmentVar, "*", doseVar)), data = df)
      surv_T_D <- survreg(as.formula(paste("Surv(", timeVar, ",", eventVar, ") ~", treatmentVar, "+", doseVar)), data = df)
      surv_D <- survreg(as.formula(paste("Surv(", timeVar, ",", eventVar, ") ~", doseVar)), data = df)
      models <- append(models, list(surv_T_x_D, surv_T_D, surv_D))
      model_names <- append(model_names, c("Treatment x Dose", "Treatment + Dose", "Dose"))
    }
    # Create an AIC table
    aic_table <- aictab(models, modnames = model_names)
    # Select the best model
    best_model_index <- which.min(aic_table$AICc)
    best_model <- models[[best_model_index]]
    list(aic_table = aic_table, best_model = best_model, data = df, timeVar = timeVar, eventVar = eventVar, treatmentVar = treatmentVar)
  })
  output$kmPlot <- renderPlot({
    analysisResults()
    df <- analysisResults()$data
    timeVar <- analysisResults()$timeVar
    eventVar <- analysisResults()$eventVar
    treatmentVar <- analysisResults()$treatmentVar
    
    # Construct the survival object directly
    surv_object <- Surv(time = df[[timeVar]], event = df[[eventVar]])
    
    # Generate survival fit for the plot
    surv_fit <- survfit(surv_object ~ df[[treatmentVar]], data = df)
    
    # Ensure strata exists even with one group
    if (length(surv_fit$strata) == 0) {
      surv_fit$strata <- as.factor(rep("All", nrow(df)))
    }
    
    # Calculate p-value using survdiff
    surv_diff <- survdiff(surv_object ~ df[[treatmentVar]], data = df)
    pvalue <- pchisq(surv_diff$chisq, length(surv_diff$n) - 1, lower.tail = FALSE)
    
    # Manually create plot data for single or multi group
    plot_data <- data.frame(
      time = surv_fit$time,
      n.risk = surv_fit$n.risk,
      n.event = surv_fit$n.event,
      n.censor = surv_fit$n.censor,
      surv = surv_fit$surv,
      strata = rep(names(surv_fit$strata), surv_fit$strata) # use strata names if multiple groups
    )
    # If there's only one group, use the column name from the formula
    if (length(surv_fit$strata) == 1) {
      plot_data$strata <- df[[treatmentVar]]
    }
    
    # Create the plot using ggplot2
    p <- ggplot(plot_data, aes(time, surv, color = strata)) +
      geom_step() +
      geom_point() +
      labs(x = timeVar, y = "Survival Probability", color = treatmentVar) +
      ggtitle("Kaplan-Meier Survival Curve") +
      xlim(0, 50) +
      theme_bw() +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent)
    
    # Facet if checkbox is checked
    if (input$facetByTreatment) {
      p <- p + facet_wrap(~strata)
    }
    
    p
  })
  
  output$summaryTable <- renderTable({
    req(analysisResults())
    aic_table <- analysisResults()$aic_table
    aic_table
  })
}

shinyApp(ui = ui, server = server)