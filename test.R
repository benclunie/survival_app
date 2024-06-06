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
      textInput("doseVar", "Dose Variable (optional)", "Dose"),
      checkboxInput("facetByTreatment", "Facet by Treatment", value = TRUE),
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
    df
  })
  
  analysisResults <- eventReactive(input$analyze, {
    df <- data()
    timeVar <- input$timeVar
    eventVar <- input$eventVar
    treatmentVar <- input$treatmentVar
    doseVar <- input$doseVar
    
    if (!(timeVar %in% names(df))) stop("Time variable not found in the data")
    if (!(eventVar %in% names(df))) stop("Event variable not found in the data")
    if (!(treatmentVar %in% names(df))) stop("Treatment variable not found in the data")
    if (doseVar != "" && !(doseVar %in% names(df))) stop("Dose variable not found in the data")
    
    df[[timeVar]] <- as.numeric(df[[timeVar]])
    df[[eventVar]] <- as.numeric(df[[eventVar]])
    df[[treatmentVar]] <- as.factor(df[[treatmentVar]])
    if (doseVar != "") {
      df[[doseVar]] <- as.factor(df[[doseVar]])
    }
    
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
    
    aic_table <- aictab(models, modnames = model_names)
    best_model_index <- which.min(aic_table$AICc)
    best_model <- models[[best_model_index]]
    
    list(aic_table = aic_table, best_model = best_model, data = df, timeVar = timeVar, eventVar = eventVar, treatmentVar = treatmentVar, doseVar = doseVar)
  })
  
  output$kmPlot <- renderPlot({
    analysisResults()
    df <- analysisResults()$data
    timeVar <- analysisResults()$timeVar
    eventVar <- analysisResults()$eventVar
    treatmentVar <- analysisResults()$treatmentVar
    doseVar <- analysisResults()$doseVar
    
    surv_object <- Surv(time = df[[timeVar]], event = df[[eventVar]])
    
    if (doseVar != "") {
      surv_fit <- survfit(surv_object ~ df[[treatmentVar]] + df[[doseVar]], data = df)
    } else {
      surv_fit <- survfit(surv_object ~ df[[treatmentVar]], data = df)
    }
    
    plot_data <- data.frame(
      time = surv_fit$time,
      n.risk = surv_fit$n.risk,
      n.event = surv_fit$n.event,
      n.censor = surv_fit$n.censor,
      surv = surv_fit$surv,
      strata = rep(names(surv_fit$strata), surv_fit$strata)
    )
    
    if (length(surv_fit$strata) == 1) {
      plot_data$strata <- df[[treatmentVar]]
    }
    
    p <- ggsurvplot(surv_fit, data = df, aes(color = strata)) +
      labs(x = timeVar, y = "Survival Probability", color = treatmentVar) +
      ggtitle("Kaplan-Meier Survival Curve") +
      theme_bw()
    
    if (input$facetByTreatment && doseVar != "") {
      p <- p + facet_wrap(~ df[[treatmentVar]] + df[[doseVar]])
    } else if (input$facetByTreatment) {
      p <- p + facet_wrap(~ df[[treatmentVar]])
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
