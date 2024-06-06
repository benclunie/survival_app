library(shiny)
library(survival)
library(survminer)
library(AICcmodavg)
library(tidyr)
library(ggplot2)
library(dplyr)
library(htmlTable)

ui <- fluidPage(
  # Include custom CSS for styling
  includeCSS("custom.css"),
  # Header elements
  tags$head(
    # Include Google font
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'https://fonts.googleapis.com/css?family=Montserrat&display=swap')
  ),
  # Title panel
  titlePanel(
    # Title panel content
    title = tags$div(
      class = "title-panel",
      # Main container
      tags$div(
        class = "main-container",
        # Left side: title, version, and subtitle
        tags$div(
          class = "left-side",
          tags$div(
            tags$span("survivoR", class = "left-side-title")
          ),
          tags$div(
            tags$span("An App for Modelling Survival Analysis for Insect Bioassays", class = "left-side-subtitle")
          )
        ),
        # Right side: logos
        tags$div(
          class = "right-side",
          # Harper Adams logo with padding
          tags$a(
            href = 'https://www.harper-adams.ac.uk/',
            tags$img(
              src = 'HAU.png',  
              height = 201,
              width = 148
            ),
            target = '_blank',
            class = 'right-side-logo'
          ),
          # White line separator
          tags$div(
            style = "display:inline-block; width:5px; background-color:white; height:148px; margin-left:20px; margin-right:20px;"
          ),
          # Entomology logo with padding
          tags$a(
            href = 'https://linktr.ee/hau_entomology',
            tags$img(
              src = 'ento.png',  # Reference to the image in the www directory
              height = 188,
              width = 188
            ),
            target = '_blank',
            class = 'right-side-logo'
          )
        )
      )
    ),
    # Window title for the browser tab
    windowTitle = 'Survival Analysis for Insect Bioassays'
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File", accept = ".csv"),
      selectInput("treatment", "Treatment Column", choices = NULL),
      div(class = "dose-container",
          div(class = "half-width",
              selectInput("dose", "Dose Column", choices = NULL)
          ),
          div(class = "half-width",
              textInput("doseUnit", "Dose Unit", value = "")
          )
      ),
      selectInput("time", "Time Column", choices = NULL),
      selectInput("status", "Status Column", choices = NULL)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Survival Curves",
                 fluidRow(
                   column(12, uiOutput("plotUI"))
                 )
        ),
        tabPanel("Model Selection", htmlOutput("modelSelection"),
                 uiOutput("aicText"))
      )
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    infile <- read.csv(input$file$datapath)
    infile
  })
  
  observe({
    req(data())
    updateSelectInput(session, "treatment", choices = names(data()))
    updateSelectInput(session, "dose", choices = names(data()))
    updateSelectInput(session, "time", choices = names(data()))
    updateSelectInput(session, "status", choices = names(data()))
  })
  
  surv_data <- reactive({
    if (is.null(input$treatment) || is.null(input$dose) || is.null(input$time) || is.null(input$status)) {
      return(NULL)
    }
    
    surv <- data()
    colnames(surv) <- make.names(colnames(surv))
    surv$Group <- with(surv, interaction(get(input$treatment), get(input$dose)))
    surv$Dose <- factor(surv[[input$dose]], levels = unique(surv[[input$dose]]))
    surv$Treatment <- factor(surv[[input$treatment]], levels = unique(surv[[input$treatment]]))
    surv <- surv[surv$Treatment != "Untreated", ]
    surv[[input$time]] <- as.numeric(surv[[input$time]])
    surv[[input$status]] <- as.numeric(surv[[input$status]])
    surv
  })
  
  observeEvent(input$plotButton, {
    output$survPlot <- renderPlot({
      surv_data_val <- surv_data()
      if (!is.null(surv_data_val)) {
        surv_fit <- survfit(Surv(surv_data_val[[input$time]], surv_data_val[[input$status]]) ~
                              Treatment + Dose, data = surv_data_val)
        
        # Extract treatment and dose from the strata names
        surv_fit_df <- data.frame(
          time = surv_fit$time,
          surv = surv_fit$surv,
          strata = rep(names(surv_fit$strata), surv_fit$strata)
        )
        
        # Separate treatment and dose in strata column
        surv_fit_df <- tidyr::separate(surv_fit_df, strata, into = c("Treatment", "Dose"), sep = ", ")
        
        # Apply custom labels if provided
        dose_unit <- if (input$doseUnit != "") paste0(" (", input$doseUnit, ")") else ""
        dose_label <- paste0("Dose", dose_unit)
        
        ggsurvplot(surv_fit, data = surv_data_val, pval = TRUE, conf.int = TRUE, 
                   risk.table = TRUE, 
                   ggtheme = theme_minimal()) +
          labs(x = "Time", y = "Survival Probability", col = dose_label) +
          facet_wrap(~ Treatment + Dose) +
          theme(strip.text = element_text(size = 14)) +
          theme(plot.title = element_text(size = 14, face = "bold"))
      }
    }, width = input$plotWidth, height = input$plotHeight)
  })
  
  output$plotUI <- renderUI({
    fluidRow(
      column(12, 
             div(class = "plot-container",
                 plotOutput("survPlot", width = "100%", height = "auto"),
                 div(class = "plot-controls",
                     fluidRow(
                       column(4, numericInput("plotWidth", "Plot Width", value = 700, min = 300)),
                       column(4, numericInput("plotHeight", "Plot Height", value = 500, min = 300)),
                       column(4, br(), actionButton("plotButton", "Generate Plot"))
                     )
                 )
             )
      )
    )
  })
  
  output$aicText <- renderUI({
    HTML("
      <p><strong>K</strong> = number of estimated parameters</p>
      <p><strong>&#916;AIC</strong> = relative differences between selected model and subsequent AIC values</p>
      <p><strong>AIC wt</strong> = relative likelihood weighting for the relevant model</p>
      <p><strong>Cum wt</strong> = cumulative AIC wt values</p>
      <p><strong>LL</strong> = log likelihood (maximum likelihood estimation)</p>
    ")
  })
}

shinyApp(ui = ui, server = server)
