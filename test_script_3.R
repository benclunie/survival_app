library(shiny)  # Load the shiny library for building web applications in R
library(survival)  # Load the survival library for survival analysis
library(survminer)  # Load the survminer library for visualizing survival analysis
library(AICcmodavg)  # Load the AICcmodavg library for model selection
library(tidyr)  # Load the tidyr library for data tidying
library(ggplot2)  # Load the ggplot2 library for data visualization
library(dplyr)  # Load the dplyr library for data manipulation
library(htmlTable)  # Load the htmlTable library for creating HTML tables

# Define the user interface (UI) of the Shiny application
ui <- fluidPage(
  includeCSS("custom.css"),  # Include custom CSS for styling
  tags$head(
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'https://fonts.googleapis.com/css?family=Montserrat&display=swap')  # Include Google font
  ),
  titlePanel(
    title = tags$div(
      class = "title-panel",
      tags$div(
        class = "main-container",
        tags$div(
          class = "left-side",
          tags$div(
            tags$span("survivoR", class = "left-side-title") 
          ),
          tags$div(
            tags$span("An App for Modelling Survival Analysis for Insect Bioassays", class = "left-side-subtitle") 
          )
        ),
        tags$div(
          class = "right-side",
          tags$a(
            href = 'https://www.harper-adams.ac.uk/',  # Link to Harper Adams University
            tags$img(
              src = 'HAU.png',  # Harper Adams logo
              height = 201,
              width = 148
            ),
            target = '_blank',
            class = 'right-side-logo'
          ),
          tags$div(
            style = "display:inline-block; width:5px; background-color:white; height:148px; margin-left:20px; margin-right:20px;"  # White line separator
          ),
          tags$a(
            href = 'https://linktr.ee/hau_entomology',  # Link to Entomology
            tags$img(
              src = 'ento.png',  # Entomology logo
              height = 188,
              width = 188
            ),
            target = '_blank',
            class = 'right-side-logo'
          )
        )
      )
    ),
    windowTitle = 'Survival Analysis for Insect Bioassays'  # Window title for the browser tab
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File", accept = ".csv"),  # File input for CSV file
      selectInput("treatment", "Treatment Column", choices = NULL),  # Drop down for selecting treatment column
      div(class = "dose-container",
          div(class = "three-quarter-width",
              selectInput("dose", "Dose Column", choices = NULL)  # Drop down for selecting dose column
          ),
          div(class = "one-quarter-width",
              textInput("doseUnit", "Dose Unit", value = "")  # Text input for dose unit
          )
      ),
      selectInput("time", "Time Column", choices = NULL),  # Drop down for selecting time column
      selectInput("status", "Status Column", choices = NULL)  # Drop down for selecting status column
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Survival Curves",
                 fluidRow(
                   column(12, uiOutput("plotUI"))  # UI output for the plot
                 )
        ),
        tabPanel("Model Selection", htmlOutput("modelSelection"),  # UI output for model selection
                 uiOutput("aicText"))  # UI output for AIC text explanation
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  data <- reactive({
    req(input$file)  # Require the file input to be available
    
    # Read the input file
    data_input <- tryCatch(
      {
        read.csv(input$file$datapath, stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = c("", "NA"))  # Read CSV file
      },
      error = function(e) {
        stop("Error reading the input file. Please ensure the file is a valid CSV file with appropriate encoding and formatting.")  # Error handling
      }
    )
    
    # Debugging: Print the structure of the data
    print(str(data_input))
    
    # Check if the input data is a data frame
    if (!is.data.frame(data_input)) {
      stop("The input file does not contain a valid data frame. Please check the file format and structure.")  # Error handling
    }
    
    # Check if the input data has any rows
    if (nrow(data_input) == 0) {
      stop("The input file is empty. Please provide a non-empty CSV file.")  # Error handling
    }
    
    # Check if the input data has any columns
    if (ncol(data_input) == 0) {
      stop("The input file does not contain any columns. Please provide a valid CSV file with at least one column.")  # Error handling
    }
    
    # Return the input data
    data_input
  })
  
  # Observe changes in data and update select inputs accordingly
  observe({
    req(data())
    updateSelectInput(session, "treatment", choices = names(data()))  # Update treatment select input
    updateSelectInput(session, "dose", choices = names(data()))  # Update dose select input
    updateSelectInput(session, "time", choices = names(data()))  # Update time select input
    updateSelectInput(session, "status", choices = names(data()))  # Update status select input
  })
  
  # Reactive function to prepare survival data
  surv_data <- reactive({
    req(input$file, input$treatment, input$dose, input$time, input$status)
    
    # Read the input file
    data_input <- tryCatch(
      {
        read.csv(input$file$datapath, stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = c("", "NA"))
      },
      error = function(e) {
        stop("Error reading the input file. Please ensure the file is a valid CSV file with appropriate encoding and formatting.")
      }
    )
    
    # Prepare the data
    colnames(data_input) <- make.names(colnames(data_input))
    data_input$Group <- with(data_input, interaction(get(input$treatment), get(input$dose)))
    data_input$Dose <- factor(data_input[[input$dose]], levels = unique(data_input[[input$dose]]))
    data_input$Treatment <- factor(data_input[[input$treatment]], levels = unique(data_input[[input$treatment]]))
    data_input[[input$time]] <- as.numeric(data_input[[input$time]])
    data_input[[input$status]] <- as.numeric(data_input[[input$status]])
    
    data_input
  })
  
  # Reactive expression for survival model fitting
  surv_model <- reactive({
    df <- surv_data()
    req(df)
    
    # Debugging: Print the structure of the data frame
    print(str(df))
    
    # Fit survival model
    fit <- survfit(Surv(df[[input$time]], df[[input$status]]) ~ Treatment + Dose, data = df)
    
    # Debugging: Print the summary of the fitted model
    print(fit)
    
    fit
  })
  
  # Render survival plot
  output$survPlot <- renderPlot({
    req(input$file, input$treatment, input$dose, input$time, input$status, input$plotButton)  # Require all necessary inputs
    
    df <- surv_data()
    
    # Fit survival model
    fit <- survfit(Surv(df[[input$time]], df[[input$status]]) ~ Treatment + Dose, data = df)
    
    # Debugging: Print the summary of the fitted model
    print(summary(fit))
    
    # Create base survival plot
    plot_surv <- ggsurvplot(
      fit,
      data = df,  # Pass the data frame directly
      fun = "event",
      legend.title = "Dose",
      conf.int = FALSE,
      palette = "jama",
      pval = FALSE
    )
    
    # Make the plot
    plot_surv$plot +
      geom_line(data = df, aes(x = !!as.name(input$time), y = !!as.name(input$status), group = Group, color = Dose), size = 0.75, alpha = 0.8) +
      xlab("Time (hours)") +
      facet_grid(~ Treatment) +
      theme_bw()
  }, width = reactive(input$plotWidth), height = reactive(input$plotHeight))
  
  observeEvent(input$plotButton, {
    traceback()  # Print traceback for debugging
  })
  
  # Plot outputs with positioning for resizing and generate button
  output$plotUI <- renderUI({
    fluidRow(
      column(12,
             div(class = "plot-container",
                 plotOutput("survPlot", width = "100%", height = "auto"),  # Plot output for survival plot
                 div(class = "plot-controls",
                     fluidRow(
                       column(4, numericInput("plotWidth", "Plot Width", value = 700, min = 300)),  # Numeric input for plot width
                       column(4, numericInput("plotHeight", "Plot Height", value = 500, min = 300)),  # Numeric input for plot height
                       column(4, br(), actionButton("plotButton", "Generate Plot"))  # Button to generate plot
                     )
                 )
             )
      )
    )
  })
  
  # Function to round numeric columns of a data frame
  round_df <- function(df, digits) {
    df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x)  # Round numeric columns
    df
  }
  
  # Model selection table output 
  output$modelSelection <- renderUI({
    surv <- surv_data()
    req(surv)
    
    surv_T_x_D <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Treatment * Dose, data = surv)  # Model with Treatment and Dose interaction
    surv_T_D <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Treatment + Dose, data = surv)  # Model with Treatment and Dose additive
    surv_T <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Treatment, data = surv)  # Model with Treatment only
    surv_D <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Dose, data = surv)  # Model with Dose only
    
    models <- list(surv_T_x_D, surv_T_D, surv_T, surv_D)  # List of models
    model_names <- c("Treatment x Dose", "Treatment + Dose", "Treatment", "Dose")  # Model names
    
    aic_table <- aictab(models, modnames = model_names)  # Create AIC table
    aic_table <- round_df(aic_table, 2)  # Round values to 2 decimal places
    
    HTML(htmlTable::htmlTable(aic_table,  # Render AIC table as HTML
                              css.class = "model-selection-table",
                              caption = "AIC Model Selection Table"))
  })
  
  output$aicText <- renderUI({
    HTML("
      <p><strong>K</strong> = number of estimated parameters</p>
      <p><strong>&#916;AIC</strong> = relative differences between selected model and subsequent AIC values</p>
      <p><strong>AIC wt</strong> = relative likelihood weighting for the relevant model</p>
      <p><strong>Cum wt</strong> = cumulative AIC wt values</p>
      <p><strong>LL</strong> = log likelihood (maximum likelihood estimation)</p>
    ")  # HTML explanation of AIC table columns
  })
}

# Run the application
shinyApp(ui = ui, server = server)

