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
          div(class = "three-quarter-width",
              selectInput("dose", "Dose Column", choices = NULL)
          ),
          div(class = "one-quarter-width",
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

# Define the server logic
server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    
    # Read the input file
    data_input <- tryCatch(
      {
        read.csv(input$file$datapath, stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = c("", "NA"))
      },
      error = function(e) {
        stop("Error reading the input file. Please ensure the file is a valid CSV file with appropriate encoding and formatting.")
      }
    )
    
    # Debugging: Print the structure of the data
    print(str(data_input))
    
    # Check if the input data is a data frame
    if (!is.data.frame(data_input)) {
      stop("The input file does not contain a valid data frame. Please check the file format and structure.")
    }
    
    # Check if the input data has any rows
    if (nrow(data_input) == 0) {
      stop("The input file is empty. Please provide a non-empty CSV file.")
    }
    
    # Check if the input data has any columns
    if (ncol(data_input) == 0) {
      stop("The input file does not contain any columns. Please provide a valid CSV file with at least one column.")
    }
    
    # Return the input data
    data_input
  })
  
  observe({
    req(data())
    updateSelectInput(session, "treatment", choices = names(data()))
    updateSelectInput(session, "dose", choices = names(data()))
    updateSelectInput(session, "time", choices = names(data()))
    updateSelectInput(session, "status", choices = names(data()))
  })
  
  surv_data <- reactive({
    req(input$treatment, input$dose, input$time, input$status)
    df <- data()
    
    # Data Preparation and tidying
    colnames(df) <- make.names(colnames(df))  # Make column names syntactically valid
    df$Group <- with(df, interaction(get(input$treatment), get(input$dose)))
    df$Dose <- factor(df[[input$dose]], levels = unique(df[[input$dose]]))
    df$Treatment <- factor(df[[input$treatment]], levels = unique(df[[input$treatment]]))
    df[[input$time]] <- as.numeric(df[[input$time]])
    df[[input$status]] <- as.numeric(df[[input$status]])
    
    # Debugging: Print the structure of the prepared data
    print(str(df))
    
    # Explicitly convert to data.frame
    df <- as.data.frame(df)
    
    # Additional debugging: print head of the dataframe
    print(head(df))
    
    df
  })
  
  plot_dimensions <- reactive({
    list(width = input$plotWidth, height = input$plotHeight)
  })
  
  output$survPlot <- renderPlot({
    req(input$plotButton)  
    df <- surv_data()
    req(df)
    
    dimensions <- plot_dimensions()
    
    # Data tidying
    df$Treatment <- factor(df$Treatment, levels = sort(unique(df$Treatment)))
    df$Dose <- factor(df$Dose)
    df$Group <- with(df, interaction(Treatment, Dose))
    
    surv_object <- Surv(df[[input$time]], df[[input$status]]) 
    
    # Debugging: Print the structure of the survival object
    print(str(surv_object))
    
    # Simplified ggsurvplot call to isolate the problem
    fit <- survfit(surv_object ~ Treatment + Dose, data = df)
    print(fit)  # Debugging: Print the fit object
    
    # Create base survival plot 
    plot_surv <- ggsurvplot(
      fit,
      data = df,  # Add the data argument
      fun = "event",  # Function to apply to the survival curve (e.g., "event", "censor", "counting", "cumhaz")
      legend.title = "Dose",  # Title for the legend
      conf.int = FALSE,
      palette = "jama",
      pval = FALSE
    )
    
    # Make the plot
    plot_surv$plot +
      geom_line(data = df, aes(x = Time, y = Status, group = Group, color = Dose), size = 0.75, alpha = 0.8) +
      xlab("Time (hours)") +
      facet_grid(~ Treatment) +
      theme_bw()
  }, width = reactive(input$plotWidth), height = reactive(input$plotHeight))
  
  observeEvent(input$plotButton, {
    traceback()  
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
  
  output$modelSelection <- renderUI({
    surv <- surv_data()
    req(surv)
    
    surv_T_x_D <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Treatment * Dose, data = surv)
    surv_T_D <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Treatment + Dose, data = surv)
    surv_T <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Treatment, data = surv)
    surv_D <- survreg(Surv(surv[[input$time]], surv[[input$status]]) ~ Dose, data = surv)
    
    models <- list(surv_T_x_D, surv_T_D, surv_T, surv_D)
    model_names <- c("Treatment x Dose", "Treatment + Dose", "Treatment", "Dose")
    
    aic_table <- aictab(models, modnames = model_names)
    aic_table <- round_df(aic_table, 2)  # Round values to 2 decimal places
    
    HTML(htmlTable::htmlTable(aic_table,
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
    ")
  })
  
  round_df <- function(df, digits) {
    df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x)
    df
  }
}

shinyApp(ui = ui, server = server)
