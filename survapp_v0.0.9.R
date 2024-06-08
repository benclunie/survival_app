library(shiny)
library(survival)
library(survminer)
library(AICcmodavg)
library(tidyr)
library(ggplot2)
library(dplyr)
library(htmlTable)

ui <- fluidPage(
  includeCSS("custom.css"),
  tags$head(
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'https://fonts.googleapis.com/css?family=Montserrat&display=swap')
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
            href = 'https://www.harper-adams.ac.uk/',
            tags$img(
              src = 'HAU.png', 
              height = 201,
              width = 148
            ),
            target = '_blank',
            class = 'right-side-logo'
          ),
          tags$div(
            style = "display:inline-block; width:5px; background-color:white; height:148px; margin-left:20px; margin-right:20px;"
          ),
          tags$a(
            href = 'https://linktr.ee/hau_entomology',
            tags$img(
              src = 'ento.png',
              height = 188,
              width = 188
            ),
            target = '_blank',
            class = 'right-side-logo'
          )
        )
      )
    ),
    windowTitle = 'Survival Analysis for Insect Bioassays'
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose CSV File", accept = ".csv"),
      selectInput("model_select", "Select Model", choices = NULL)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions",
                 fluidRow(
                   column(12, uiOutput("instructionsUI"))
                 )
        ),
        tabPanel("Survival Curves",
                 fluidRow(
                   column(12, uiOutput("plotUI"))
                 )
        ),
        tabPanel("Model Selection", htmlOutput("modelSelection"), 
                 uiOutput("aicText"))
      )
    )
  ),
  tags$footer(
    class = "footer",
    HTML("&copy; Ben Clunie and Joe Roberts 2024")
  )
)  

server <- function(input, output, session) {
  surv_data <- reactive({
    req(input$file)
    data <- read.csv(input$file$datapath, stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = c("", "NA"))
    data$Group <- paste(data$Treatment, data$Dose, sep = "_")
    return(data)
  })
  
  observe({
    req(surv_data())
    updateSelectInput(session, "model_select", choices = c("Treatment x Dose", "Treatment + Dose", "Treatment", "Dose"))
  })
  
  selected_model <- reactive({
    req(surv_data())
    surv <- surv_data()
    surv_T_x_D <- survreg(Surv(Time, Status) ~ Treatment * Dose, data = surv)
    surv_T_D <- survreg(Surv(Time, Status) ~ Treatment + Dose, data = surv)
    surv_T <- survreg(Surv(Time, Status) ~ Treatment, data = surv)
    surv_D <- survreg(Surv(Time, Status) ~ Dose, data = surv)
    
    models <- list(surv_T_x_D, surv_T_D, surv_T, surv_D)
    names(models) <- c("Treatment x Dose", "Treatment + Dose", "Treatment", "Dose")
    
    models[[input$model_select]]
  })
  
  output$survivalPlot <- renderPlot({
    surv <- surv_data()  # Ensure the survival data is scoped here
    model <- selected_model()
    
    pred <- data.frame(predict(model, type = "quantile", p = seq(0.01, 0.99, by = 0.01)))
    pred$Group <- surv$Group
    pred$Treatment <- surv$Treatment
    pred$Dose <- surv$Dose
    
    pred2 <- array(dim = c(nlevels(as.factor(surv$Treatment)) * nlevels(as.factor(surv$Dose)), ncol(pred)))
    for (i in 1:ncol(pred2)) {
      pred2[, i] <- tapply(pred[, i], pred$Group, unique)
    }
    
    y_val <- seq(0.99, 0.01, by = -0.01)
    pred2 <- data.frame(t(pred2[, 1:length(y_val)]), y_val)
    pred2_long <- gather(pred2, key = "Group", value = "time", -y_val)
    pred2_long$time <- as.numeric(pred2_long$time)
    pred2_long$Group <- rep(levels(as.factor(surv$Group)), each = length(y_val))
    pred2_long$Treatment <- rep(rep(levels(as.factor(surv$Treatment)), each = length(y_val)), nlevels(as.factor(surv$Dose)))
    pred2_long$Dose <- rep(rep(levels(as.factor(surv$Dose)), each = length(y_val) * nlevels(as.factor(surv$Treatment))))
    pred2_long <- na.omit(pred2_long)
    
    surv_fit <- survfit(Surv(Time, Status) ~ Treatment + Dose, data = surv)
    plot_surv <- ggsurvplot(fit = surv_fit, data = surv, conf.int = FALSE, col = "Dose")
    
    plot_surv$plot +
      geom_line(data = pred2_long, aes(x = time, y = y_val, group = Group, colour = Dose), size = 0.75, alpha = 0.8) +
      xlab("Time (hours)") +
      facet_grid(~Treatment) +
      theme_bw()
  }, width = reactive(input$plotWidth), height = reactive(input$plotHeight))
  
  output$plotUI <- renderUI({
    fluidRow(
      column(12,
             div(class = "plot-container",
                 plotOutput("survivalPlot", width = "100%", height = "auto"),
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
  
  round_df <- function(df, digits) {
    df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x)
    df
  }
  
  output$modelSelection <- renderUI({
    surv <- surv_data()
    req(surv)
    
    surv_T_x_D <- survreg(Surv(surv$Time, surv$Status) ~ Treatment * Dose, data = surv)
    surv_T_D <- survreg(Surv(surv$Time, surv$Status) ~ Treatment + Dose, data = surv)
    surv_T <- survreg(Surv(surv$Time, surv$Status) ~ Treatment, data = surv)
    surv_D <- survreg(Surv(surv$Time, surv$Status) ~ Dose, data = surv)
    
    models <- list(surv_T_x_D, surv_T_D, surv_T, surv_D)
    model_names <- c("Treatment x Dose", "Treatment + Dose", "Treatment", "Dose")
    
    aic_table <- aictab(models, modnames = model_names)
    aic_table <- round_df(aic_table, 2)
    
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
  
  output$instructionsUI <- renderUI({
    HTML("
      <h3>Instructions</h3>
      <p>Add your instructions here.</p>
    ")
  })
}

shinyApp(ui = ui, server = server)
