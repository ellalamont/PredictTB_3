library(shiny)

# dev.off() # Sometimes need this to get the shiny to show the pheatmap?

library(devtools)
library(DuffyTools)

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# Define UI ----
ui <- fluidPage(

  titlePanel("Forest Plots"),
  
  fluidRow(
    column(width = 3,
           
           # Need this extra thing for printing the console messages
           tags$style(HTML("
           pre {
           white-space: pre-wrap;
           word-break: break-word;
           margin: 0; padding: 0;}")),
           
           # Dropdown for selecting with DEG file to use
           selectInput("my_DEG_file",
                       label = "DEG file",
                       choices = df_names),
           
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           
           # Adjust text size
           numericInput("my_text_size", 
                        label = "Text size", 
                        value = 1.1, min = 0.5, max = 2, step = 0.1),
           
           hr(),  # nice horizontal line
           h6("Server Output"),
           verbatimTextOutput("console_output")  # To show what the server prints
           ),
    
    column(width = 9, # Max is 12...
           plotOutput("forest_plot", height = "650px"))
           # uiOutput("dynamic_forest_plot")) # To make the Forest plot change length based on number of gene sets
    )
  )


# Define server logic ----
server <- function(input, output, session){
  
  # Store console messages
  console_messages <- reactiveVal("")  # To store logs
  
  # To make the Forest plot change length based on number of gene sets
  # output$dynamic_forest_plot <- renderUI({
  #   n_GeneSets <- length(allGeneSetList[[input$my_GeneSetSource]])
  #   base_height <- 200 # minimum height (pixels)
  #   height_per_gene <- 1.5 # how many pixels per gene set
  #   plot_height <- max(base_height, n_GeneSets * height_per_gene)
  #   plotOutput("forest_plot", height = paste0(plot_height, "px"))
  # })
  
  # Make the Forest Plot
  output$forest_plot <- renderPlot({
    
    split_name <- strsplit(input$my_DEG_file, "\\.ComparedTo\\.")[[1]] # Separate the df_name 
    Right_DataSet <- split_name[1]
    Left_DataSet  <- split_name[2]
    
    captured_output <- capture.output({ # Need this to capture the message printed to the console
      plotGeneSetForest(file = list_dfs_f2[[input$my_DEG_file]],
                        geneSets = allGeneSetList[[input$my_GeneSetSource]],
                        main = input$my_DEG_file,
                        # max.show = 30,
                        min.genes.per.set = 3,
                        left.label = paste0(Left_DataSet), 
                        right.label = paste0(Right_DataSet),
                        xRange = 3, # Changes how far out the log2fold change axis goes
                        text.cex = input$my_text_size, pt.cex = 1.75, lwd = 4.2) 
      })
    
    # To save what the server prints
    console_messages(paste(captured_output, collapse = "\n"))
  })
  
  # Output the messages
  output$console_output <- renderText({
    console_messages()
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)


