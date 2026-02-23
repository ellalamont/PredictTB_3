# 1/12/26
# 2/20/26: Adjusted the p-value and log2fold change options
# 2/20/26: Making things complicated by adding a summary table so we see sample sizes

library(shiny)
library(gmodels)
library(DT)

source("Import_DEG_sets.R") # for list_dfs_f - All FDR-adjusted p-values and with Log2Fold >1 or >2 (Rv genes only)
source("Import_GeneSets.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=10))



# Define UI ----
ui <- fluidPage(
  titlePanel("Run1-4 >60% Txn Coverage Volcanos"),
  
  fluidRow(
    
    column(width = 4,
           selectInput("my_comparison",
                       label = "Volcano plots",
                       choices = df_names,
                       width = "100%"),
           
           # Add buttons for changing the log2fold threshold
           radioButtons("log2fc_threshold", 
                        label = "log2 fold change threshold", 
                        choices = c("log2FC > abs(1)" = "1", "log2FC > abs(2)" = "2", "log2FC > abs(2.5)" = "3"),
                        selected = "1",
                        inline = TRUE),
          
           fluidRow(
             column(width = 3,
                    textInput("my_GeneID", 
                              label = "Gene ID (Will label gene orange)",
                              value = "Rv...")
             )
           ),
           # Add checkbox to toggle gene set points
           checkboxInput("show_gene_set", label = "Show gene sets", value = FALSE),
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           # Dropdown for selecting the gene set within the chosen rda file.
           # Start with no selection
           selectInput("my_GeneSet",
                       label = "Gene Set (Will label genes yellow)",
                       choices = NULL),
           
           h5("Sample Size Filters"),
           selectInput("filter_week",
                       label = "Week",
                       choices = unique(my_pipeSummary$Week),
                       selected = NULL, multiple = T),
           selectInput("filter_outcome",
                       label = "Outcome",
                       choices = unique(my_pipeSummary$Outcome),
                       selected = NULL, multiple = T),
           checkboxGroupInput("breakdown_vars", 
                              label = "Column Options",
                              choices = c("Arm" = "Arm", "Lineage" = "main_lineage"),
                              selected = NULL)
    ),
    
    column(width = 8,
           plotlyOutput("volcano_plot",
                        width = "90%", height = "600px"),
           
           hr(),
           h4("Sample Counts"),
           DT::DTOutput("sample_summary_table")
    ),
  )
  
)

# Define server logic ----
server <- function(input, output, session) {
  
  # When a new gene set source is selected, update the gene set dropdown
  observeEvent(input$my_GeneSetSource, {
    updateSelectInput(session, "my_GeneSet",
                      choices = names(allGeneSetList[[input$my_GeneSetSource]]),
                      selected = NULL)
  })
  
  # Volcano Plot
  output$volcano_plot <- renderPlotly({
    
    # Add data for labelling a single gene
    single_gene <- list_dfs_f2[[input$my_comparison]] %>% 
      filter(GENE_ID == input$my_GeneID)
    
    # Add data for labelling a gene set
    gene_set <- list_dfs_f2[[input$my_comparison]] %>%
      filter(GENE_ID %in% allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]])
    
    # Make new data name for plotting
    plot_data <- list_dfs_f2[[input$my_comparison]]
    
    # Choose DE column and DE_labels column based on selected threshold
    if (input$log2fc_threshold == "1") {
      de_col <- "DE1"
      label_col <- "DE1_labels"
      p_value <- "AVG_PVALUE"
      log2fc_cutoff <- 1
    } else if (input$log2fc_threshold == "2") {
      de_col <- "DE2"
      label_col <- "DE2_labels"
      p_value <- "AVG_PVALUE"
      log2fc_cutoff <- 2
    } else if (input$log2fc_threshold == "3") {
      de_col <- "DE2.5"
      label_col <- "DE2.5_labels"
      p_value <- "AVG_PVALUE"
      log2fc_cutoff <- 2.5
    }
    
    
    # Make the Volcano Plot
    my_volcano <- plot_data %>%
      
      ggplot(aes(x = LOG2FOLD, y = -log10(.data[[p_value]]), 
                 col = .data[[de_col]], label = .data[[label_col]], 
                 text = GENE_ID, label2 = GENE_NAME, label3 = PRODUCT)) + 
      geom_point(alpha = 0.7) + 
      
      # Add a differently colored point
      geom_point() +
      # Conditionally add the gene set points (yellow) based on checkbox value
      { if(input$show_gene_set) 
        geom_point(data = gene_set, color = "yellow", aes(col = .data[[de_col]], label = .data[[label_col]], text = GENE_ID))
        else 
          NULL } +
      geom_point(data = single_gene, color = "yellow", aes(col = .data[[de_col]], label = .data[[label_col]], text = GENE_ID)) + 
      
      labs(title = input$my_comparison) + 
      geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), col = "grey", linetype = "dashed") + 
      geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
      scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
    # geom_label_repel(max.overlaps = 10) # Can do geom_text_repel or geom_label_rebel
    
    # Determine the max and min axes values for labeling 
    plot_build <- ggplot_build(my_volcano)
    y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
    x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
    x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
    
    # Add the gene number annotations
    text_up <- list_dfs_f2[[input$my_comparison]] %>% filter(.data[[de_col]] == "significant up") %>% nrow()
    text_down <- list_dfs_f2[[input$my_comparison]] %>% filter(.data[[de_col]] == "significant down") %>% nrow()
    my_volcano_annotated <- my_volcano +
      annotate("text", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
      annotate("text", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
    
    final_plot <- my_volcano_annotated + my_plot_themes 
    final_plot
  })
  
  # Reactive Table of sample sizes
  summary_table <- reactive({
    
    df <- my_pipeSummary %>%
      filter(!is.na(Patient)) %>%
      filter(Txn_Coverage_f >= 60)
    
    # Apply filters only if selected
    if (!is.null(input$filter_week))
      df <- df %>% filter(Week %in% input$filter_week)
    if (!is.null(input$filter_outcome))
      df <- df %>% filter(Outcome %in% input$filter_outcome)
    
    # Group only if a filter is provided for Arm and Lineage
    grouping_vars <- c("Week", "Outcome")
    if (length(input$breakdown_vars) > 0)
      grouping_vars <- c(grouping_vars, input$breakdown_vars)
    
    df %>%
      group_by(across(all_of(grouping_vars))) %>%
      summarize(N_samples = n(), .groups = "drop")
  })
  
  output$sample_summary_table <- DT::renderDT({
    DT::datatable(
      summary_table(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)