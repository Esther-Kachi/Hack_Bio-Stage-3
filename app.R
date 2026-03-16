#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)



library(shiny)
library(dplyr)
library(ggplot2)

# BOILERPLATE FUNCTIONS 

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], 
               max = 255, 
               alpha = (100 - percent) * 255 / 100, 
               names = name)
  invisible(t.col)
}

scale_0_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(50, length(x)))
  (x - rng[1]) / (rng[2] - rng[1]) * 100
}

compute_gene_stats <- function(expr_mat, meta_df, target_ct) {
  in_cells <- meta_df$cell_id[meta_df$cell_type == target_ct]
  out_cells <- meta_df$cell_id[meta_df$cell_type != target_ct]
  
  xin <- expr_mat[in_cells, , drop = FALSE]
  xout <- expr_mat[out_cells, , drop = FALSE]
  
  det_in <- colMeans(xin > 0)
  det_out <- colMeans(xout > 0)
  mean_in <- colMeans(xin)
  mean_out <- colMeans(xout)
  diff <- mean_in - mean_out
  
  return(data.frame(
    gene = colnames(expr_mat),
    det_in = det_in,
    det_out = det_out,
    mean_in = mean_in,
    mean_out = mean_out,
    diff = diff,
    stringsAsFactors = FALSE
  ))
}

pick_marker_gene <- function(gene_stats_df) {
  gene_stats_df <- gene_stats_df %>% 
    arrange(desc(diff), desc(det_in))
  gene_stats_df$gene[1]
}

# DATA LOADING & ALIGNMENT 

# Load datasets
expr_mat <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/expression_matrix.csv", row.names = 1)
meta_df <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/cell_metadata.csv")
umap_df <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/umap_coordinates.csv")

# Ensure alignment of Cell IDs
common_cells <- intersect(rownames(expr_mat), intersect(meta_df$cell_id, umap_df$cell_id))
expr_mat <- expr_mat[common_cells, ]
meta_df <- meta_df[meta_df$cell_id %in% common_cells, ]
meta_df <- meta_df[match(common_cells, meta_df$cell_id), ]
umap_df <- umap_df[umap_df$cell_id %in% common_cells, ]
umap_df <- umap_df[match(common_cells, umap_df$cell_id), ]

# Define user interface
ui <- fluidPage(
  titlePanel("Cell Type Snapshot Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("cell_type", "Select Cell Type/Cluster:", 
                  choices = unique(meta_df$cell_type)),
      hr(),
      wellPanel(
        h4("Marker Gene Info"),
        verbatimTextOutput("marker_info")
      )
    ),
    
    mainPanel(
      plotOutput("umap_plot", height = "500px"),
      hr(),
      h4("Gene Statistics Table"),
      tableOutput("stats_table")
    )
  )
)

# Define the Server
server <- function(input, output) {
  
  # Reactive expression to compute stats when cell type changes
  stats_reactive <- reactive({
    compute_gene_stats(expr_mat, meta_df, input$cell_type)
  })
  
  # Reactive expression to pick the best marker
  marker_reactive <- reactive({
    pick_marker_gene(stats_reactive())
  })
  
  # Output: Marker Text Info
  output$marker_info <- renderText({
    best_gene <- marker_reactive()
    stats <- stats_reactive() %>% filter(gene == best_gene)
    paste0("Selected Cell Type: ", input$cell_type, "\n",
           "Best Marker Gene: ", best_gene, "\n",
           "Specificity Score (diff): ", round(stats$diff, 4))
  })
  
  # Output: UMAP Plot
  output$umap_plot <- renderPlot({
    best_gene <- marker_reactive()
    # Get expression for the best gene and scale
    gene_expr <- expr_mat[, best_gene]
    scaled_expr <- scale_0_100(gene_expr)
    
    plot_data <- cbind(umap_df, Expression = scaled_expr)
    
    ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Expression)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_gradient(low = "lightgrey", high = "blue") +
      theme_minimal() +
      labs(title = paste("Marker Signal:", best_gene),
           subtitle = paste("Targeted Cluster:", input$cell_type),
           color = "Scaled Expr (0-100)")
  })
  
  # Output: Stats Table
  output$stats_table <- renderTable({
    stats_reactive() %>%
      arrange(desc(diff)) %>%
      head(10) # Showing top 10 for clarity
  })
}

# Run the application
shinyApp(ui = ui, server = server)
