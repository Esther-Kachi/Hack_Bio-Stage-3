# Hack_Bio-Stage-3

BIO-DATA-VISUALIZATION-R-STAGE-THREE
Cell Type Snapshot Explorer (scRNA‑seq) AUTHOR: ESTHER AZUNNA A lightweight Shiny app for exploring UMAP embeddings, marker gene expression, and gene specificity scores across cell types.

How to Run the App

Install required packages:install.packages(c("shiny", "ggplot2", "dplyr", "readr"))
Run directly from the project folder: shiny::runApp()
Or, launch the deployed app via the shinyapps.io link: 

Gene Specificity Score Each gene is ranked by a simple specificity score: diff = mean_in − mean_out Where: mean_in: average expression of the gene within the selected cell type mean_out: average expression outside that cell type A larger diff means the gene is more specific (highly expressed inside the group and low outside).

Marker Gene Selection The app automatically identifies a marker gene for the selected cell type by: 1.Computing diff for every gene 2.Ranking genes by descending diff 3.Selecting the top gene (highest specificity score) as the marker This marker gene is then used to color the UMAP plot and summarize expression patterns.

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
