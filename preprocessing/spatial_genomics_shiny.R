require(shiny)
require(bslib)




cluster_page <- tabPanel(
  title = "Cluster",
  titlePanel("if-pan-annaradli-spatial"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("resolution", "Select resolution cluster", value = 0.3, min = 0.05, max = 3)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('UMAP', plotOutput("umap_plot"), plotOutput('umap_sample')),
        tabPanel('Cluster', plotOutput("cluster_plot", width = "800px", height = "600px"), 
                 uiOutput("return_cluster"),
                 plotOutput("interest_cluster_plot", width = "800px", height = "600px"),
                 actionButton("find_markers", "Find markers"),
                 DT::DTOutput("markers_interest_cluster"))
      )
    )
  )
)

feature_page <- tabPanel(
  title = "Feature",
  titlePanel("if-pan-annaradli-spatial"),
    sidebarLayout(
      sidebarPanel(
        selectInput('select_gene', 'Select gene', unique(colfilt.info.peaks$gene_name)),
        uiOutput("return_select_peak"),
        sliderInput("size_spot", "Select size spot parametr", value = 1.3, min = 0.5, max = 2)
      ),
      mainPanel(
        plotOutput("spatial_feature_plot", width = "800px", height = "600px"),
        uiOutput("return_cluster_feature"),
        plotOutput("spatial_feature_plot_cluster", width = "800px", height = "600px"),
        DT::DTOutput("peaks_table_gene")
      )
    )
)


ui <- navbarPage(
  title = "IF-PAN", 
  theme = bs_theme(version = 4, bootswatch = "cosmo"),
  cluster_page, 
  feature_page
)


server <- function(input, output, 
                   session) {
  
  data_cluster <- reactive(FindClusters(integrated.analysis, resolution = input$resolution))
  tmp.bcs_merge <- reactive(data_cluster() %>%.$seurat_clusters %>% 
                              as.data.frame() %>% 
                              rename(cluster = ".") %>%
                              rownames_to_column(var = "sample.barcode") %>%
                              separate("sample.barcode", c("sample", "barcode"), sep = "_") %>% 
                              left_join(., bcs_merge, by = c("barcode", "sample")))
  
  # plot for cluster
  output$cluster_plot <- renderPlot({
    plot_clusters(data_cluster = data_cluster(),
                  size = 1.3)
  })
  
  # select number cluster to plot interest cluster
  output$return_cluster <- renderUI({
    sliderInput('cluster', 'Select cluster',
                value = 1, min = 0, 
                max = data_cluster() %>% levels() %>% 
                  as.numeric() %>% max())
  })
  
  # plot for interest cluster
  output$interest_cluster_plot <- renderPlot({
    plot_interest_cluster(data_cluster = data_cluster(), 
                          interest_cluster = input$cluster,
                          size = 1.3)
  })
  
  # select peak to create plot for interest feature
  output$return_select_peak <- renderUI({
    selectInput("select_peak", "Select peak", colfilt.info.peaks %>% select(gene_name, peak_id)  %>%
                  filter(gene_name == input$select_gene) %>% 
                  select(peak_id) %>%
                  mutate(peak_id = str_replace_all(peak_id, "_", "-")))
  })
  
  # Find markers for cluster
  find_markers <- eventReactive(input$find_markers, {
    FindMarkers(data_cluster(), ident.1 = input$cluster, min.pct = 0.25) %>%
      rownames_to_column(., "peak_id") %>%
      mutate(peak_id = str_replace_all(peak_id, "-", "_")) %>%
      left_join(., colfilt.info.peaks[, c(4, 16)], by = "peak_id") %>% 
      select(gene_name, peak_id, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
  })
  
  # create table markers for interest cluster
  output$markers_interest_cluster <- DT::renderDT({
    find_markers()
  })
  
  # table for interest gene
  output$peaks_table_gene <- DT::renderDT({
    colfilt.info.peaks %>% 
      filter(gene_name == input$select_gene) %>%
      magrittr::set_colnames(colfilt.info.peaks %>% colnames() %>% str_replace_all(.,"_", " ")) %>%
      select(-c('chromosome2', 'start gene', 'end gene',
                'gene id'))
      
  })
  
  # plot for interest feature
  output$spatial_feature_plot <- renderPlot({
    plot_feature(data = tmp.bcs_merge(), peak_id = input$select_peak,  
                         size = input$size_spot)
  })
  
  output$return_cluster_feature <- renderUI({
    sliderInput('cluster_feature', 'Select clusrer',
                value = 1, min = 0, 
                max = data_cluster() %>% levels() %>% 
                  as.numeric() %>% max())
  })
  
  # plot for interest feature for interest cluster
  output$spatial_feature_plot_cluster <- renderPlot({
    plot_feature_cluster(data_cluster = tmp.bcs_merge(),
                         peak_id = input$select_peak,
                         size = input$size_spot,
                         interest_cluster = input$cluster_feature)
  })
  
  
  output$umap_plot <- renderPlot({
    umap.p1 <- DimPlot(data_cluster(), reduction = "umap", group.by = "sample")
    umap.p2 <- DimPlot(data_cluster(), reduction = "umap", label = TRUE, repel = TRUE)
    umap.p1 + umap.p2
  })
  
  output$umap_sample <- renderPlot({
    DimPlot(data_cluster(), reduction = "umap", split.by = "sample")
  })
  
  
  output$tsne_sample <- renderPlot({
    DimPlot(data_tsne(), reduction = "tsne", split.by = "sample")
  })
  
}

shinyApp(ui = ui, server = server)

