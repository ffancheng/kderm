library(shiny)
library(shinydashboard)
library(plotly)
library(maniTools)
library(tidyverse)
library(dimRed)
# library(reticulate)
# library(here)
library(viridis)
library(hdrcde)
library(igraph)
library(matrixcalc)
# library(akima)
# library(car)
library(ggforce)
# library(ks)
library(patchwork)
library(copula)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(123)


# dr_demo <- function(sim_data, algor, k, d, kernel = 'rbfdot') {
#   if (is.null(algor) | is.null(sim_data) | is.null(k) | is.null(d)) return(NULL)
#   if (is.na(algor)) return(NULL)
#   algor <- toupper(algor)
#   if (algor == "PCA") {
#     # PCA (on centered and scaled data)
#     run_time <- system.time({
#       pca_dr <- sim_data$data %>% center_and_standardise() %>% prcomp()
#       proj_data <- sim_data$data %*% pca_dr$rotation[,1:2]
#     })
#   }

#   # MDS
#   if (algor == "MDS")
#     run_time <- system.time({ proj_data <- cmdscale(dist(sim_data$data), k = d) })

#   # Isomap
#   if (algor == "ISOMAP")
#     run_time <- system.time({ proj_data <- RDRToolbox::Isomap(sim_data$data, dims = d, k = k)$dim2 })

#   # LLE
#   if (algor == "LLE")
#     run_time <- system.time({ proj_data <- LLE2(sim_data$data, dim = d, k = k) })

#   if (algor == "DIFFUSIONMAP")
#     run_time <- system.time({ proj_data <- diffusionMap::diffuse(dist(sim_data$data), neigen = d)$X })

#   # t-SNE
#   if (algor == "TSNE")
#     run_time <- system.time({ proj_data <- tsne::tsne(sim_data$data, k = d) })

#   # KernelPCA
#   if (algor == "KPCA")
#     run_time <- system.time({ proj_data <- kernlab::kpca(sim_data$data, kernel = kernel, features = d)@pcv })

#   # SPE
#   if (algor == "SPE")
#     run_time <- system.time({ proj_data <- spe::spe(sim_data$data, edim = d)$x })

#   # Laplacian Eigenmaps
#   if (algor == "LE")
#     run_time <- system.time({ proj_data <- Laplacian_Eigenmaps(sim_data$data, k = k, d = d)$eigenvectors })

#   # HessianLLE
#   if (algor == 'HLLE')
#     run_time <- system.time({ proj_data <- Hessian_LLE(sim_data$data, k = k, d = d)$projection })

#   # LTSA
#   if (algor == 'LTSA')
#     run_time <- system.time({ proj_data <- Local_TSA(sim_data$data, k = k, d = d) })

#   p1 <- plotly_2D(proj_data, sim_data$colors)
#   plot_title <- paste(algor, "embedding. Time taken: ", round(run_time[[1]], 3), "s.", sep = "")
#   p1 <- layout(p1, title = plot_title)
#   list(p1, run_time)
# }




ml_outlier <- function(x, s = 2, k = min(10, nrow(x)), radius = 0, 
                       adjacency = NULL, affinity = NULL,
                       method, 
                       annmethod = c("kdtree", "annoy", "hnsw"),
                       eps = 0, nt = 50, nlinks = 16, ef.construction = 200, ef.search = 10,
                       distance = c("euclidean", "manhattan")[1], diag = FALSE,
                       treetype = c("kd", "bd")[1],
                       searchtype = c("standard", "priority", "radius")[3],
                       perplexity = round(k/3), theta = 0.5, # t-SNE
                       invert.h = TRUE,
                       n.plot = 10, ell_size = 1, n.grid = 10, noutliers = 10,
                       prob = c(1, 50, 99),
                       ...) {
  
  if (is.null(x) | is.null(method)) return(NULL)
  if(is.null(x$colors)) x$colors <- x$data[,3] # use 3rd column as colors
  
  run_time <- system.time(
    metriclearn <- metricML(x$data, s, k, radius,
                            method = method, 
                            invert.h = invert.h, eps = eps, nt = nt, nlinks = nlinks,
                            annmethod = annmethod, distance = distance, 
                            treetype = treetype,
                            searchtype = searchtype)
  )[[1]]
  
  p_emb <- plot_embedding(metriclearn, color = x$colors) +
    labs(x = "", y = "", color = "") +
    plot_ellipse(metriclearn, add = TRUE, n.plot = n.plot,
                 color = blues9[5], fill = blues9[5], alpha = 0.1,
                 ell_size = ell_size) +
    ggtitle(paste0(substring(method, 4), " 2D embedding. Time taken: ", round(run_time, 3), "s.", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))

  fn <- metriclearn$embedding
  Rn <- metriclearn$rmetric
  E1 <- fn[,1]; E2 <- fn[,2]
  # prob <- c(1, 50, 99)
  f_vkde <- vkde2d(x = E1, y = E2, h = Rn, n = n.grid) # estimated densities with variable bandwidth
  fxy <- hdrcde:::interp.2d(f_vkde$x, f_vkde$y, f_vkde$z, x0 = E1, y0 = E2) # linear interpolation
  
  p_vkde <- plot_outlier(x = metriclearn, n.grid = n.grid, prob = prob, noutliers = noutliers, scales = 1, ell_size = ell_size)
  p_hdr <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers)
  p_hdr$p <- p_hdr$p +
    plot_ellipse(metriclearn, n.plot = n.plot, add = TRUE, ell_size = ell_size)
  
  return(list(metriclearn = metriclearn, p_emb = p_emb, p_vkde = p_vkde, p_hdr = p_hdr, run_time))
  
}


t1 <- list(
  family = "Times New Roman",
  size = 16,
  color = "Black")


# Define UI for application that plot 2d metadata, 3d dataset, and 2d embedding
ui <- dashboardPage(

  dashboardHeader(title = "Outlier detection via manifold learning",
                  titleWidth = 350,
                  tags$li(a(href = 'https://github.com/ffancheng/kderm',
                            icon("github"),
                            title = "Github Repo"),
                          class = "dropdown"
                  )
  ),

  dashboardSidebar(width = 150,
                   sidebarMenu(
                     menuItem("Manifold learning", tabName = "manifoldlearning", icon = icon("table")),
                     menuItem("Comparison", tabName = "comparison", icon = icon("chart-bar"))
                   )
  ),

  dashboardBody(
    tabItems(
      tabItem(
        tabName = "manifoldlearning",
        fluidRow(
          column(8,
                 fluidRow(
                   column(6, plotOutput("plot_2d")),
                   column(6, plotlyOutput("plot_3d"))
                 ),
                 br(),
                 fluidRow(
                   column(4, plotlyOutput("plot_embed")),
                   column(4, plotlyOutput("plot_outlier")),
                   column(4, plotlyOutput("plot_vkde"))
                 ) 
          ),
          
          column(4,
                 wellPanel(
                   style = "background-color: #00CC96;",
                   h4("Manifold"),
                   ## Load manifold from files
                   # fluidRow(
                   #   column(4,
                   #          fileInput("file_input", "Load File",
                   #                    accept = c(
                   #                      'text/csv', 'text/comma-separated-values',
                   #                      'text/tab-separated-values', 'text/plain', '.csv', '.tsv')
                   #          ),
                   #          checkboxInput('header', 'Header', TRUE)
                   #   ),
                   #   column(4, radioButtons('sep', 'Separator',
                   #                          c(Comma=',', Semicolon=';', Tab='\t'),
                   #                          ',')
                   #   ),
                   #   column(4, radioButtons('quote', 'Quote',
                   #                          c(None='', 'Double Quote'='"', 'Single Quote'="'"),
                   #                          '"')
                   #   )
                   # ),
                   # # textOutput("file_text"),
                   # hr(),
                   
                   # Generate data from different mappings
                   fluidRow(
                     column(4,
                            selectInput("data_input", "Examples",
                                        choices = c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve"),
                                        # c("Swiss Roll", "Swiss Hole", "Corner Planes",   
                                        #           "Punctured Sphere", "Twin Peaks", "Clusters",
                                        #           "Toroidal Helix", "Gaussian"),
                                        selected = "Swiss Roll")
                     ),
                     column(3, numericInput("num_pts", "No. of Points", 1000)),
                     # column(3, uiOutput("data_par")),
                     column(5, selectInput("meta_input", "Meta data type", 
                                           choices = c("copula", "gaussian"),
                                           selected = "gaussian"))
                     # column(3, submitButton("Load Example"))  # press button when changing any input
                   ),
                   
                 ),
                 
                 wellPanel(
                   style = "background-color: #636EFA;",
                   h4("Manifold learning parameters"),                            
                   fluidRow(
                     column(3, numericInput("d", "Target Dimension d", 2, min = 2, max = 2)),
                     column(3, selectInput("distance", "Distance", choices = c("euclidean", "manhattan"), selected = "euclidean")),
                     column(3, selectInput("searchtype", "Search Type", choices = list("KNN" = "standard", "Radius" = "radius"), selected = "radius")), 
                     column(3, uiOutput("k_radius")),
                     # column(4, numericInput("k", "Nearest Neighbors k", 8, min = 1)),  # add radius search
                     # column(4, selectInput("kernel", "Kernel",                    ## TODO
                     #                       choices = c("rbfdot", #"polydot", "vanilladot", "tanhdot",
                     #                                   "laplacedot", "besseldot", "anovadot"#, "splinedot"
                     #                       )))
                   ),
                   # textOutput("comment_text")
                   h5("Approximate nearest neighbor searching"),                                
                   fluidRow(
                     column(6, selectInput("ann", label = "Method", choices = list("k-d trees" = "kdtree", "Annoy" = "annoy", "HNSW" = "hnsw")), seleted = "k-d trees"),
                     column(6, uiOutput("annpar")),   
                   )
                 ),
                 
                 wellPanel(
                   style = "background-color: #4d94ff;",
                   h4("Manifold learning algorithm"),                                  
                   fluidRow(
                     radioButtons("algor", label = "",
                                  choices = list("ISOMAP" = "annIsomap",
                                                 "Locally Linear Embedding (LLE)" = "annLLE",
                                                 "Laplacian Eigenmaps" = "annLaplacianEigenmaps",
                                                 "Hessian LLE" = "annHLLE", 
                                                 "t-SNE" = "anntSNE",
                                                 "UMAP" = "annUMAP"),
                                  inline = TRUE)
                   ),
                   # textOutput("plot_text")  ## TODO: not showing the number
                 ),
                 
                 wellPanel(
                   style = "background-color: #ff6666;",
                   h4("Embedding and outlier plot parameters"),
                   fluidRow(
                     column(3, numericInput("n.plot", "Number of ellipses", 10, min = 0, step = 1)),
                     column(3, numericInput("ell_size", "Ellipse size", 0, min = 0, max = 100)),
                     column(3, numericInput("noutliers", "Number of outliers", 10, min = 1, max = 1e8, step = 1)),
                     column(3, numericInput("n.grid", "Number of grid points for KDE", 10, min = 1, max = 500, step = 5)),
                   )
                 ),
                 
                 # wellPanel(
                 #   textOutput("metric.learn")
                 # )
          )       
        )
      ),
      
      
      # Second tab
      tabItem(
        tabName = "comparison",
        fluidRow(
          column(12, offset = 0,
          checkboxGroupInput("algor_group", label = h3("Algorithms"),
                             choices = list("ISOMAP" = "annIsomap",
                                            "Locally Linear Embedding (LLE)" = "annLLE",
                                            "Laplacian Eigenmaps" = "annLaplacianEigenmaps",
                                            "Hessian LLE" = "annHLLE", 
                                            "t-SNE" = "anntSNE",
                                            "UMAP" = "annUMAP"),
                             inline = TRUE)
        )),
        fluidRow(
          column(1, actionButton("button", "Update")),
          column(11, textOutput("info_text"))
        ),
        br(),
        
        fluidRow( uiOutput("plot_first_row") ),
        fluidRow( uiOutput("plot_second_row") ),
        fluidRow( uiOutput("plot_third_row") ),
        fluidRow( uiOutput("plot_fourth_row") ),
        fluidRow( uiOutput("plot_fifth_row") ),
        fluidRow( uiOutput("plot_sixth_row") )
      )

    )
  )
)




server <- shinyServer(function(input, output, session) {
  
  # # read data from input file and assign colors
  # data_from_file <- reactive({
  #   inFile <- input$file_input
  #   if (is.null(inFile)) return（NULL）
  #   sim_data <- read.csv( inFile$datapath, header = input$header,
  #                         sep = input$sep, quote = input$quote )
  #   if (ncol(sim_data) >= 4) {
  #     colors <- sim_data[,4]
  #   } else {
  #     colors <- z # sim_data[,3], use 3rd dimension as colors  ## TODO, check z
  #   }
  #   list(data = as.matrix(sim_data), colors = colors)
  # })
  
  # reduce_to_3d <- reactive({
  #   sim_data <- data_from_file()
  #   sim_data$data <- sim_data$data[,1:3]
  #   sim_data
  # })
  
  sim_data <- reactive({ mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input) })

  # DR_data <- reactiveValues(simulation = sim_data())
  # total_time <- reactiveValues(time_taken = NULL)
  metriclearn_data <- reactive({
    ml_outlier(x = sim_data(), s = input$d, 
              k = input$searchtype_parameter,
              radius = input$searchtype_parameter,
              method = input$algor, 
              invert.h = TRUE, eps = input$ann_parameter,
              nt = input$ann_parameter, nlinks = input$ann_parameter,
              annmethod = input$ann, distance = input$distance, 
              treetype = "kd",
              searchtype = input$searchtype,
              n.plot = input$n.plot,
              n.grid = input$n.grid,
              noutliers = input$noutliers, 
              ell_size = input$ell_size)
  })
  
  ## Used for maniTools package to simulate data
  # output$data_par <- renderUI({
  #   data_param_label <- switch(input$data_input,
  #                              "Swiss Roll" = "Height",
  #                              "Swiss Hole" = "Height",
  #                              "Corner Planes" = "Angles",
  #                              "Punctured Sphere" = "Z scale",
  #                              "Twin Peaks" = "Z scale",
  #                              "Clusters" = "Number of clusters",
  #                              "Toroidal Helix" = "Sample rate",
  #                              "Gaussian" = "Sigma")
  #   initial_value <- switch(input$data_input,
  #                           "Swiss Roll" = 1,
  #                           "Swiss Hole" = 1,
  #                           "Corner Planes" = 45,
  #                           "Punctured Sphere" = 1,
  #                           "Twin Peaks" = 1,
  #                           "Clusters" = 3,
  #                           "Toroidal Helix" = 1,
  #                           "Gaussian" = 1)
  #   numericInput("data_parameter", data_param_label, value = initial_value)
  # })

  output$annpar <- renderUI({
    ann_param_label <- switch(input$ann,
                              "kdtree" = "epsilon",
                              "annoy" = "ntrees",
                              "hnsw" = "nlinks")
    ann_initial_value <- switch(input$ann,
                                "kdtree" = 0,
                                "annoy" = 50,
                                "hnsw" = 16)
    numericInput("ann_parameter", label = ann_param_label, value = ann_initial_value)
  })
  
  output$k_radius <- renderUI({
    searchtype_label <- switch(input$searchtype,
                              "standard" = "No. of nearest neighbors K",
                              "radius" = "Search radius")
    searchtype_initial_value <- switch(input$searchtype,
                                       "KNN" = 10,
                                       "radius" = 10)
    numericInput("searchtype_parameter", label = searchtype_label, value = searchtype_initial_value)
  })
  
  # First tab ================================================================================
  output$plot_3d <- renderPlotly({
    # if (!is.null(data_from_file())) {
    #   sim_data <- reduce_to_3d()
    # } else {
      # sim_data <- mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input)
      # DR_data$simulation <- sim_data
    # }
    if (is.null(sim_data()$data) | (ncol(sim_data()$data) < 3)) {
      plotly_empty(type = "scatter", mode = "markers")
    } else {
      # plotly_3D(sim_data()) %>% layout(title = paste("3D", input$data_input, "mapping data"))
      plot_ly( data.frame(sim_data()$data), x = ~ x, y = ~ y, z = ~ z, color = sim_data()$colors, colors = scales::hue_pal()(4),
               type = "scatter3d", mode = "markers", marker = list(size = 4) ) %>% 
        layout( title = list(text = paste("3D", input$data_input, "mapping data"), font = t1) )
    }
  })
  
  output$plot_2d <- renderPlot({
    # if (!is.null(data_from_file())) {
    #   sim_data <- reduce_to_3d()
    # } else {
      # sim_data <- mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input)
      # DR_data$simulation <- sim_data
    # }
    if (is.null(sim_data()$data) | (ncol(sim_data()$data) < 3)) {
      # plotly_empty(type = "scatter", mode = "markers")
      ggplot()
    } else {
      viridis_color <- if(input$meta_input != "gaussian") scale_color_viridis(direction = 1) else NULL # scale_color_manual(values = c('#00CC96', '#EF553B', '#636EFA', '#AB63FA')) # Not matching the default plotly colors
      sim_data()$metadata %>%
        as_tibble() %>% 
        ggplot(aes(x, y, color = sim_data()$colors)) +
        geom_point(alpha = 0.8) +
        viridis_color + # use four gaussian index for colors
        # scale_color_viridis() +
        labs(color = "", title = paste("Meta data for", input$data_input, "mapping"), color = "") +
        theme(plot.title = element_text(hjust = 0.5))
      # plotly_2D(sim_data()$metadata, colors = sim_data()$colors) # if change to plotlyOutput and renderPlotly
    }
  })
  
  output$plot_outlier <- renderPlot({
    # if (!is.null(data_from_file())) {
    #   sim_data <- reduce_to_3d()
    # } else {
    #   sim_data <- mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input)
    #   DR_data$simulation <- sim_data
    # }
    if (is.null(sim_data()$data) | (ncol(sim_data()$data) < 3)) {
      # plotly_empty(type = "scatter", mode = "markers")
      ggplot()
    } else {
    # res <- dr_demo(sim_data = DR_data$simulation, algor = "isomap",    # TODO
    #                k = input$searchtype_parameter, d = input$d, kernel = input$kernel)
    # total_time$time_taken <- res[[2]]
    # res[[1]]
    
    # res <- ml_outlier(x = sim_data(), s = input$d, 
    #                   k = input$searchtype_parameter,
    #                   radius = input$searchtype_parameter,
    #                   method = input$algor, 
    #                   invert.h = TRUE, eps = input$ann_parameter,
    #                   nt = input$ann_parameter, nlinks = input$ann_parameter,
    #                   annmethod = input$ann, distance = input$distance, 
    #                   treetype = "kd",
    #                   searchtype = input$searchtype,
    #                   n.plot = input$n.plot,
    #                   n.grid = input$n.grid,
    #                   noutliers = input$noutliers, 
    #                   ell_size = input$ell_size)
    
    # metriclearn_data$metric <- res$metriclearn
    # total_time$time_taken <- res$run_time # reactiveValues

    # ((res$p_hdr$p + res$p_vkde$p) + coord_fixed()) + plot_layout(guides = 'collect')
    
    }
  })

  output$plot_embed <- renderPlotly({
    ggplotly(metriclearn_data()$p_emb) %>%
      layout(title = list(font = t1) )
  })

  output$plot_outlier <- renderPlotly({
    ggplotly(metriclearn_data()$p_hdr$p) %>%
      layout(title = list(text = "Outliers with fixed bandwidth", font = t1) )
  })

  output$plot_vkde <- renderPlotly({
    ggplotly(metriclearn_data()$p_vkde$p) %>%
      layout(title = list(text = "Outliers with variable bandwidth", font = t1) )
  })
  

  # output$file_text <- renderText({"Plot first 3 dimensions only.
  #    Use the 4th dimension as colors if any; otherwise, the 3rd dimension is used."})
  # output$comment_text <- renderText({"The target dimension is fixed at 2."})
  output$plot_text <- renderPrint({
    cat("Time taken:", metriclearn_data()$run_time, "s. \n")
  })
  
  
  # Second tab ================================================================================
  output$info_text <- renderText({"(Note: some algorithms may take long to run (e.g. Isomap and t-SNE),
     please avoid clicking the 'Update' button while the calculation is being performed.)"})
  # output$c_plot_1 <- renderPlotly({
  #   # if (!is.null(data_from_file())) {
  #   #   sim_data <- reduce_to_3d()
  #   # } else {
  #   #   sim_data <- mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input)
  #   #   DR_data$simulation <- sim_data
  #   # }
  #   plotly_3D(sim_data())
  # })
  
  output$plot_first_row <- renderUI({
    plot_output_list <- lapply(1:2, function(i) {
      column(6, plotlyOutput(paste0("c_plot_", i + 1)))
    })
    do.call(tagList, plot_output_list)
  })
  output$plot_second_row <- renderUI({
    plot_output_list <- lapply(3:5, function(i) {
      column(4, plotlyOutput(paste0("c_plot_", i + 1)))
    })
    do.call(tagList, plot_output_list)
  })
  output$plot_third_row <- renderUI({
    plot_output_list <- lapply(6:8, function(i) {
      column(4, plotlyOutput(paste0("c_plot_", i + 1)))
    })
    do.call(tagList, plot_output_list)
  })
  output$plot_fourth_row <- renderUI({
    plot_output_list <- lapply(9:11, function(i) {
      column(4, plotlyOutput(paste0("c_plot_", i + 1)))
    })
    do.call(tagList, plot_output_list)
  })


  observeEvent(input$button, {
    algor_list <- input$algor_group
    for (i in 1:6) {
      local({
        local_i <- i + 1
        output[[paste0("c_plot_", local_i)]] <-
          renderPlotly({
            if ((local_i - 1) %in% seq_along(algor_list)) {
              # if (!is.null(data_from_file())) {
              #   DR_data$simulation <- data_from_file()
              # }
              # res <- dr_demo(sim_data(), algor = algor_list[local_i - 1],
              #                k = input$k, d = input$d, kernel = input$kernel)
              # total_time$time_taken <- res[[2]]
              # res[[1]]
              res <- ml_outlier(x = sim_data(), s = input$d, 
                                k = input$searchtype_parameter,
                                radius = input$searchtype_parameter,
                                method = algor_list[local_i - 1], #input$algor, 
                                invert.h = TRUE, eps = input$ann_parameter,
                                nt = input$ann_parameter, nlinks = input$ann_parameter,
                                annmethod = input$ann, distance = input$distance, 
                                treetype = "kd",
                                searchtype = input$searchtype,
                                n.plot = input$n.plot,
                                n.grid = input$n.grid,
                                noutliers = input$noutliers, 
                                ell_size = input$ell_size)

            } else {
              plotly_empty(type = "scatter", mode = "markers")
            }
          })
      })
    }
  })
})


# Run the application
shinyApp(ui = ui, server = server)
